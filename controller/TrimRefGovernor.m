classdef TrimRefGovernor < handle
    % TRIMREFGOVERNOR scalar reference governor for the trim-schedule
    % transition, following Gilbert & Kolmanovsky, "Nonlinear tracking
    % control in the presence of state and control constraints: a
    % generalized reference governor", Automatica 38 (2002) 2063-2073.
    %
    % Signals (all commands are a forward-speed value in ft/s):
    %   r(t)  scheduled command   - what the T_seg plan asks for
    %   v(t)  applied command     - what RSLQR actually receives
    %   x(t)  vehicle state
    %
    % Governor law, Eq. (5):        v(t) = v(t-1) + kappa(r(t) - v(t-1))
    % Scalar selection, Eq. (10):   kappa = max{k in [0,1] : V(x, v(t-1)+k(r-v(t-1))) <= 0}
    % Small-increment rules:        Eqs. (13)-(15), via delta and eps
    %
    % The value function V is built as in Sec. 4 / Eq. (23) of that paper:
    % hold the candidate command, roll the ACTUAL closed loop forward T
    % seconds, and take the worst constraint value along the way. This is
    % what makes assumption (A5) - "V(x,v) <= 0 implies holding v keeps the
    % constraints satisfied" - true by construction. Using the raw BRT as V
    % instead does NOT satisfy (A5): the BRT certifies reachability under the
    % OPTIMAL control, while the loop actually flies RSLQR, so states near the
    % BRT boundary can be left behind even with the command frozen (measured:
    % 5.3 s outside both BRTs after a freeze).
    %
    % Only the BRT term is used as a constraint for now; more (altitude,
    % attitude, actuator saturation) drop straight into worst_violation().

    properties
        traj        % TrimScheduleTrajectory.build() output (the schedule)
        brtV        % BRTValue - continuous V_BRT(x, v)
        plant       % LpC_GUAM used for the rollout (snapshot/restore)
        lin         % ClosedLoopModel or [] - when set, predict() propagates
                    % the linearised closed loop instead of stepping the
                    % plant. Same Eq. (23) construction, same worst-case over
                    % [0,T); only how the trajectory is obtained changes.
        dt

        T           % prediction horizon [s]
        N           % rollout steps = T/dt
        delta       % min meaningful command move [ft/s]   Eq. (13)
        eps_margin  % depth required for a small move       Eq. (14)
        margin      % constraint tightening: admissible means V <= -margin,
                    % not V <= 0. Needed when predict() uses a model that is
                    % not the real loop (lin ~= []): the prediction error has
                    % to be paid for somewhere, and tightening the constraint
                    % by a bound on that error restores (A5). 0 for the exact
                    % nonlinear rollout.
        M           % kappa grid resolution; see the default below

        v           % currently applied command [ft/s]
        pos_n       % governed north position reference [ft]
        n_call      % diagnostics
        n_rollout
        last_info
    end

    methods
        function obj = TrimRefGovernor(traj, brtV, plant, dt, opts)
            if nargin < 5 || isempty(opts), opts = struct(); end
            obj.traj  = traj;
            obj.brtV  = brtV;
            obj.plant = plant;
            obj.dt    = dt;

            obj.T          = getfield_default(opts, 'T',     6.0);
            obj.delta      = getfield_default(opts, 'delta', 0.3);
            obj.eps_margin = getfield_default(opts, 'eps',   0.02);
            % kappa is a FRACTION of (r - v), so the smallest step the grid can
            % take is |r - v|/M ft/s. An impatient schedule makes |r - v| large,
            % and with M = 10 that smallest step is already too big to certify:
            % kappa collapses to 0 and the command never advances. M = 40 clears
            % it (T_seg = 0.5 s stalls at 8.5 ft/s with M = 10, completes in
            % 51.6 s with M = 40). Eq. (10) assumes the exact kappa*, so this is
            % a discretisation limit, not a property of the governor.
            obj.M          = getfield_default(opts, 'M',     40);
            obj.lin        = getfield_default(opts, 'lin',   []);
            obj.margin     = getfield_default(opts, 'margin', 0);
            obj.N          = round(obj.T / dt);
            obj.reset();
        end

        function reset(obj)
            obj.v         = obj.traj.trim_lon(1, 1);   % first trim speed
            obj.pos_n     = 0;
            obj.n_call    = 0;
            obj.n_rollout = 0;
            obj.last_info = [];
        end

        % ---------------------------------------------------------------
        % Algorithm 1 - constraint-violation prediction, Eq. (23)
        % ---------------------------------------------------------------
        function V = predict(obj, v_cand, early_exit)
            % V = worst constraint violation over [0, T) when the command is
            % HELD at v_cand. V <= 0 means the candidate is admissible.
            %
            % early_exit: when only the sign is needed, stop as soon as V > 0.
            % The value itself is only required for the Eq. (14) test.
            if nargin < 3, early_exit = true; end
            if ~isempty(obj.lin)
                V = obj.predict_linear(v_cand, early_exit);
                return;
            end
            obj.n_rollout = obj.n_rollout + 1;

            snap = obj.plant.saveState();      % rewind point
            p    = obj.pos_n;
            V    = -inf;

            for i = 1:obj.N
                % constraint value at the CURRENT step, before advancing
                V = max(V, obj.worst_violation(v_cand));
                if early_exit && V > -obj.margin
                    break;
                end
                % one simulated step with the command held at v_cand
                p   = p + v_cand * obj.dt;
                ref = struct('pos', [p; 0; -80], 'vel', [v_cand; 0; 0], ...
                             'chi', 0, 'chi_dot', 0);
                obj.plant.step(ref);
            end

            obj.plant.restoreState(snap);      % the real vehicle never moved
        end

        function V = predict_linear(obj, v_cand, early_exit)
            % Same question as predict(), answered by propagating
            %     d(k+1) = c(v) + Phi(v) d(k),    d = xi - xi_e(v)
            % instead of stepping the plant. The plant is only read, never
            % advanced, so no snapshot is needed.
            obj.n_rollout = obj.n_rollout + 1;
            xi = obj.lin.pack(obj.plant, [obj.pos_n; 0; obj.lin.Z_REF]);
            [P, xe, cc] = obj.lin.at(v_cand);
            d = xi - xe;
            sel = [4 6 11 8];                  % u, w, q, theta within xi
            V = -inf;
            for i = 1:obj.N
                V = max(V, obj.brtV.value(xe(sel) + d(sel), v_cand));
                if early_exit && V > -obj.margin
                    break;
                end
                d = cc + P*d;
            end
        end

        % ---------------------------------------------------------------
        % Algorithm 2 - scalar selection, Eqs. (10), (12)-(15)
        % ---------------------------------------------------------------
        function [ref, info] = step(obj, t)
            obj.n_call = obj.n_call + 1;
            r = obj.schedule(t);
            v_prev = obj.v;

            gate = 'advance';  kappa = 1;
            V_sched = NaN;  V_hold = NaN;
            cands = zeros(2, 0);        % [kappa_try; V_try] actually evaluated

            if abs(r - v_prev) < 1e-12
                kappa = 1;  gate = 'on-schedule';
            else
                % (1st) is the schedule itself admissible?
                thr = -obj.margin;
                V_sched = obj.predict(r);
                cands(:, end+1) = [1; V_sched];
                if V_sched <= thr
                    kappa = 1;  gate = 'advance';
                else
                    % (2nd) is holding admissible?  Eq. (12)
                    V_hold = obj.predict(v_prev, false);
                    cands(:, end+1) = [0; V_hold];
                    if V_hold > thr
                        kappa = 0;  gate = 'hold-infeasible';
                    else
                        % (3rd) line search, Eq. (10)
                        kappa = 0;  gate = 'hold';
                        for j = obj.M-1 : -1 : 1
                            k_try = j / obj.M;
                            V_try = obj.predict(v_prev + k_try*(r - v_prev));
                            cands(:, end+1) = [k_try; V_try]; %#ok<AGROW>
                            if V_try <= thr
                                kappa = k_try;  gate = 'partial';
                                break;
                            end
                        end
                        % (4th) small-increment suppression, Eqs. (13)-(15)
                        if kappa * abs(r - v_prev) < obj.delta && ...
                                V_hold > -obj.eps_margin
                            kappa = 0;  gate = 'suppressed';
                        end
                    end
                end
            end

            obj.v     = v_prev + kappa * (r - v_prev);
            obj.pos_n = obj.pos_n + obj.v * obj.dt;   % integrate v, not r

            ref = struct('pos', [obj.pos_n; 0; -80], 'vel', [obj.v; 0; 0], ...
                         'chi', 0, 'chi_dot', 0);
            info = struct('kappa', kappa, 'gate', gate, 'v', obj.v, 'r', r, ...
                          'v_prev', v_prev, 'V_sched', V_sched, ...
                          'V_hold', V_hold, 'cands', cands, ...
                          'n_rollout', obj.n_rollout);
            obj.last_info = info;
        end

        function ref = hold(obj)
            % Re-issue the current command without re-solving. The governor
            % is decimated (a command does not need a 100 Hz update, and each
            % solve costs several rollouts), but the position reference must
            % still integrate at the simulation rate.
            obj.pos_n = obj.pos_n + obj.v * obj.dt;
            ref = struct('pos', [obj.pos_n; 0; -80], 'vel', [obj.v; 0; 0], ...
                         'chi', 0, 'chi_dot', 0);
        end

        function r = schedule(obj, t)
            % r(t): the forward speed the nominal T_seg plan asks for.
            r = interp1(obj.traj.time, obj.traj.lon(1, :), ...
                        min(max(t, obj.traj.time(1)), obj.traj.time(end)));
        end
    end

    methods (Access = private)
        function h = worst_violation(obj, v_cand)
            % Constraint value at the plant's current state, for command
            % v_cand. Negative = margin, positive = violation.
            % Extra constraints go here as extra max() terms.
            x_lon = obj.plant.state([4 6 11 8]);
            h = obj.brtV.value(x_lon, v_cand);          % h1: BRT departure
        end
    end
end
