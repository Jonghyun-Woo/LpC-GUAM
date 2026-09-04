classdef TrimPointGovernor < handle
    % TRIMPOINTGOVERNOR trajectory time-scaling reference governor (option 2)
    % on the trim-schedule reference trajectory, with a PREDICTIVE liveness
    % gate.
    %
    % Modeled on the dual-stage machine reference governor of Di Cairano,
    % Goldsmith et al. 2015 (Garone survey Sec. 5.2): the reference PATH is a
    % fixed sequence of points {r_i} that must not be modified; each sampling
    % period the governor only chooses HOW MANY points to advance,
    %
    %   n(t) = max n in {0, 0.5, 1, ..., n_max}
    %   s.t.  reference-dependent constraint + admissibility hold,
    %
    % which time-scales the ideal trajectory.
    %
    % Constraints (our liveness adaptation):
    %   proximity : |u_ref(p_cand) - u_vehicle| <= s_f
    %               (the dual-stage ||p_s - v|| <= s_f constraint: the
    %               commanded point may never run away from the vehicle)
    %   predictive: at the candidate advance rate n, the time to the node is
    %               t_rem/n, and the forecast value at that deadline
    %                   V_next + dV/dt * (t_rem/n - t_mar)
    %               must be negative. Slower n = later deadline = easier.
    %               This is the same LivenessPredictor used by option 1.
    %   node gate : crossing trim node k+1 requires V_{k+1}(x) < -eps_ready
    %               at that instant (hard backstop on the certificate)
    %   pace cap  : governed progress may never run AHEAD of the nominal
    %               schedule position p_nom(t)
    %   hold rule : outside BOTH scheduled BRTs, or a predicted RE-EXIT
    %               -> n = 0 (recursive feasibility: the held command is a
    %               trim-line point the closed loop converges back to)
    %
    % Anti-deadlock: a predicted LATE ENTRY never drives n to 0. Freezing
    % mid-segment parks the command at a trim-line point that is itself
    % outside BRT_{k+1}, so the vehicle would converge there and never
    % enter. Late entry therefore floors the advance at n_min.
    %
    % Units: ft/s, rad/s, rad. Mission frame rules follow
    % ReferenceTrajectory('trim_schedule'): vel(3) = 0, pos(3) = -80 ft.

    properties
        traj          % TrimScheduleTrajectory.build() output (schedule)
        brt           % 1 x n_trim cell of 4D BRT value arrays
        gv            % 1 x 4 cell of BRT grid vectors
        dt            % control sample time [s]
        s_f           % proximity cap on |u_ref - u_vehicle| [ft/s]
        n_max         % max points advanced per period (1 = nominal rate)
        n_min         % floor on the advance under a late-entry forecast
        eps_ready     % readiness margin on V_next at the node gate (>= 0)

        pred           % LivenessPredictor
        use_predictor  % false -> present-tense gate only
        predict_entry  % act on a LATE-ENTRY forecast (default OFF, see note)
        predict_reexit % act on a RE-EXIT forecast

        progress      % governed progress in [1, n_trim]
        pos_n         % governed north position [ft]
        n_hold        % samples with n = 0 (diagnostics)
        n_slow        % samples with 0 < n < 1 (diagnostics)
        n_catchup     % samples with n > 1 (diagnostics)
        last_seg
        last_info
    end

    methods
        function obj = TrimPointGovernor(traj, brt, gv, dt, s_f, n_max, eps_ready, opts)
            % Defaults from closed-loop tuning (2026-07-27):
            %  - s_f = 6 clears the ~4 ft/s natural tracking lag during the
            %    nominal ramp (s_f = 4 throttles constantly).
            %  - n_max = 1: catch-up DISABLED. Sustained >1x-nominal
            %    acceleration demands pitch the vehicle down and bleed
            %    altitude, which the (u,w,q,theta) BRT cannot see, and the
            %    mission ends in a dive (verified with n_max = 1.5 and 3).
            if nargin < 5 || isempty(s_f),       s_f = 6.0;      end
            if nargin < 6 || isempty(n_max),     n_max = 1.0;    end
            if nargin < 7 || isempty(eps_ready), eps_ready = 0.05; end
            if nargin < 8 || isempty(opts),      opts = struct(); end
            obj.traj = traj;  obj.brt = brt;  obj.gv = gv;  obj.dt = dt;
            obj.s_f = s_f;    obj.n_max = n_max;  obj.eps_ready = eps_ready;

            % predict_entry defaults OFF: acting on a late-entry forecast
            % means slowing the command, but the vehicle only enters
            % BRT_{k+1} because the command marches toward trim k+1, so
            % slowing it also delays the entry. Measured on the closed loop
            % (tests/compare_governor_gates.m): +2.5 s mission time for
            % almost no gain in certified time. Re-exit prediction is the
            % half that pays off.
            obj.n_min          = getfield_default(opts, 'n_min', 0.5);
            obj.use_predictor  = getfield_default(opts, 'use_predictor', true);
            obj.predict_entry  = getfield_default(opts, 'predict_entry', false);
            obj.predict_reexit = getfield_default(opts, 'predict_reexit', true);
            obj.pred = LivenessPredictor( ...
                getfield_default(opts, 'fd_window', 0.2), ...
                getfield_default(opts, 't_margin',  0.3), ...
                getfield_default(opts, 'eps_v',     0.0));
            obj.reset();
        end

        function reset(obj)
            obj.progress  = 1;
            obj.pos_n     = 0;
            obj.n_hold    = 0;
            obj.n_slow    = 0;
            obj.n_catchup = 0;
            obj.last_seg  = [];
            obj.last_info = [];
            obj.pred.reset();
        end

        function [ref, info] = step(obj, state, t)
            x_lon = state([4, 6, 11, 8]);
            n_tr  = obj.traj.n_trim;

            k  = min(floor(obj.progress), n_tr - 1);
            V_cur  = obj.V_of(k,     x_lon);
            V_next = obj.V_of(k + 1, x_lon);

            % --- predictive forecast --------------------------------------
            if isempty(obj.last_seg) || obj.last_seg ~= k
                obj.pred.reset();
                obj.last_seg = k;
            end
            t_rem = (k + 1 - obj.progress) * obj.traj.T_seg(k);
            pr = obj.pred.update(t, V_next, t_rem);

            % --- choose the largest admissible advance n ------------------
            n_sel = 0;  gate = 'no-admissible-advance';
            dp    = obj.dt / obj.traj.T_seg(k);
            if t < obj.traj.t_node(1) || obj.progress >= n_tr
                n_sel = 0;  gate = 'hold-phase';
            elseif V_cur >= 0 && V_next >= 0
                n_sel = 0;  gate = 'outside-both-hold';
            elseif obj.use_predictor && obj.predict_reexit ...
                    && pr.ready_now && ~pr.ready_pred
                n_sel = 0;  gate = 'predict-reexit-freeze';
            else
                % nominal schedule position at this wall-clock time
                p_nom = interp1(obj.traj.time, obj.traj.progress, ...
                                min(t, obj.traj.time(end)));
                for n = obj.n_max:-0.5:0
                    p_cand = min(obj.progress + n * dp, n_tr);
                    if ~obj.admissible(p_cand, p_nom, state, V_next), continue; end
                    % predictive: does this advance rate meet the deadline?
                    if obj.use_predictor && obj.predict_entry ...
                            && n > 0 && V_next >= 0
                        t_eff = max(t_rem / n - obj.pred.t_margin, 0);
                        if V_next + pr.dVdt * t_eff >= -obj.pred.eps_v
                            continue;
                        end
                    end
                    n_sel = n;
                    if n == 0,      gate = 'proximity-hold';
                    elseif n < 1,   gate = 'slow';
                    elseif n == 1,  gate = 'nominal';
                    else,           gate = 'catch-up';
                    end
                    break;
                end

                % Anti-deadlock floor: a late-entry forecast must never stop
                % the command completely (see class header).
                if n_sel == 0 && V_next >= 0
                    p_min = min(obj.progress + obj.n_min * dp, n_tr);
                    if obj.admissible(p_min, p_nom, state, V_next)
                        n_sel = obj.n_min;  gate = 'predict-late-slow';
                    end
                end
                obj.progress = min(obj.progress + n_sel * dp, n_tr);
            end

            obj.n_hold    = obj.n_hold    + (n_sel == 0);
            obj.n_slow    = obj.n_slow    + (n_sel > 0 && n_sel < 1);
            obj.n_catchup = obj.n_catchup + (n_sel > 1);

            % --- governed reference sample --------------------------------
            u_ref = obj.u_of(obj.progress);
            obj.pos_n = obj.pos_n + u_ref * obj.dt;
            ref = struct('pos', [obj.pos_n; 0; -80], ...
                         'vel', [u_ref; 0; 0], ...
                         'chi', 0, 'chi_dot', 0);

            info = struct('n', n_sel, 'kappa', min(n_sel, 1), 'gate', gate, ...
                          'progress', obj.progress, 'seg', k, ...
                          'V_cur', V_cur, 'V_next', V_next, 'u_ref', u_ref, ...
                          't_rem', t_rem, 'V_hat', pr.V_hat, 'dVdt', pr.dVdt);
            obj.last_info = info;
        end
    end

    methods (Access = private)
        function tf = admissible(obj, p_cand, p_nom, state, V_next)
            % pace cap: never ahead of the nominal T_seg clock
            if p_cand > p_nom + 1e-12
                tf = false; return;
            end
            % node gate: may not cross an uncertified node (+eps so a
            % progress sitting exactly ON a node counts as already past it)
            if floor(p_cand + 1e-12) > floor(obj.progress + 1e-12) ...
                    && ~(V_next < -obj.eps_ready)
                tf = false; return;
            end
            % proximity: command may not run ahead of the vehicle
            if obj.u_of(p_cand) - state(4) > obj.s_f
                tf = false; return;
            end
            tf = true;
        end

        function u = u_of(obj, p)
            ks = min(floor(p), obj.traj.n_trim - 1);
            fr = p - ks;
            u  = obj.traj.trim_lon(1, ks) + ...
                 fr * (obj.traj.trim_lon(1, ks + 1) - obj.traj.trim_lon(1, ks));
        end

        function V = V_of(obj, k, x_lon)
            xc = x_lon(:)' - obj.traj.trim_lon(:, k)';
            for d = 1:4
                xc(d) = min(max(xc(d), obj.gv{d}(1)), obj.gv{d}(end));
            end
            V = interpn(obj.gv{1}, obj.gv{2}, obj.gv{3}, obj.gv{4}, ...
                        obj.brt{k}, xc(1), xc(2), xc(3), xc(4), 'linear');
        end
    end
end
