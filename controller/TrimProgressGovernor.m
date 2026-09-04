classdef TrimProgressGovernor < handle
    % TRIMPROGRESSGOVERNOR scalar reference governor (option 1) on the
    % trim-schedule reference trajectory, with a PREDICTIVE liveness gate.
    %
    % Classic scalar/dynamic reference governor (Garone, Di Cairano,
    % Kolmanovsky 2017, Eq. (9); Kapasouris 1990, Gilbert et al. 1994-95):
    %
    %   v(t) = v(t-1) + kappa(t) * (r(t) - v(t-1)),   kappa in [0, 1]
    %
    % where kappa is maximized subject to the governed command keeping the
    % state constraint-admissible. Here the "reference" is the PROGRESS
    % along the trim schedule and admissibility is expressed through the
    % LON BRT value functions.
    %
    % PREDICTIVE gate (LivenessPredictor): instead of only asking "am I
    % inside BRT_{k+1} right now", the governor forecasts the value at the
    % handover deadline,
    %
    %   V_hat = V_next + dV/dt * (t_rem - t_mar),
    %
    % and intervenes as soon as V_hat >= 0. Two failure modes need OPPOSITE
    % responses, which is why they are handled separately:
    %
    %   LATE ENTRY (V_next >= 0, closing too slowly)
    %       -> SLOW DOWN, do not freeze. Freezing the command mid-segment
    %          parks it at a trim-line point that is itself outside
    %          BRT_{k+1}, so the vehicle converges there and can never
    %          enter (deadlock). Advancing at kappa* stretches the deadline
    %          while still pulling the command toward trim k+1.
    %   RE-EXIT (V_next < 0 but rising)
    %       -> FREEZE. The vehicle is diverging from the trim line on a
    %          transient; holding the command lets the closed loop settle
    %          back to a trim point, which is the center of the BRT.
    %
    % The present-tense node gate is KEPT as a hard backstop: a handover is
    % never executed unless V_next < -eps_ready at that instant.
    %
    % Units: ft/s, rad/s, rad. Mission frame rules follow
    % ReferenceTrajectory('trim_schedule'): vel(3) = 0, pos(3) = -80 ft.

    properties
        traj          % TrimScheduleTrajectory.build() output (schedule)
        brt           % 1 x n_trim cell of 4D BRT value arrays
        gv            % 1 x 4 cell of BRT grid vectors
        eps_ready     % readiness margin on V_next at the node gate (>= 0)
        dt            % control sample time [s]

        pred           % LivenessPredictor
        use_predictor  % false -> present-tense gate only (baseline behaviour)
        predict_entry  % act on a LATE-ENTRY forecast (see note below)
        predict_reexit % act on a RE-EXIT forecast
        kappa_min      % floor on kappa in late-entry mode (anti-deadlock)

        progress      % governed progress in [1, n_trim]
        pos_n         % governed north position [ft]
        n_hold        % samples with kappa == 0 (diagnostics)
        n_slow        % samples with 0 < kappa < 1 (diagnostics)
        last_seg      % segment of the previous call (anchor-change detect)
        last_info
    end

    methods
        function obj = TrimProgressGovernor(traj, brt, gv, dt, eps_ready, opts)
            if nargin < 5 || isempty(eps_ready), eps_ready = 0.05; end
            if nargin < 6 || isempty(opts), opts = struct(); end
            obj.traj      = traj;
            obj.brt       = brt;
            obj.gv        = gv;
            obj.dt        = dt;
            obj.eps_ready = eps_ready;

            obj.use_predictor  = getfield_default(opts, 'use_predictor', true);
            obj.predict_entry  = getfield_default(opts, 'predict_entry', false);
            obj.predict_reexit = getfield_default(opts, 'predict_reexit', true);
            obj.kappa_min      = getfield_default(opts, 'kappa_min', 0.1);
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
            obj.last_seg  = [];
            obj.last_info = [];
            obj.pred.reset();
        end

        function [ref, info] = step(obj, state, t)
            % One governor update.
            %   state : 12x1 vehicle state [rn re rd u v w phi th psi p q r]
            %   t     : current time [s]
            %   ref   : struct pos/vel/chi/chi_dot for RSLQR.control
            x_lon = state([4, 6, 11, 8]);          % [u; w; q; theta]
            n_tr  = obj.traj.n_trim;

            k  = min(floor(obj.progress), n_tr - 1);   % current segment
            kn = k + 1;
            V_cur  = obj.V_of(k,  x_lon);
            V_next = obj.V_of(kn, x_lon);

            % --- predictive forecast --------------------------------------
            % Anchor change invalidates the finite-difference history.
            if isempty(obj.last_seg) || obj.last_seg ~= k
                obj.pred.reset();
                obj.last_seg = k;
            end
            t_rem = (kn - obj.progress) * obj.traj.T_seg(k);
            pr = obj.pred.update(t, V_next, t_rem);

            % --- kappa selection ------------------------------------------
            in_hold_phase = (t < obj.traj.t_node(1)) || (obj.progress >= n_tr);

            if in_hold_phase
                kappa = 1;  gate = 'hold-phase';
            elseif V_cur >= 0 && V_next >= 0
                kappa = 0;  gate = 'outside-both-hold';
            elseif obj.next_step_crosses_node() && ~(V_next < -obj.eps_ready)
                kappa = 0;  gate = 'node-gate-not-ready';
            elseif obj.use_predictor && ~pr.ready_pred && pr.ready_now ...
                    && obj.predict_reexit
                kappa = 0;  gate = 'predict-reexit-freeze';
            elseif obj.use_predictor && ~pr.ready_pred && ~pr.ready_now ...
                    && obj.predict_entry
                kappa = max(pr.kappa_star, obj.kappa_min);
                gate  = 'predict-late-slow';
            else
                kappa = 1;  gate = 'advance';
            end

            % --- progress update ------------------------------------------
            if t >= obj.traj.t_node(1) && obj.progress < n_tr
                dp = obj.dt / obj.traj.T_seg(k);
                obj.progress = min(obj.progress + kappa * dp, n_tr);
            end
            obj.n_hold = obj.n_hold + (kappa == 0);
            obj.n_slow = obj.n_slow + (kappa > 0 && kappa < 1);

            % --- governed reference sample --------------------------------
            ks   = min(floor(obj.progress), n_tr - 1);
            frac = obj.progress - ks;
            u_ref = obj.traj.trim_lon(1, ks) + ...
                    frac * (obj.traj.trim_lon(1, ks + 1) - obj.traj.trim_lon(1, ks));
            obj.pos_n = obj.pos_n + u_ref * obj.dt;

            ref = struct('pos', [obj.pos_n; 0; -80], ...
                         'vel', [u_ref; 0; 0], ...
                         'chi', 0, 'chi_dot', 0);

            info = struct('kappa', kappa, 'gate', gate, 'progress', obj.progress, ...
                          'seg', k, 'V_cur', V_cur, 'V_next', V_next, ...
                          'u_ref', u_ref, 't_rem', t_rem, ...
                          'V_hat', pr.V_hat, 'dVdt', pr.dVdt, ...
                          'kappa_star', pr.kappa_star);
            obj.last_info = info;
        end
    end

    methods (Access = private)
        function tf = next_step_crosses_node(obj)
            % Would a full-rate step cross the next trim node?
            k  = min(floor(obj.progress), obj.traj.n_trim - 1);
            dp = obj.dt / obj.traj.T_seg(k);
            tf = (obj.progress + dp) >= k + 1;
        end

        function V = V_of(obj, k, x_lon)
            % BRT value of anchor k at absolute lon state (clamped to grid).
            xc = x_lon(:)' - obj.traj.trim_lon(:, k)';
            for d = 1:4
                xc(d) = min(max(xc(d), obj.gv{d}(1)), obj.gv{d}(end));
            end
            V = interpn(obj.gv{1}, obj.gv{2}, obj.gv{3}, obj.gv{4}, ...
                        obj.brt{k}, xc(1), xc(2), xc(3), xc(4), 'linear');
        end
    end
end
