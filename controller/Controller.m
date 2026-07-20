classdef Controller < handle
    % Vehicle controller: owns the baseline RSLQR and the safety (liveness)
    % filter, and finishes the full control computation. control() returns the
    % absolute effector commands (engine, surface) that LpC_GUAM hands to the
    % plant.
    %
    % Per-step pipeline:
    %   1. baseline_controller.control   -> nominal effector perturbation + trim U0
    %   2. build the current/next longitudinal BRT anchor frames
    %   3. safety_filter.filter          -> liveness-projected perturbation
    %   4. baseline_controller.total_cmd -> absolute engine/surface commands

    properties
        controller_config       % ControllerConfig (hub-owned)
        baseline_controller     % RSLQR (gain-scheduled control + allocation)
        safety_filter           % LivenessFilter (longitudinal liveness), or [] to bypass
        safety_filter_wh_anchor % WH anchor [ft/s] for the LON BRT scheduling
    end

    methods
        function obj = Controller(config, dt)
            % config : ControllerConfig (provides .rslqr and .filter)
            % dt     : sim timestep [s] (servo-compensator discretization)
            obj.controller_config   = config;
            obj.baseline_controller = RSLQR(config.rslqr, dt);

            % Longitudinal liveness filter (production axis). Loads BRT value
            % functions from FilterConfig.tables_dir_default; passes through
            % when tables are absent or (uh,wh) is outside coverage. Its UH/WH
            % breakpoints come from the baseline controller's trim table.
            filterCfg = config.filter;
            obj.safety_filter = LivenessFilter('lon', filterCfg.mode, ...
                FilterConfig.tables_dir_default, ...
                obj.baseline_controller.UH, obj.baseline_controller.WH);

            obj.safety_filter_wh_anchor = filterCfg.wh_anchor;
            if isempty(obj.safety_filter_wh_anchor)
                obj.safety_filter_wh_anchor = obj.baseline_controller.WH(3); % default anchor
            end
        end

        function [engine, surface] = control(obj, state, ref)
            % Full closed-loop control. Returns absolute effector commands.
            %
            % The safety-filter mode selects what runs after the baseline
            % allocation:
            %   'off'   -> baseline controller only (no liveness filter)
            %   'blend' -> CBF-like blending filter after the baseline
            %   'lr'    -> least-restrictive filter after the baseline
            % 'blend'/'lr' share the same post-baseline projection; the filter
            % law itself is realized inside LivenessFilter.filter per its mode.
            [perturb_cmd, U0] = obj.baseline_controller.control(state, ref);

            if isempty(obj.safety_filter)
                mode = 'off';
            else
                mode = lower(obj.safety_filter.mode);
            end

            switch mode
                case 'off'
                    % Baseline only: nothing added after the nominal allocation.
                case {'blend', 'lr'}
                    [c_brt_info, n_brt_info] = obj.build_brt_frames(state, U0);

                    % Map the nominal 11-effector perturbation into each anchor frame.
                    u_lon_nom = perturb_cmd(1:11);
                    c_brt_info.u_anchor_nom = u_lon_nom + c_brt_info.dtrim;
                    n_brt_info.u_anchor_nom = u_lon_nom + n_brt_info.dtrim;

                    [u_anchor_f, info] = obj.safety_filter.filter(c_brt_info, n_brt_info, state);
                    perturb_cmd(1:11)  = u_anchor_f - info.target_dtrim;
                otherwise
                    error('Controller:mode', ...
                        'Unknown safety-filter mode "%s" (expected off, blend, or lr).', mode);
            end

            [engine, surface] = obj.baseline_controller.total_cmd(perturb_cmd, U0);
        end

        function [c_brt_info, n_brt_info] = build_brt_frames(obj, state, U0)
            % Build the current/next longitudinal BRT anchor frames the liveness
            % filter selects between. The anchor UH breakpoints bracket the
            % current body-x velocity; WH is fixed to safety_filter_wh_anchor.
            UH = obj.baseline_controller.UH;

            uhA = max(UH(1), min(UH(end), state(4)));
            [~, current_id] = min(abs(UH - uhA));
            next_id = min(length(UH), current_id + 1);

            c_brt_info = obj.make_brt_frame(state, U0, current_id, 'current');
            n_brt_info = obj.make_brt_frame(state, U0, next_id,    'next');
        end

        function info = make_brt_frame(obj, state, U0, uh_id, name)
            % Assemble one BRT anchor frame at trim-table UH index uh_id and the
            % fixed WH anchor. u_anchor_nom is filled in by control() once the
            % nominal allocation is known.
            bc    = obj.baseline_controller;
            u_idx = obj.safety_filter.spec.U0_idx;

            uhA = bc.UH(uh_id);
            whA = obj.safety_filter_wh_anchor;

            [X0a, U0a] = bc.interp_xu0(uhA, whA);         % anchor trim state/input
            Ap_a = bc.interp_mtrx(bc.LON.Ap, uhA, whA);   % linear dynamics at anchor trim
            Bp_a = bc.interp_mtrx(bc.LON.Bp, uhA, whA);

            % Perturbation state relative to the anchor trim [u; w; q; theta].
            x_anchor = [state(4)  - X0a(1); ...
                        state(6)  - X0a(3); ...
                        state(11) - X0a(5); ...
                        state(8)  - X0a(11)];
            dtrim = U0(u_idx) - U0a(u_idx);               % 11x1 effector trim offset

            info = struct( ...
                'uhA',          uhA, ...
                'whA',          whA, ...
                'X0a',          X0a, ...
                'U0a',          U0a, ...
                'Ap_a',         Ap_a, ...
                'Bp_a',         Bp_a, ...
                'x_anchor',     x_anchor, ...
                'dtrim',        dtrim, ...
                'u_anchor_nom', [], ...
                'idx',          uh_id, ...
                'name',         name);
        end

        function [x0, engine0, surface0] = initial_condition(obj, refTraj)
            % Trim initial condition for the plant at the first reference point
            % (as GUAM's setupTrim does for the initial reference velocity).
            % Returns the 12x1 rigid state and the actuator trims the driver
            % injects into LpC_GUAM.reset.
            uh0 = refTraj.vel(1, 1);
            wh0 = refTraj.vel(3, 1);
            [X0, U0] = obj.baseline_controller.interp_xu0(uh0, wh0);

            x0        = zeros(12, 1);
            x0(1:3)   = refTraj.pos(:, 1);   % NED position
            x0(4:6)   = X0(1:3);             % body velocity at trim
            x0(7:9)   = X0(10:12);           % Euler angles at trim
            x0(10:12) = X0(4:6);             % body rates at trim

            % Zero effector perturbation -> total_cmd yields the trim actuator
            % positions in plant form ([lift1..8; pusher], [LA RA LE RE RUD]).
            [engine0, surface0] = obj.baseline_controller.total_cmd(zeros(13, 1), U0);
        end

        function [X0, U0] = interp_xu0(obj, uh, wh)
            % Delegate trim lookup to the baseline controller.
            [X0, U0] = obj.baseline_controller.interp_xu0(uh, wh);
        end

        function reset(obj)
            obj.baseline_controller.reset();
            if ~isempty(obj.safety_filter)
                obj.safety_filter.reset_counters();
            end
        end
    end
end
