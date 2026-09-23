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
        safety_filter           % LivenessFilter (per-axis liveness), or [] to bypass
        safety_filter_wh_anchor % WH anchor [ft/s] for the BRT scheduling
    end

    methods
        function obj = Controller(config, dt)
            % config : ControllerConfig (provides .rslqr and .filter)
            % dt     : sim timestep [s] (servo-compensator discretization)
            obj.controller_config   = config;
            obj.baseline_controller = RSLQR(config, dt);

            % Liveness filter over filterCfg.axes ({'lon'} or {'lon','lat'}).
            % Loads BRT value functions from FilterConfig.tables_dir_default;
            % passes through when tables are absent or (uh,wh) is outside
            % coverage. UH/WH breakpoints come from the baseline trim table.
            filterCfg = config.filter;
            obj.safety_filter = LivenessFilter(filterCfg.axes, filterCfg.mode, ...
                                               FilterConfig.tables_dir_default, ...
                                               obj.baseline_controller.UH, obj.baseline_controller.WH, ...
                                               filterCfg.use_disturbance);

            obj.safety_filter_wh_anchor = filterCfg.wh_anchor;
            if isempty(obj.safety_filter_wh_anchor)
                obj.safety_filter_wh_anchor = obj.baseline_controller.WH(3); % default anchor
            end
        end

        function [engine, surface] = control(obj, state, ref)
            % Full closed-loop control. Returns absolute effector commands.
            %
            % Safety-filter modes: 'off' -> baseline only; 'blend' -> CBF-like
            % liveness projection over the filter's active axes (one constraint
            % per axis, solved jointly over the shared 13-effector input).
            [perturb_cmd, U0] = obj.baseline_controller.control(state, ref);

            if isempty(obj.safety_filter)
                mode = 'off';
            else
                mode = lower(obj.safety_filter.mode);
            end

            if ~strcmp(mode, 'off')
                % Build current/next anchor frames for each active axis and
                % project the nominal 13-effector perturbation onto the tube.
                axes = obj.safety_filter.axes;
                frames = struct();
                for i = 1:numel(axes)
                    ax = axes{i};
                    [frames.(ax).curr, frames.(ax).next] = obj.build_frames(ax, state);
                end

                [u_all_f, info] = obj.safety_filter.filter(frames, perturb_cmd(1:13), U0);
                perturb_cmd(1:13) = u_all_f - info.target_dtrim;
            end

            [engine, surface] = obj.baseline_controller.total_cmd(perturb_cmd, U0);
        end

        function [curr, next] = build_frames(obj, ax, state)
            % Current/next BRT anchor frames for axis ax. The anchor UH
            % breakpoints bracket the current body-x velocity; WH is fixed to
            % safety_filter_wh_anchor.
            UH = obj.baseline_controller.UH;
            uhA = max(UH(1), min(UH(end), state(4)));
            [~, current_id] = min(abs(UH - uhA));
            next_id = min(length(UH), current_id + 1);

            curr = obj.make_frame(ax, state, current_id, [ax '_current']);
            next = obj.make_frame(ax, state, next_id,    [ax '_next']);
        end

        function info = make_frame(obj, ax, state, uh_id, name)
            % One BRT anchor frame for axis ax at trim-table UH index uh_id and
            % the fixed WH anchor. The perturbation state and reduced dynamics
            % are selected per axis from FilterConfig.axisSpec (state_rows /
            % trim_rows) and the RSLQR LON/LAT reduced model.
            spec     = FilterConfig.axisSpec(ax);
            baseline = obj.baseline_controller;

            uhA = baseline.UH(uh_id);
            whA = obj.safety_filter_wh_anchor;

            [X0a, U0a] = baseline.interp_xu0(uhA, whA);
            dyn  = baseline.(upper(ax));                        % LON or LAT reduced model
            Ap_a = baseline.interp_mtrx(dyn.Ap, uhA, whA);
            Bp_a = baseline.interp_mtrx(dyn.Bp, uhA, whA);

            x_anchor = state(spec.state_rows) - X0a(spec.trim_rows);

            info = struct( ...
                'uhA',      uhA, ...
                'whA',      whA, ...
                'X0a',      X0a, ...
                'U0a',      U0a, ...
                'Ap_a',     Ap_a, ...
                'Bp_a',     Bp_a, ...
                'x_anchor', x_anchor(:), ...
                'idx',      uh_id, ...
                'name',     name, ...
                'axis',     ax);
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
