classdef LivenessFilter < handle
    % HJ-reachability liveness filter behind a nominal controller.
    %
    % Each active axis (obj.axes) adds one CBF-like liveness constraint on the
    % shared 13-effector perturbation u; blend mode projects the nominal command:
    %   min ||u - u_nom||^2
    %   s.t. beta_a'*u <= -gamma*(V_a + live_margin) - alpha_a - s_a,  lb<=u<=ub
    % with alpha_a = gradV_a'*(Ap_a*x_a), beta_a = (gradV_a'*Bp_a)'. s_a is the
    % optional worst-case model-mismatch term (use_disturbance). Units: ft/s, rad/s, rad.

    properties
        mode
        gamma
        live_margin

        axes            % active axes, e.g. {'lon'} or {'lon','lat'}
        value_function  % struct keyed by axis; each carries its .axis_spec
        combined_spec   % shared 13-effector input spec

        use_disturbance % include the worst-case mismatch term s_a
        dist            % struct keyed by axis: .beta (15x4xUH), .half (4x1)

        n_calls
        n_active
        last_info
    end

    methods
        function obj = LivenessFilter(axes, mode, tables_dir, uh_breakpoint, wh_breakpoint, use_disturbance)
            if ischar(axes) || isstring(axes), axes = cellstr(axes); end
            obj.axes        = lower(axes(:))';
            obj.mode        = lower(mode);
            obj.gamma       = FilterConfig.gamma;
            obj.live_margin = FilterConfig.live_margin;

            if nargin < 3 || isempty(tables_dir)
                tables_dir = FilterConfig.tables_dir_default;
            end
            if nargin < 6 || isempty(use_disturbance), use_disturbance = false; end

            obj.value_function = struct();
            for i = 1:numel(obj.axes)
                axis = obj.axes{i};
                obj.value_function.(axis) = ValueFunction(FilterConfig.axisSpec(axis), ...
                                                          tables_dir, uh_breakpoint, wh_breakpoint);
            end
            obj.combined_spec = FilterConfig.combinedSpec();

            obj.use_disturbance = logical(use_disturbance);
            obj.dist = struct();
            if obj.use_disturbance
                obj.load_disturbance();
            end

            obj.n_calls  = 0;
            obj.n_active = 0;
            obj.last_info = [];
        end

        function load_disturbance(obj)
            % Load the per-axis model-mismatch envelope (quadratic fit). Disables
            % the term (with a warning) if the file or an axis' fit is missing.
            qf_path = FilterConfig.disturbance_quadfit_default;
            if ~isfile(qf_path)
                warning('LivenessFilter:disturbance', ...
                    'Disturbance quadfit "%s" not found; using nominal filter.', qf_path);
                obj.use_disturbance = false;
                return;
            end
            qf = load(qf_path);
            for i = 1:numel(obj.axes)
                axis = obj.axes{i};
                if ~isfield(qf, ['beta_' axis]) || ~isfield(qf, ['half_' axis])
                    warning('LivenessFilter:disturbance', ...
                        'Quadfit missing "%s" axis; using nominal filter.', axis);
                    obj.use_disturbance = false;
                    return;
                end
                obj.dist.(axis).beta = qf.(['beta_' axis]);      % 15 x 4 x UH
                obj.dist.(axis).half = qf.(['half_' axis])(:);   % 4 x 1
            end
        end

        function reset_counters(obj)
            obj.n_calls  = 0;
            obj.n_active = 0;
        end

        function [u, info] = filter(obj, frames, u_all_nom, U0)
            % frames    : struct keyed by axis; .curr / .next BRT anchor frames.
            % u_all_nom : 13x1 nominal effector perturbation (mission frame).
            % U0        : 13x1 mission trim input.
            % Returns the projected 13x1 u (anchor frame); caller maps back with
            % u - info.target_dtrim.
            obj.n_calls = obj.n_calls + 1;
            input_spec  = obj.combined_spec;
            active_axes = obj.axes;
            num_axes    = numel(active_axes);

            % Evaluate curr/next anchor candidates. UH scheduling is
            % axis-independent, so the curr-vs-next decision is shared.
            curr_by_axis = struct();
            next_by_axis = struct();
            same_anchor      = true;
            transition_ready = true;
            for i = 1:num_axes
                axis     = active_axes{i};
                value_fn = obj.value_function.(axis);
                curr = obj.evaluate_anchor(frames.(axis).curr, value_fn, [axis '_current']);
                next = obj.evaluate_anchor(frames.(axis).next, value_fn, [axis '_next']);
                curr_by_axis.(axis) = curr;
                next_by_axis.(axis) = next;
                same_anchor      = same_anchor      && (curr.idx == next.idx);
                transition_ready = transition_ready && (next.valid && next.V < 0);
            end
            use_next = same_anchor || transition_ready;

            target = struct();
            all_valid  = true;
            all_inside = true;
            for i = 1:num_axes
                axis = active_axes{i};
                if use_next, target.(axis) = next_by_axis.(axis); else, target.(axis) = curr_by_axis.(axis); end
                all_valid  = all_valid  && target.(axis).valid;
                all_inside = all_inside && target.(axis).inside_grid;
            end

            % Shared anchor trim (identical across axes at this anchor).
            U0a          = target.(active_axes{1}).U0a;
            trim_offset  = U0(input_spec.U0_idx) - U0a(input_spec.U0_idx);
            u_nom        = u_all_nom(:) + trim_offset;
            [lb, ub]     = LivenessFilter.input_bounds(U0a, input_spec);
            u0           = min(max(u_nom, lb), ub);

            V_all = struct();
            for i = 1:num_axes, V_all.(active_axes{i}) = target.(active_axes{i}).V; end

            info = struct( ...
                'target_dtrim', trim_offset, 'active', false, 'mode', obj.mode, ...
                'ok', false, 'inside_grid', all_inside, ...
                'V', target.(active_axes{1}).V, 'V_all', V_all, 'dVdt', NaN, 'rhs', NaN, ...
                'dist', 0, 'dist_all', struct(), ...
                'u_nom', u_nom, 'u0', u0, 'u', u_nom, 'du', 0, ...
                'sat_clip', norm(u0 - u_nom), 'command_changed', 0, ...
                'solver', 'none', 'quadprog_exitflag', NaN, 'lb', lb, 'ub', ub);

            if strcmp(obj.mode, 'off')
                u = u_nom; info.u = u; obj.last_info = info; return;
            end
            if ~strcmp(obj.mode, 'blend')
                error('LivenessFilter:mode', ...
                    'Liveness filter supports only ''blend'' or ''off'' (got "%s").', obj.mode);
            end

            % No guarantee unless every active axis is valid and in-grid.
            if ~all_valid
                u = u0;
                info.solver = 'no-coverage-box-clipped';
                info.u = u; info.du = norm(u - u0); info.command_changed = norm(u - u_nom);
                obj.last_info = info; return;
            end
            info.ok = true;

            % One liveness constraint per axis, over the shared input.
            Aineq     = zeros(num_axes, input_spec.nu);
            bineq     = zeros(num_axes, 1);
            alpha_vec = zeros(num_axes, 1);
            dist_vec  = zeros(num_axes, 1);
            for i = 1:num_axes
                axis = active_axes{i};
                tgt  = target.(axis);
                Bp_shared = zeros(size(tgt.Ap_a, 1), input_spec.nu);
                Bp_shared(:, input_spec.cols.(axis)) = tgt.Bp_a;
                alpha_vec(i) = tgt.gradV(:)' * (tgt.Ap_a * tgt.x_anchor(:));
                Aineq(i, :)  = (tgt.gradV(:)' * Bp_shared);
                if obj.use_disturbance
                    dbar = obj.disturbance_bound(axis, tgt.x_anchor, tgt.idx);
                    dist_vec(i) = abs(tgt.gradV(:))' * dbar;
                end
                bineq(i)     = -obj.gamma * (tgt.V + obj.live_margin) - alpha_vec(i) - dist_vec(i);
            end
            info.rhs = -obj.gamma * (target.(active_axes{1}).V + obj.live_margin);
            dist_all = struct();
            for i = 1:num_axes, dist_all.(active_axes{i}) = dist_vec(i); end
            info.dist     = dist_vec(1);
            info.dist_all = dist_all;

            if all(Aineq * u0 <= bineq + 1e-10)
                u = u0;
                info.active = false;
                info.solver = 'inactive-clipped-nominal';
            else
                [u, exitflag] = LivenessFilter.solve_blend_qp(u_nom, Aineq, bineq, lb, ub);
                % [u, exitflag] = LivenessFilter.solve_blend_qp(u_nom, Aineq, bineq);
                info.active = true;
                info.solver = 'quadprog';
                info.quadprog_exitflag = exitflag;
                if isempty(u) || exitflag <= 0
                    u = u0;
                    info.active = false;
                    info.solver = 'quadprog-failed-box-clipped';
                end
            end

            dVdt_all = struct();
            for i = 1:num_axes, dVdt_all.(active_axes{i}) = alpha_vec(i) + Aineq(i, :) * u; end
            info.dVdt     = dVdt_all.(active_axes{1});
            info.dVdt_all = dVdt_all;
            if info.active, obj.n_active = obj.n_active + 1; end
            info.du = norm(u - u0);
            info.command_changed = norm(u - u_nom);
            info.u = u;
            obj.last_info = info;
        end

        function candidate = evaluate_anchor(~, frame, value_fn, name)
            % Evaluate one BRT anchor frame: V, gradV and coverage flags.
            x = frame.x_anchor(:);
            inside_grid = all(x >= value_fn.grid_min(:) - 1e-12) && ...
                          all(x <= value_fn.grid_max(:) + 1e-12);

            [V, gradV, ok] = value_fn.query(x, frame.uhA, frame.whA);
            if isempty(gradV) || numel(gradV) ~= 4
                gradV = NaN(4, 1);
            end

            candidate = frame;
            candidate.name = name;
            candidate.x_anchor = x;
            candidate.V = V;
            candidate.gradV = gradV(:);
            candidate.ok = ok;
            candidate.inside_grid = inside_grid;
            candidate.valid = ok && inside_grid && isfinite(V);
        end

        function dbar = disturbance_bound(obj, axis, x, idx)
            % Per-channel worst-case mismatch magnitude (4x1) from the quadratic
            % envelope at UH index idx. Zeros outside the fit coverage.
            d = obj.dist.(axis);
            if idx < 1 || idx > size(d.beta, 3)
                dbar = zeros(4, 1);
                return;
            end
            z    = x(:) ./ d.half;                       % normalized deviation
            phi  = LivenessFilter.design_quad(z(:)');    % 1 x 15
            dbar = max((phi * d.beta(:, :, idx))', 0);   % 4 x 1
        end
    end

    methods (Static)
        function Phi = design_quad(Z)
            % Quadratic feature row(s): [1, z1..z4, z1^2..z4^2, pairwise] (15 terms).
            % Must match tests/diag_model_residual.m.
            z1 = Z(:,1); z2 = Z(:,2); z3 = Z(:,3); z4 = Z(:,4); o = ones(size(z1));
            Phi = [o, z1, z2, z3, z4, z1.^2, z2.^2, z3.^2, z4.^2, ...
                   z1.*z2, z1.*z3, z1.*z4, z2.*z3, z2.*z4, z3.*z4];
        end

        function [lb, ub] = input_bounds(U0, spec)
            % Per-effector perturbation bounds about the trim (from spec.U0_idx):
            %   lb_i = max(phys_lb, trim_i - Delta) - trim_i
            %   ub_i = min(phys_ub, trim_i + Delta) - trim_i
            trim = U0(spec.U0_idx);
            trim = trim(:);

            lift_lb = FilterConfig.rpm2radps(FilterConfig.Pi_lift_rpm(1));
            lift_ub = FilterConfig.rpm2radps(FilterConfig.Pi_lift_rpm(2));
            push_lb = FilterConfig.rpm2radps(FilterConfig.Pi_push_rpm(1));
            push_ub = FilterConfig.rpm2radps(FilterConfig.Pi_push_rpm(2));
            surf_lb = deg2rad(FilterConfig.srf_deg(1));
            surf_ub = deg2rad(FilterConfig.srf_deg(2));

            Dlift = FilterConfig.rpm2radps(FilterConfig.Delta_lift_RPM);
            Dpush = FilterConfig.rpm2radps(FilterConfig.Delta_push_RPM);
            Dsurf = deg2rad(FilterConfig.Delta_surf_Deg);

            nu = spec.nu;
            lb = zeros(nu, 1);
            ub = zeros(nu, 1);

            for i = 1:nu
                switch char(spec.effector_type(i))
                    case 'lift'
                        pl = lift_lb; pu = lift_ub; D = Dlift;
                    case 'push'
                        pl = push_lb; pu = push_ub; D = Dpush;
                    case 'surf'
                        pl = surf_lb; pu = surf_ub; D = Dsurf;
                    otherwise
                        error('LivenessFilter:effType', ...
                              'Unknown effector type "%s".', spec.effector_type(i));
                end

                lb(i) = max(pl, trim(i) - D) - trim(i);
                ub(i) = min(pu, trim(i) + D) - trim(i);
            end
        end

        function [u, exitflag] = solve_blend_qp(u_nom, Aineq, bineq, lb, ub)
            % min 0.5*||u-u_nom||^2  s.t.  Aineq*u <= bineq,  lb <= u <= ub
            u_nom = u_nom(:);
            lb    = lb(:);
            ub    = ub(:);
            nu    = numel(u_nom);

            H = eye(nu);
            f = -u_nom;

            opts = optimoptions('quadprog', ...
                                'Display', 'off', ...
                                'Algorithm', 'interior-point-convex', ...
                                'OptimalityTolerance', 1e-10, ...
                                'ConstraintTolerance', 1e-10, ...
                                'StepTolerance', 1e-12, ...
                                'MaxIterations', 200);

            [u, ~, exitflag] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], opts);
        end
    end
end
