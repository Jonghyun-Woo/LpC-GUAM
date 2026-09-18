classdef LivenessFilter < handle
    % HJ-reachability liveness filter operating behind a nominal controller.
    %
    % Handles one or more axes uniformly (obj.axes, e.g. {'lon'} or {'lon','lat'}).
    % Each active axis contributes ONE liveness constraint on the shared
    % 13-effector perturbation input u (combined_spec ordering = perturb_cmd):
    %
    %   dV_a/dt = alpha_a + beta_a'*u,  alpha_a = gradV_a'*(Ap_a*x_a),
    %                                   beta_a  = (gradV_a'*Bp_a_ext)'
    %
    % Blend mode solves the CBF-like projection with MATLAB quadprog:
    %
    %   min_u 0.5*||u - u_nom||^2
    %   s.t.  beta_a'*u <= -gamma*(V_a + live_margin) - alpha_a - s_a   (each axis a)
    %         lb <= u <= ub
    %
    % Each axis' Bp is embedded into the shared input space via combined_spec.cols,
    % so the shared lift columns are arbitrated jointly by the single QP.
    %
    % s_a is the worst-case model-mismatch disturbance term (use_disturbance),
    % matching the disturbed BRT generation:
    %   s_a = max_{|d_c|<=dbar_c} gradV_a'd = sum_c |gradV_a,c(x)| * dbar_a,c(x),
    % dbar_a(x) the per-channel quadratic envelope (FilterConfig quadfit). s_a=0
    % recovers the nominal (optimistic) filter.
    %
    % Units: ft/s, rad/s, rad.

    properties
        mode
        gamma
        live_margin

        axes            % cell list of active axes (e.g. {'lon'} or {'lon','lat'})
        value_function  % struct keyed by axis; each ValueFunction carries its .axis_spec
        combined_spec   % shared 13-effector input spec (FilterConfig.combinedSpec)

        use_disturbance % include the worst-case mismatch term s_a in each constraint
        dist            % struct keyed by axis: .beta (15x4xUH), .half (4x1)

        n_calls
        n_active
        last_info
    end

    methods
        function obj = LivenessFilter(axes, mode, tables_dir, uh_breakpoint, wh_breakpoint, use_disturbance)
            % axes : axis name or cell list. Each named axis gets its own value
            %        function; the liveness constraints are always solved jointly
            %        over the shared 13-effector input.
            % use_disturbance : add the worst-case mismatch term s_a (default false).
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
                ax   = obj.axes{i};
                spec = FilterConfig.axisSpec(ax);
                obj.value_function.(ax) = ValueFunction(spec, tables_dir, uh_breakpoint, wh_breakpoint);
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
            % Load the per-axis model-mismatch disturbance envelope (quadratic
            % fit) used in BRT generation. Disables the term (with a warning) if
            % the file or an axis' fit is missing.
            qf_path = FilterConfig.disturbance_quadfit_default;
            if ~isfile(qf_path)
                warning('LivenessFilter:disturbance', ...
                    'Disturbance quadfit "%s" not found; using nominal filter.', qf_path);
                obj.use_disturbance = false;
                return;
            end
            qf = load(qf_path);
            for i = 1:numel(obj.axes)
                ax = obj.axes{i};
                if ~isfield(qf, ['beta_' ax]) || ~isfield(qf, ['half_' ax])
                    warning('LivenessFilter:disturbance', ...
                        'Quadfit missing "%s" axis; using nominal filter.', ax);
                    obj.use_disturbance = false;
                    return;
                end
                obj.dist.(ax).beta = qf.(['beta_' ax]);      % 15 x 4 x UH
                obj.dist.(ax).half = qf.(['half_' ax])(:);   % 4 x 1
            end
        end

        function reset_counters(obj)
            obj.n_calls  = 0;
            obj.n_active = 0;
        end

        function [u, info] = filter(obj, frames, u_all_nom, U0)
            % frames    : struct keyed by axis; frames.(ax).curr / .next are BRT
            %             anchor frames (uhA, whA, X0a, U0a, Ap_a, Bp_a, x_anchor, idx).
            % u_all_nom : 13x1 nominal combined effector perturbation (mission frame).
            % U0        : 13x1 mission trim input.
            % Returns the projected 13x1 u (anchor frame) and diagnostics; the
            % caller maps back with u - info.target_dtrim.
            obj.n_calls = obj.n_calls + 1;
            cspec = obj.combined_spec;
            ax_list = obj.axes;
            n_ax = numel(ax_list);

            % Evaluate curr/next candidates per axis. The UH scheduling (idx) is
            % axis-independent, so the curr-vs-next anchor decision is shared.
            curr_cand = struct();
            next_cand = struct();
            same_idx = true;
            transition_ready = true;
            for i = 1:n_ax
                ax = ax_list{i};
                vf = obj.value_function.(ax);
                cc = obj.eval_candidate(frames.(ax).curr, vf, vf.axis_spec, [ax '_current']);
                cn = obj.eval_candidate(frames.(ax).next, vf, vf.axis_spec, [ax '_next']);
                curr_cand.(ax) = cc;
                next_cand.(ax) = cn;
                same_idx = same_idx && (cc.idx == cn.idx);
                transition_ready = transition_ready && (cn.valid && cn.V < 0);
            end
            use_next = same_idx || transition_ready;

            target = struct();
            all_valid  = true;
            all_inside = true;
            for i = 1:n_ax
                ax = ax_list{i};
                if use_next, target.(ax) = next_cand.(ax); else, target.(ax) = curr_cand.(ax); end
                all_valid  = all_valid  && target.(ax).valid;
                all_inside = all_inside && target.(ax).inside_grid;
            end

            % Shared anchor trim (identical across axes at this anchor).
            U0a       = target.(ax_list{1}).U0a;
            dtrim_all = U0(cspec.U0_idx) - U0a(cspec.U0_idx);
            u_nom     = u_all_nom(:) + dtrim_all;
            [lb, ub]  = LivenessFilter.input_bounds(U0a, cspec);
            u0        = min(max(u_nom, lb), ub);

            V_all = struct();
            for i = 1:n_ax, V_all.(ax_list{i}) = target.(ax_list{i}).V; end

            info = struct( ...
                'target_dtrim', dtrim_all, 'active', false, 'mode', obj.mode, ...
                'ok', false, 'inside_grid', all_inside, ...
                'V', target.(ax_list{1}).V, 'V_all', V_all, 'dVdt', NaN, 'rhs', NaN, ...
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

            % No formal guarantee unless every active axis is valid and in-grid.
            if ~all_valid
                u = u0;
                info.solver = 'no-coverage-box-clipped';
                info.u = u; info.du = norm(u - u0); info.command_changed = norm(u - u_nom);
                obj.last_info = info; return;
            end
            info.ok = true;

            % One liveness constraint per axis, over the shared input.
            Aineq     = zeros(n_ax, cspec.nu);
            bineq     = zeros(n_ax, 1);
            alpha_vec = zeros(n_ax, 1);
            dist_vec  = zeros(n_ax, 1);
            for i = 1:n_ax
                ax = ax_list{i};
                t = target.(ax);
                Bp_ext = zeros(size(t.Ap_a, 1), cspec.nu);
                Bp_ext(:, cspec.cols.(ax)) = t.Bp_a;
                alpha_vec(i) = t.gradV(:)' * (t.Ap_a * t.x_anchor(:));
                Aineq(i, :)  = (t.gradV(:)' * Bp_ext);
                if obj.use_disturbance
                    dbar = obj.disturbance_bound(ax, t.x_anchor, t.idx);
                    dist_vec(i) = abs(t.gradV(:))' * dbar;
                end
                bineq(i)     = -obj.gamma * (t.V + obj.live_margin) - alpha_vec(i) - dist_vec(i);
            end
            info.rhs = -obj.gamma * (target.(ax_list{1}).V + obj.live_margin);
            dist_all = struct();
            for i = 1:n_ax, dist_all.(ax_list{i}) = dist_vec(i); end
            info.dist     = dist_vec(1);
            info.dist_all = dist_all;

            if all(Aineq * u0 <= bineq + 1e-10)
                u = u0;
                info.active = false;
                info.solver = 'inactive-clipped-nominal';
            else
                [u, exitflag] = LivenessFilter.solve_blend_qp(u_nom, Aineq, bineq, lb, ub);
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
            for i = 1:n_ax, dVdt_all.(ax_list{i}) = alpha_vec(i) + Aineq(i, :) * u; end
            info.dVdt     = dVdt_all.(ax_list{1});
            info.dVdt_all = dVdt_all;
            if info.active, obj.n_active = obj.n_active + 1; end
            info.du = norm(u - u0);
            info.command_changed = norm(u - u_nom);
            info.u = u;
            obj.last_info = info;
        end

        function candidate = eval_candidate(~, brt_info, vf, spec, name)
            x = brt_info.x_anchor(:);

            inside_grid = all(x >= spec.grid_min(:) - 1e-12) && ...
                          all(x <= spec.grid_max(:) + 1e-12);

            [V, gradV, ok] = vf.query(x, brt_info.uhA, brt_info.whA);

            if isempty(gradV) || numel(gradV) ~= 4
                gradV = NaN(4, 1);
            end

            candidate = brt_info;
            candidate.name = name;
            candidate.x_anchor = x;
            candidate.V = V;
            candidate.gradV = gradV(:);
            candidate.ok = ok;
            candidate.inside_grid = inside_grid;
            candidate.valid = ok && inside_grid && isfinite(V);
        end

        function dbar = disturbance_bound(obj, ax, x, idx)
            % Per-channel worst-case mismatch magnitude dbar_c(x) (4x1) from the
            % quadratic envelope at UH index idx. Returns zeros outside the fit
            % coverage. x is the deviation-from-trim state (center = 0).
            d = obj.dist.(ax);
            if idx < 1 || idx > size(d.beta, 3)
                dbar = zeros(4, 1);
                return;
            end
            z   = x(:) ./ d.half;                       % normalized deviation
            phi = LivenessFilter.design_quad(z(:)');    % 1 x 15
            dbar = max((phi * d.beta(:, :, idx))', 0);  % 4 x 1
        end
    end

    methods (Static)
        function Phi = design_quad(Z)
            % Quadratic feature row(s) matching tests/diag_model_residual.m:
            %   [1, z1..z4, z1^2..z4^2, pairwise products]  (15 terms).
            z1 = Z(:,1); z2 = Z(:,2); z3 = Z(:,3); z4 = Z(:,4); o = ones(size(z1));
            Phi = [o, z1, z2, z3, z4, z1.^2, z2.^2, z3.^2, z4.^2, ...
                   z1.*z2, z1.*z3, z1.*z4, z2.*z3, z2.*z4, z3.*z4];
        end

        function [lb, ub] = input_bounds(U0, spec)
            % Per-effector perturbation bounds:
            %   lb_i = max(phys_lb_i, trim_i - Delta_i) - trim_i
            %   ub_i = min(phys_ub_i, trim_i + Delta_i) - trim_i
            % trim inputs are picked from RSLQR U0 by spec.U0_idx.
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
            % Solve:  min 0.5*||u-u_nom||^2  s.t.  Aineq*u <= bineq,  lb <= u <= ub
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
