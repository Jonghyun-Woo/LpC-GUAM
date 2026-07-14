classdef LivenessFilter < handle
    % HJ-reachability liveness filter operating behind a nominal controller.
    %
    % Linear-rate only version:
    %   dV/dt = alpha + beta'*u
    %   alpha = gradV'*(A*x)
    %   beta  = (gradV'*B)'
    %
    % Blend mode solves the CBF-like projection with MATLAB quadprog:
    %
    %   min_u 0.5*||u-u_nom||^2
    %   s.t.  beta'*u <= -gamma*(V + live_margin) - alpha
    %         lb <= u <= ub
    %
    % Units: ft/s, rad/s, rad.

    properties
        spec
        mode
        gamma
        eps_band
        live_margin

        lut

        n_calls
        n_active
        last_info
    end

    methods
        function obj = LivenessFilter(channel, mode, tables_dir, uh_bp, wh_bp)
            obj.spec        = LivenessConfig.channelSpec(channel);
            obj.mode        = lower(mode);
            obj.gamma       = LivenessConfig.gamma;
            obj.eps_band    = LivenessConfig.eps_band;
            obj.live_margin = LivenessConfig.live_margin;

            if nargin < 3 || isempty(tables_dir)
                tables_dir = LivenessConfig.tables_dir_default;
            end

            obj.lut = ValueFunctionLUT(obj.spec, tables_dir, uh_bp, wh_bp);

            obj.n_calls  = 0;
            obj.n_active = 0;
            obj.last_info = [];
        end

        function reset_counters(obj)
            obj.n_calls  = 0;
            obj.n_active = 0;
        end

        function [u, info] = filter(obj, current_brt_info, next_brt_info, x_full) %#ok<INUSD>
            % info_struct.uhA  : BRT table scheduling variables
            % info_struct.whA  : BRT table scheduling variables
            % info_struct.X0a  : 12x1 trim state for input-bound computation
            % info_struct.U0a  : 13x1 trim input for input-bound computation
            % info_struct.Ap_a : linear reduced dynamics matrices
            % info_struct.Bp_a : linear reduced dynamics matrices
            % info_struct.x_anchor : : 4x1 perturbation state
            % info_struct.dtrim : trim perturbation from mission frame (Nominal Controller)
            % info_struct.u_anchor_nom : nu x1 nominal effector perturbation
            % info_struct.idx : BRT index
            
            obj.n_calls = obj.n_calls + 1;

            % Trim evaluation
            c = obj.eval_brt_candidate(current_brt_info, 'current');
            n = obj.eval_brt_candidate(next_brt_info, 'next');
            
            % Target Trim selection (Next 기준으로 수행)
            transition_ready = n.valid && n.V < 0;

            if c.idx == n.idx
                target = n;
            elseif transition_ready
                target = n;
            else
                target = c;
            end

            x     = target.x_anchor(:);
            uh    = target.uhA;
            wh    = target.whA;
            U0    = target.U0a;
            A     = target.Ap_a;
            B     = target.Bp_a;
            u_nom = target.u_anchor_nom(:);
            dtrim = target.dtrim(:);

            % Input bounds are computed early
            [lb, ub] = LivenessFilter.input_bounds(U0, obj.spec);
            u0 = min(max(u_nom, lb), ub);
            
            % Initialization info struct
            info = struct( ...
                'target_dtrim', dtrim, ...
                'active', false, ...
                'mode', obj.mode, ...
                'ok', false, ...
                'inside_grid', target.inside_grid, ...
                'V', NaN, ...
                'dVdt', NaN, ...
                'rhs', NaN, ...
                'u_nom', u_nom, ...
                'u0', u0, ...
                'u', u_nom, ...
                'du', 0, ...
                'sat_clip', norm(u0 - u_nom), ...
                'command_changed', 0, ...
                'solver', 'none', ...
                'quadprog_exitflag', NaN ...
            );

            % ---------------------------------------------------------------------
            % OFF mode: preserve nominal input in the selected frame.
            % ---------------------------------------------------------------------
            if strcmp(obj.mode, 'off')
                u = u_nom;
                info.u = u;
                info.command_changed = norm(u - u_nom);
                obj.last_info = info;
                return;
            end

            % ---------------------------------------------------------------------
            % Invalid selected BRT frame.
            % No formal BRT/CBF guarantee. Use box-clipped fallback.
            % ---------------------------------------------------------------------
            if ~target.inside_grid
                u = u0;
                info.solver = 'outside-grid-box-clipped';
                info.u = u;
                info.du = norm(u - u0);
                info.command_changed = norm(u - u_nom);
                obj.last_info = info;
                return;
            end

            if ~target.ok
                u = u0;
                info.solver = 'no-table-coverage-box-clipped';
                info.u = u;
                info.du = norm(u - u0);
                info.command_changed = norm(u - u_nom);
                obj.last_info = info;
                return;
            end
            
            % ---------------------------------------------------------------------
            % Valid selected BRT frame.
            % --------------------------------------------------------------------- 
            V     = target.V;
            gradV = target.gradV(:);

            info.V  = V;
            info.ok = true;
            info.inside_grid = true;

            alpha = gradV(:)' * (A * x);
            beta  = (gradV(:)' * B)';

            info.alpha     = alpha;
            info.beta_norm = norm(beta);

            switch lower(obj.mode)
                case 'lr'
                    % Linear least-restrictive style backup.
                    if V < -(obj.live_margin + obj.eps_band)
                        u = u0;
                        info.active = false;
                        info.solver = 'lr-inactive';
                    else
                        u = LivenessFilter.minimize_linear_over_box(beta, lb, ub);
                        info.active = true;
                        info.solver = 'lr-bangbang';
                    end

                    info.rhs = 0;
                    info.dVdt = alpha + beta' * u;

                case 'blend'
                    % CBF-like smooth blending filter.
                    rhs = -obj.gamma * (V + obj.live_margin);
                    b   = rhs - alpha;

                    info.rhs = rhs;
                    info.b   = b;

                    info.dV0 = alpha + beta' * u0;
                    info.dV_lin_nom = info.dV0;

                    u_min_rate = LivenessFilter.minimize_linear_over_box(beta, lb, ub);
                    info.dV_lin_min = alpha + beta' * u_min_rate;
                    info.feasible = (info.dV_lin_min <= rhs + 1e-9);

                    % If box-clipped nominal already satisfies the constraint,
                    % do not count this as active QP intervention.
                    if beta' * u0 <= b + 1e-10
                        u = u0;
                        info.active = false;
                        info.solver = 'inactive-clipped-nominal';
                    else
                        [u, exitflag] = LivenessFilter.solve_blend_quadprog(u_nom, beta, b, lb, ub);

                        info.active = true;
                        info.solver = 'quadprog';
                        info.quadprog_exitflag = exitflag;

                        if isempty(u) || exitflag <= 0
                            u = u_min_rate;
                            info.solver = 'quadprog-failed-minrate-backup';
                            info.feasible = false;
                        end
                    end

                    info.dVdt = alpha + beta' * u;

                otherwise
                    error('LivenessFilter:mode', ...
                        'Unknown mode "%s" (expected off, lr, or blend).', obj.mode);
            end

            if info.active
                obj.n_active = obj.n_active + 1;
            end

            info.du = norm(u - u0);
            info.command_changed = norm(u - u_nom);
            info.u = u;

            obj.last_info = info;
        end

        function cand = eval_brt_candidate(obj, brt_info, name)
            x = brt_info.x_anchor(:);

            inside_grid = all(x >= obj.spec.grid_min(:) - 1e-12) && ...
                          all(x <= obj.spec.grid_max(:) + 1e-12);

            [V, gradV, ok] = obj.lut.query(x, brt_info.uhA, brt_info.whA);

            if isempty(gradV) || numel(gradV) ~= 4
                gradV = NaN(4, 1);
            end

            cand = brt_info;
            cand.name = name;
            cand.x_anchor = x;
            cand.V = V;
            cand.gradV = gradV(:);
            cand.ok = ok;
            cand.inside_grid = inside_grid;
            cand.valid = ok && inside_grid && isfinite(V);
        end

        function s = strip_candidate(~, cand)
            % Lightweight diagnostic copy for last_info.
            s = struct();
            s.name = cand.name;
            s.idx = cand.idx;
            s.uhA = cand.uhA;
            s.whA = cand.whA;
            s.V = cand.V;
            s.ok = cand.ok;
            s.inside_grid = cand.inside_grid;
            s.valid = cand.valid;
            s.x_anchor = cand.x_anchor;
        end
    end

    methods (Static)
        function [lb, ub] = input_bounds(U0, spec)
            % Per-effector perturbation bounds:
            %
            %   lb_i = max(phys_lb_i, trim_i - Delta_i) - trim_i
            %   ub_i = min(phys_ub_i, trim_i + Delta_i) - trim_i
            %
            % trim inputs are picked from RSLQR U0 by spec.U0_idx.

            trim = U0(spec.U0_idx);
            trim = trim(:);

            lift_lb = LivenessConfig.rpm2radps(LivenessConfig.Pi_lift_rpm(1));
            lift_ub = LivenessConfig.rpm2radps(LivenessConfig.Pi_lift_rpm(2));
            push_lb = LivenessConfig.rpm2radps(LivenessConfig.Pi_push_rpm(1));
            push_ub = LivenessConfig.rpm2radps(LivenessConfig.Pi_push_rpm(2));
            surf_lb = deg2rad(LivenessConfig.srf_deg(1));
            surf_ub = deg2rad(LivenessConfig.srf_deg(2));

            Dlift = LivenessConfig.rpm2radps(LivenessConfig.Delta_lift_RPM);
            Dpush = LivenessConfig.rpm2radps(LivenessConfig.Delta_push_RPM);
            Dsurf = deg2rad(LivenessConfig.Delta_surf_Deg);

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

        function [u, exitflag] = solve_blend_quadprog(u_nom, beta, b, lb, ub)
            % Solve:
            %
            %   min_u 0.5*||u-u_nom||^2
            %   s.t.  beta'*u <= b
            %         lb <= u <= ub
            %
            % quadprog form:
            %
            %   min 0.5*u'*H*u + f'*u

            u_nom = u_nom(:);
            beta  = beta(:);
            lb    = lb(:);
            ub    = ub(:);

            nu = numel(u_nom);

            H = eye(nu);
            H = 0.5 * (H + H');
            f = -u_nom;

            Aineq = beta(:)';
            bineq = b;

            opts = optimoptions('quadprog', ...
                                'Display', 'off', ...
                                'Algorithm', 'interior-point-convex', ...
                                'OptimalityTolerance', 1e-10, ...
                                'ConstraintTolerance', 1e-10, ...
                                'StepTolerance', 1e-12, ...
                                'MaxIterations', 200);

            [u, ~, exitflag] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], opts);
        end

        function u = minimize_linear_over_box(a, lb, ub)
            % Solve:
            %
            %   min a'*u
            %   s.t. lb <= u <= ub
            %
            % Componentwise solution.

            a  = a(:);
            lb = lb(:);
            ub = ub(:);

            u = (a >= 0) .* lb + (a < 0) .* ub;
        end
    end
end