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

        % Compatibility fields.
        % Existing scripts may still set f.rate_mode or inject fdyn through
        % RSLQR.set_liveness_dynamics(). This version ignores both.
        rate_mode
        fdyn

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

            % Kept for backward compatibility only.
            obj.rate_mode = 'linear';
            obj.fdyn      = [];

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

        function [u, info] = filter(obj, x, u_nom, A, B, U0, uh, wh, x_full) %#ok<INUSD>
            % x     : 4x1 perturbation state
            % u_nom : nu x1 nominal effector perturbation
            % A,B   : linear reduced dynamics matrices
            % U0    : 13x1 trim input for input-bound computation
            % uh,wh : BRT table scheduling variables

            if nargin < 9
                x_full = []; %#ok<NASGU>
            end

            x     = x(:);
            u_nom = u_nom(:);

            obj.n_calls = obj.n_calls + 1;

            % Input bounds are computed early
            [lb, ub] = LivenessFilter.input_bounds(U0, obj.spec);
            u0 = min(max(u_nom, lb), ub);

            inside_grid = all(x >= obj.spec.grid_min(:) - 1e-12) && ...
                          all(x <= obj.spec.grid_max(:) + 1e-12);

            info = struct( ...
                'V', NaN, ...
                'active', false, ...
                'mode', obj.mode, ...
                'rate_mode', 'linear', ...
                'ok', false, ...
                'inside_grid', inside_grid, ...
                'dVdt', NaN, ...
                'rhs', NaN, ...
                'dV0', NaN, ...
                'du', 0, ...
                'sat_clip', norm(u0 - u_nom), ...
                'command_changed', 0, ...
                'u_nom', u_nom, ...
                'u0', u0, ...
                'u', u_nom, ...
                'alpha', NaN, ...
                'beta_norm', NaN, ...
                'b', NaN, ...
                'dV_lin_nom', NaN, ...
                'dV_lin_min', NaN, ...
                'feasible', NaN, ...
                'solver', 'none', ...
                'quadprog_exitflag', NaN, ...
                'quadprog_output', [] ...
            );

            % OFF mode: preserve old pass-through behavior.
            if strcmp(obj.mode, 'off')
                u = u_nom;
                info.u = u;
                obj.last_info = info;
                return;
            end

            if ~inside_grid
                u = u_nom;
                info.solver = 'outside-grid-pass-through';
                info.u = u;
                info.command_changed = norm(u - u_nom);
                obj.last_info = info;
                return;
            end

            [V, gradV, ok] = obj.lut.query(x, uh, wh);
            info.V  = V;
            info.ok = ok;

            if ~ok
                u = u_nom;
                info.solver = 'no-table-coverage-pass-through';
                info.u = u;
                info.command_changed = norm(u - u_nom);
                obj.last_info = info;
                return;
            end

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
                    rhs = -obj.gamma * (V + obj.live_margin);
                    b   = rhs - alpha;

                    info.rhs = rhs;
                    info.b   = b;

                    info.dV0 = alpha + beta' * u0;
                    info.dV_lin_nom = info.dV0;

                    u_min_rate = LivenessFilter.minimize_linear_over_box(beta, lb, ub);
                    info.dV_lin_min = alpha + beta' * u_min_rate;
                    info.feasible = (info.dV_lin_min <= rhs + 1e-9);

                    [u, exitflag, output] = LivenessFilter.solve_blend_quadprog( ...
                        u_nom, beta, b, lb, ub);

                    info.active = true;
                    info.solver = 'quadprog';
                    info.quadprog_exitflag = exitflag;
                    info.quadprog_output = output;

                    if isempty(u) || exitflag <= 0
                        u = u_min_rate;
                        info.solver = 'quadprog-failed-minrate-backup';
                        info.feasible = false;
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

        function [u, exitflag, output] = solve_blend_quadprog(u_nom, beta, b, lb, ub)
            % Solve:
            %
            %   min_u 0.5*||u-u_nom||^2
            %   s.t.  beta'*u <= b
            %         lb <= u <= ub
            %
            % quadprog form:
            %
            %   min 0.5*u'*H*u + f'*u

            if exist('quadprog', 'file') ~= 2
                error('LivenessFilter:quadprogMissing', ...
                      ['quadprog was not found. Install MATLAB Optimization Toolbox ', ...
                       'or use the previous toolbox-free LivenessFilter.']);
            end

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

            [u, ~, exitflag, output] = quadprog( ...
                H, f, Aineq, bineq, [], [], lb, ub, [], opts);
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