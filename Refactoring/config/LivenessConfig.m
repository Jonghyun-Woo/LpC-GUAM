classdef LivenessConfig
    % Configuration constants and per-channel spec for the longitudinal
    % (and, for verification only, lateral) HJ-reachability liveness filter.
    %
    % Mirrors the helperOC BRT-generation constants (Config.m, GUAM_Config.m)
    % so that the runtime filter's control-authority box matches the box the
    % value functions were computed under.
    %
    % UNITS: grids, states and inputs are ft/s, rad/s, rad. The m/deg scaling
    % found in helperOC visualize_*_tube.m is VISUALIZATION-ONLY and must NOT
    % be applied to the value-function grid.
    %
    % Reference: docs/refs/liveness-filter.md, helperOC Config.m / GUAM_Config.m.

    properties (Constant)
        % --- Filter parameters (channel-independent) ---
        gamma    = 5.0;     % smooth-blending CBF rate (paper recommends high gamma; tune post-integ)
        eps_band = 1e-3;    % LR boundary band: treat V >= -eps_band as boundary/outside (default-live)
        % dV/dt evaluation for the blend trigger/QP:
        %   'linear'    : dV/dt = alpha + beta'*u  (alpha=gradV'*A*x, beta=(gradV'*B)')
        %   'nonlinear' : dV/dt = gradV' * f_plant(x_full, u)  via an injected
        %                 plant-rate closure (see RSLQR.set_liveness_dynamics).
        % The linear model sign-flips at large perturbation (w>~12 ft/s), blinding
        % the trigger; 'nonlinear' evaluates the true rate. Falls back to 'linear'
        % if no closure / full state is available.
        rate_mode = 'linear';
        % Conservative live-set margin c >= 0: tighten the enforced live set to
        % {V <= -c} so the filter engages before the true boundary (V=0),
        % holding V <= -c. Constraint becomes dV/dt <= -gamma*(V+c); LR boundary
        % becomes V >= -(c+eps_band). c=0 reproduces the untightened filter.
        live_margin = 0.0;

        % --- Perturbation caps from trim (must match BRT generation) ---
        Delta_lift_RPM = 300;   % max lift-rotor speed change from trim (RPM)
        Delta_push_RPM = 300;   % max pusher-rotor speed change from trim (RPM)
        Delta_surf_Deg = 15;    % max control-surface deflection change from trim (deg)

        % --- Physical (absolute) control bounds ---
        Pi_lift_rpm = [0, 1600];    % lift-rotor speed physical bounds (RPM)
        Pi_push_rpm = [0, 2000];    % pusher-rotor speed physical bounds (RPM)
        srf_deg     = [-30, 30];    % control-surface deflection physical bounds (deg)

        % Default directory scanned for BRT value-function .mat files.
        tables_dir_default = 'tables/BRT';
    end

    methods (Static)
        function spec = channelSpec(channel)
            % Return the per-channel spec struct used by ValueFunctionLUT and
            % LivenessFilter. All grid data in ft/s, rad/s, rad.
            %
            %   channel : 'lon' (production) | 'lat' (verification-only)
            %
            % Fields:
            %   grid_min, grid_max : 4x1 grid corner (grid_min = -grid_max)
            %   grid_num           : 4x1 grid node counts N
            %   target_ub          : 4x1 target box upper corner (lb = -ub)
            %   brt_prefix         : BRT file-name prefix
            %   nu                 : number of physical effectors (filter dim)
            %   U0_idx             : indices into RSLQR U0 (13x1) picking the
            %                        per-effector trim inputs in filter order
            switch lower(channel)
                case 'lon'
                    spec.grid_max   = [16; 33; 1.5; 0.75];
                    spec.grid_num   = [21; 41; 61; 31];
                    spec.target_ub  = [3; 1.5; 0.05; 0.05];
                    spec.brt_prefix = 'GUAM_LON_BRT';
                    spec.nu         = 11;                    % [Pi_1..8, Pi_9(pusher), delta_e, delta_f]
                    spec.U0_idx     = [5:12, 13, 3, 1];      % lift1..8, pusher, elevator, flap
                    % Effector type per filter slot: 8x lift, 1x push, 2x surface
                    spec.effector_type = [repmat("lift", 1, 8), "push", "surf", "surf"];
                case 'lat'
                    spec.grid_max   = [13; 3; 0.35; 0.55];
                    spec.grid_num   = [21; 101; 21; 21];
                    spec.target_ub  = [3; 0.05; 0.05; 0.05];
                    spec.brt_prefix = 'GUAM_LAT_BRT';
                    spec.nu         = 10;                    % [Pi_1..8, delta_a, delta_r]
                    spec.U0_idx     = [5:12, 2, 4];          % lift1..8, aileron, rudder
                    spec.effector_type = [repmat("lift", 1, 8), "surf", "surf"];
                otherwise
                    error('LivenessConfig:channelSpec', ...
                          'Unknown channel "%s" (expected ''lon'' or ''lat'').', channel);
            end
            spec.channel  = lower(channel);
            spec.grid_min = -spec.grid_max;
        end

        function out = rpm2radps(rpm)
            % Convert rotor speed from RPM to rad/s (matches GUAM_Config).
            out = rpm * 2 * pi / 60;
        end
    end
end
