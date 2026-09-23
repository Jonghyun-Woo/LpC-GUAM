classdef FilterConfig < handle
    % Configuration constants and per-axis spec for the longitudinal
    % (and, for verification only, lateral) HJ-reachability liveness filter.
    %
    % handle class: runtime knobs (mode/wh_anchor) set on the hub
    % (cfg.controller.filter) are shared by reference with the consumers.

    properties
        % --- Runtime knobs (set on the hub before LpC_GUAM construction) ---
        mode      = 'blend';         % 'blend' | 'off'  (liveness filter mode)
        wh_anchor = 0;               % BRT scheduling WH anchor [ft/s]; 0 = WH2
                                     %   (level, wh=0). []=Controller default WH(3)
        axes      = {'lon', 'lat'};  % active filter axes; one liveness constraint per
                                     %   axis, solved as a single QP over the shared
                                     %   13-effector input. {'lon'} = single-axis.
        use_disturbance = true;      % add the model-mismatch disturbance support term
                                     %   s_a = sum_c |gradV_c| * dbar_c(x) to each axis'
                                     %   CBF constraint, matching the disturbed BRT
                                     %   generation. false = nominal (optimistic) filter.
    end

    properties (Constant)
        % --- Filter parameters (axis-independent) ---
        gamma    = 5.1;     % smooth-blending CBF rate (paper recommends high gamma; tune post-integ)
        eps_band = 1e-3;    % LR boundary band: treat V >= -eps_band as boundary/outside (default-live)
        live_margin = 0.0;  % Conservative live-set margin c >= 0

        % --- Perturbation caps from trim (must match BRT generation yml) ---
        Delta_lift_RPM = 100;   % lift-rotor speed change from trim (RPM); aero-valid cap
        Delta_push_RPM = 300;   % pusher-rotor speed change from trim (RPM)
        Delta_surf_Deg = 15;    % control-surface deflection change from trim (deg)

        % --- Physical (absolute) control bounds ---
        Pi_lift_rpm = [0, 1600];    % lift-rotor speed physical bounds (RPM)
        Pi_push_rpm = [0, 2000];    % pusher-rotor speed physical bounds (RPM)
        srf_deg     = [-30, 30];    % control-surface deflection physical bounds (deg)

        % Parent directory holding per-axis BRT value-function subfolders
        % (<AXIS>_BRT), each with converged GUAM_<AXIS>_BRT_UH*_WH*.mat files.
        tables_dir_default = 'reachable_data/guam_output';

        % State-dependent model-mismatch disturbance envelope (quadratic fit),
        % as used in the BRT generation. Fields: beta_<axis> (15x4xUH),
        % half_<axis> (1x4), UH_LIST. See tests/diag_model_residual.m.
        disturbance_quadfit_default = 'reachable_data/mc_verify/guam_disturbance_quadfit.mat';
    end

    methods
        function obj = FilterConfig(overrides)
            if nargin < 1 || isempty(overrides), overrides = struct(); end
            obj.mode      = getfield_default(overrides, 'filter_mode', obj.mode);
            obj.wh_anchor = getfield_default(overrides, 'filter_wh_anchor', obj.wh_anchor);
            obj.axes      = getfield_default(overrides, 'filter_axes', obj.axes);
            obj.use_disturbance = getfield_default(overrides, 'filter_use_disturbance', obj.use_disturbance);
        end
    end

    methods (Static)
        function spec = axisSpec(axis)
            % Return the per-axis spec struct used by ValueFunction and
            % LivenessFilter. All grid data in ft/s, rad/s, rad.
            %
            %   axis : 'lon' | 'lat'
            %
            % Fields:
            %   grid_min, grid_max : 4x1 grid corner (grid_min = -grid_max).
            %                        Fallback only: ValueFunction reads the
            %                        authoritative grid from the BRT .mat files.
            %   grid_num           : 4x1 grid node counts N (fallback only)
            %   target_ub          : 4x1 target box upper corner (lb = -ub)
            %   brt_prefix         : BRT file-name prefix
            %   nu                 : number of physical effectors (filter dim)
            %   U0_idx             : indices into RSLQR U0 (13x1) picking the
            %                        per-effector trim inputs in filter order
            %   state_rows         : plant-state (12x1) rows giving [x1;x2;x3;x4]
            %   trim_rows          : X0 trim-state (12x1) rows giving [x1;x2;x3;x4]
            switch lower(axis)
                case 'lon'
                    spec.grid_max   = [16.4042; 22.9659; 0.6109; 0.4363];
                    spec.grid_num   = [33; 51; 49; 33];
                    spec.target_ub  = [3; 3; 0.10; 0.10];
                    spec.brt_prefix = 'GUAM_LON_BRT';
                    spec.nu         = 11;                    % [Pi_1..8, Pi_9(pusher), delta_e, delta_f]
                    spec.U0_idx     = [5:12, 13, 3, 1];      % lift1..8, pusher, elevator, flap
                    % Effector type per filter slot: 8x lift, 1x push, 2x surface
                    spec.effector_type = [repmat("lift", 1, 8), "push", "surf", "surf"];
                    spec.state_rows = [4; 6; 11; 8];         % plant state -> [u; w; q; theta]
                    spec.trim_rows  = [1; 3; 5; 11];         % X0 trim state -> [u; w; q; theta]
                case 'lat'
                    spec.grid_max   = [16.4042; 1.5708; 0.6109; 0.7854];
                    spec.grid_num   = [33; 85; 49; 39];
                    spec.target_ub  = [3; 0.10; 0.10; 0.10];
                    spec.brt_prefix = 'GUAM_LAT_BRT';
                    spec.nu         = 10;                    % [Pi_1..8, delta_a, delta_r]
                    spec.U0_idx     = [5:12, 2, 4];          % lift1..8, aileron, rudder
                    spec.effector_type = [repmat("lift", 1, 8), "surf", "surf"];
                    spec.state_rows = [5; 10; 12; 7];        % plant state -> [v; p; r; phi]
                    spec.trim_rows  = [2; 4; 6; 10];         % X0 trim state -> [v; p; r; phi]
                otherwise
                    error('FilterConfig:axisSpec', ...
                        'Unknown axis "%s" (expected ''lon'' or ''lat'').', axis);
            end
            spec.axis  = lower(axis);
            spec.grid_min = -spec.grid_max;
        end

        function cspec = combinedSpec()
            % Unified 13-effector perturbation input shared by the joint lon+lat
            % liveness filter, ordered exactly as the RSLQR perturb_cmd:
            %   [lift1..8, pusher, flap, aileron, elevator, rudder]
            % Reuses input_bounds via U0_idx/effector_type/nu. cols.<axis>
            % place each axis' Bp columns (in that axis' own filter order) into
            % this shared 13-column input space via cols.<axis>:
            %   lon Bp order [lift1..8, pusher, elevator, flap]
            %   lat Bp order [lift1..8, aileron, rudder]
            cspec.nu            = 13;
            cspec.effector_type = [repmat("lift", 1, 8), "push", "surf", "surf", "surf", "surf"];
            cspec.U0_idx        = [5:12, 13, 1, 2, 3, 4];   % trim from U0=[flap;ail;ele;rud;eng1..9]
            cspec.cols.lon      = [1:8, 9, 12, 10];
            cspec.cols.lat      = [1:8, 11, 13];
        end

        function out = rpm2radps(rpm)
            % Convert rotor speed from RPM to rad/s (matches GUAM_Config).
            out = rpm * 2 * pi / 60;
        end
    end
end
