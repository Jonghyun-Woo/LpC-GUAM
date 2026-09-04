classdef RSLQR < handle
    % Robust Servo-mechanism LQR controller for Lift+Cruise transition.
    %
    % Reference: Acheson et al., "Examination of Unified Control Approaches
    % Incorporating Generalized Control Allocation", AIAA 2021-0999
    %
    % Two-phase usage:
    %   Design  (offline) : obj = RSLQR(aircraft, XEQ_TABLE, FreeVar_Table, Trans_Table, rho, grav)
    %   Runtime (online)  : [lon_gains, lat_gains] = obj.interpolate(UH_curr, WH_curr)
    %
    % Depends on (must be on MATLAB path):
    %   RSLQRConfig, poly_aero_wrapper_Mod_du (offline gain design only)

    properties
        vehicleCfg  % VehicleConfig (mass/inertia; offline gain design only)
        rslqrCfg    % RSLQRConfig (weights and dimension constants)
        dt          % servo-compensator discretization step [s] (injected from SimConfig)
        UH          % N_trim x 1  body-frame x-velocity (u) breakpoints [ft/s]
        WH          % M_trim x 1  body-frame z-velocity (w) breakpoints [ft/s]
                    % (verified: WH == XU0(3)=body w exactly; NOT heading/inertial.
                    %  At cruise pitch, body w>0 is ~level flight, not descent.)
        XU0
        LON       % 1x(N*M*L) struct array — longitudinal gains per schedule point
        LAT       % 1x(N*M*L) struct array — lateral gains per schedule point

        ctrl_lon_i  % Longitudinal integral term of controller
        ctrl_lat_i  % Lateral integral term of controller

        phi_v
        theta_v

        % Outer position-loop gain, Vb_cmd = ref_vel - Vb_0 + k_pos*e_pos_body.
        % It used to be the literal 0.1 in perturb_lin_ctrl; the default here
        % is that same 0.1, so nothing changes unless a caller sets it. Exposed
        % only so the slow altitude/speed oscillation can be tested against it:
        % if its period moves with this gain the loop is responsible, if the
        % period stays put the mode belongs to the airframe.
        k_pos = 0.1

        % Level-flight vertical-speed feedforward. Body w for level flight is
        % u*tan(theta), which is 0 at hover and ~13 ft/s at cruise, but the
        % commanded w perturbation is only k_pos*e_pos_body(3). With no other
        % source, the position error itself has to supply that 7-15 ft/s, and
        % since the outer loop is proportional-only it can only do that by
        % standing off the reference altitude: measured 151 ft of offset at
        % 100 ft/s. Supplying the term directly removes the reason for the
        % offset to exist.
        %   0 : off (shipped behaviour, nothing changes)
        %   1 : feedforward from the commanded speed and the trim attitude
        %   2 : same but using the flown attitude (becomes feedback on theta)
        % Mode 1 uses the trim attitude, which the vehicle is up to 14 deg away
        % from mid-transition, so it feeds the wrong correction and drives the
        % mission altitude to -90 ft. Mode 2 reads the attitude actually flown
        % and is what is used: standing altitude offset 103 -> 0 ft at a fixed
        % 100 ft/s, mission position error 11.5 -> 9.6 ft. Note this makes the
        % term feedback on theta rather than a true feedforward; it does not
        % diverge on this mission but its stability is not proven.
        w_ff_mode = 2
    end

    methods

        function obj = RSLQR(cfg, dt)
            % Load the precomputed gain-scheduled RSLQR design
            % (trim table + lon/lat gain tables over the UH x WH envelope).
            % trim_table_Poly_ConcatVer4p0.mat must be on the MATLAB path.
            %
            % Injected config (from the central hub):
            %   rslqrCfg  : RSLQRConfig  (weights/limits)
            %   dt        : sim timestep [s] (servo-compensator discretization)

            obj.vehicleCfg  = VehicleConfig();  % mass/inertia for offline gain design
            obj.rslqrCfg    = cfg.rslqr;
            obj.dt  = dt;
            S       = load('trim_table_Poly_ConcatVer4p0.mat');
            obj.XU0 = S.XU0_interp;
            obj.UH  = S.UH;
            obj.WH  = S.WH;

            % Gain source: keep the stored gains, or redesign them from the
            % current Q/R/W weights.
            if obj.rslqrCfg.update_gains
                S = obj.set_gains();
            end

            obj.LON.Ki  = S.Ki_lon_interp;
            obj.LON.Kx  = S.Kx_lon_interp;
            obj.LON.Kv  = S.Kv_lon_interp;
            obj.LON.F   = S.F_lon_interp;
            obj.LON.G   = S.G_lon_interp;
            obj.LON.C   = S.C_lon_interp;
            obj.LON.Cv  = S.Cv_lon_interp;
            obj.LON.W   = S.W_lon_interp;
            obj.LON.B   = S.B_lon_interp;
            obj.LON.Ap  = S.Ap_lon_interp;
            obj.LON.Bp  = S.Bp_lon_interp;
            
            obj.LAT.Ki  = S.Ki_lat_interp;
            obj.LAT.Kx  = S.Kx_lat_interp;
            obj.LAT.Kv  = S.Kv_lat_interp;
            obj.LAT.F   = S.F_lat_interp;
            obj.LAT.G   = S.G_lat_interp;
            obj.LAT.C   = S.C_lat_interp;
            obj.LAT.Cv  = S.Cv_lat_interp;
            obj.LAT.W   = S.W_lat_interp;
            obj.LAT.B   = S.B_lat_interp;
            obj.LAT.Ap  = S.Ap_lat_interp;
            obj.LAT.Bp  = S.Bp_lat_interp;

            obj.ctrl_lon_i = zeros([3, 1]);
            obj.ctrl_lat_i = zeros([3, 1]);

            obj.phi_v = 0;
            obj.theta_v = 0;
        end

        function [iU, kU, iW, kW] = lookup_breakpt(obj, uh_curr, wh_curr)
            % Clamp (uh_curr, wh_curr) to breakpoint bounds and return the
            % lower-cell indices and normalised fractional positions within
            % the enclosing breakpoint cell.
            %
            % Inputs
            %   uh_curr : scalar  current heading-frame x-velocity [ft/s]
            %   wh_curr : scalar  current heading-frame z-velocity [ft/s]
            %
            % Outputs
            %   iU : lower UH breakpoint index  (iU+1 is the upper index)
            %   kU : fractional position in [0,1] along the UH interval
            %   iW : lower WH breakpoint index  (iW+1 is the upper index)
            %   kW : fractional position in [0,1] along the WH interval
            uh_curr = max(obj.UH(1),  min(obj.UH(end),  uh_curr));
            wh_curr = max(obj.WH(1),  min(obj.WH(end),  wh_curr));
    
            iU = min(find(obj.UH <= uh_curr, 1, 'last'), length(obj.UH) - 1);
            iW = min(find(obj.WH <= wh_curr, 1, 'last'), length(obj.WH) - 1);
    
            kU = (uh_curr - obj.UH(iU)) / (obj.UH(iU+1) - obj.UH(iU));
            kW  = (wh_curr - obj.WH(iW)) / (obj.WH(iW+1) - obj.WH(iW));
        end

        function [X0, U0] = interp_xu0(obj, uh_curr, wh_curr)
            % Bilinearly interpolate the stored trim table obj.XU0 at
            % (uh_curr, wh_curr) and split the result into state and input.
            %
            % obj.XU0 is 25 x N_uh x N_wh:
            %   rows 1-12  : trim state vector  [u; v; w; p; q; r; ax; ay; az; phi; theta; psi]
            %   rows 13-25 : trim input vector  [flap; aileron; elevator; rudder; lift rotor speeds (1:8); pusher rotor speed (9)]
            %
            % Inputs
            %   uh_curr : scalar  current heading-frame x-velocity [ft/s]
            %   wh_curr : scalar  current heading-frame z-velocity [ft/s]
            %
            % Outputs
            %   X0 : 12x1  trim state vector at (uh_curr, wh_curr)
            %   U0 : 13x1  trim input vector at (uh_curr, wh_curr)
            [iU, kU, iW, kW] = obj.lookup_breakpt(uh_curr, wh_curr);
            xu0 = (1 - kU) * (1 - kW) * obj.XU0(:, iU,   iW  ) + ...
                       kU  * (1 - kW) * obj.XU0(:, iU+1, iW  ) + ...
                  (1 - kU) *      kW  * obj.XU0(:, iU,   iW+1) + ...
                       kU  *      kW  * obj.XU0(:, iU+1, iW+1);

            X0 = xu0(1:12);
            U0 = xu0(13:25);
        end

        function matrix_out = interp_mtrx(obj, matrix_table, uh_curr, wh_curr)
            % Bilinearly interpolate a matrix-valued schedule table at
            % (uh_curr, wh_curr).  Used to interpolate gain matrices (Ki,
            % Kx, B, etc.) stored as m x n x N_uh x N_wh arrays.
            %
            % Inputs
            %   matrix_table : m x n x N_uh x N_wh  precomputed matrix table
            %   uh_curr      : scalar  current heading-frame x-velocity [ft/s]
            %   wh_curr      : scalar  current heading-frame z-velocity [ft/s]
            %
            % Output
            %   matrix_out : m x n  interpolated matrix at (uh_curr, wh_curr)
            [iU, kU, iW, kW] = obj.lookup_breakpt(uh_curr, wh_curr);
            matrix_out = (1 - kU) * (1 - kW) * matrix_table(:, :, iU,   iW  ) + ...
                              kU  * (1 - kW) * matrix_table(:, :, iU+1, iW  ) + ...
                         (1 - kU) *      kW  * matrix_table(:, :, iU,   iW+1) + ...
                              kU  *      kW  * matrix_table(:, :, iU+1, iW+1);
        end

        function [e_pos_body, e_chi] = guidance_error(obj, state, ref)
            % In RBD, state is defined as [rn, re, rd, u, v, w, phi, theta, psi, p, q, r]' (12x1)
            ref_pos = ref.pos;
            ref_chi = ref.chi;

            phi     = state(7);
            theta   = state(8);
            psi     = state(9);
            i2b_rotm = obj.rotm_i2b(phi, theta, psi);

            e_pos_body = i2b_rotm * (ref_pos - state(1:3));

            s_chi_ref = sin(ref_chi);
            c_chi_ref = cos(ref_chi);
            s_psi = sin(psi);
            c_psi = cos(psi);

            y = (s_chi_ref * c_psi) - (c_chi_ref * s_psi);
            x = (c_chi_ref * c_psi) + (s_chi_ref * s_psi);
            
            e_chi = atan2(y, x);
        end

        function [X0, U0, ctrl_lon, ctrl_lat] = scheduled_params(obj, u_cmd, w_cmd)
            [X0, U0] = obj.interp_xu0(u_cmd, w_cmd);

            ctrl_lon.Ki  = obj.interp_mtrx(obj.LON.Ki, u_cmd, w_cmd);
            ctrl_lon.Kx  = obj.interp_mtrx(obj.LON.Kx, u_cmd, w_cmd);
            ctrl_lon.Kv  = obj.interp_mtrx(obj.LON.Kv, u_cmd, w_cmd);
            ctrl_lon.F   = obj.interp_mtrx( obj.LON.F, u_cmd, w_cmd);
            ctrl_lon.G   = obj.interp_mtrx( obj.LON.G, u_cmd, w_cmd);
            ctrl_lon.C   = obj.interp_mtrx( obj.LON.C, u_cmd, w_cmd);
            ctrl_lon.Cv  = obj.interp_mtrx(obj.LON.Cv, u_cmd, w_cmd);
            ctrl_lon.W   = obj.interp_mtrx( obj.LON.W, u_cmd, w_cmd);
            ctrl_lon.B   = obj.interp_mtrx( obj.LON.B, u_cmd, w_cmd);
            ctrl_lon.Ap  = obj.interp_mtrx(obj.LON.Ap, u_cmd, w_cmd);
            ctrl_lon.Bp  = obj.interp_mtrx(obj.LON.Bp, u_cmd, w_cmd);

            ctrl_lat.Ki  = obj.interp_mtrx(obj.LAT.Ki, u_cmd, w_cmd);
            ctrl_lat.Kx  = obj.interp_mtrx(obj.LAT.Kx, u_cmd, w_cmd);
            ctrl_lat.Kv  = obj.interp_mtrx(obj.LAT.Kv, u_cmd, w_cmd);
            ctrl_lat.F   = obj.interp_mtrx( obj.LAT.F, u_cmd, w_cmd);
            ctrl_lat.G   = obj.interp_mtrx( obj.LAT.G, u_cmd, w_cmd);
            ctrl_lat.C   = obj.interp_mtrx( obj.LAT.C, u_cmd, w_cmd);
            ctrl_lat.Cv  = obj.interp_mtrx(obj.LAT.Cv, u_cmd, w_cmd);
            ctrl_lat.W   = obj.interp_mtrx( obj.LAT.W, u_cmd, w_cmd);
            ctrl_lat.B   = obj.interp_mtrx( obj.LAT.B, u_cmd, w_cmd);
            ctrl_lat.Ap  = obj.interp_mtrx(obj.LAT.Ap, u_cmd, w_cmd);
            ctrl_lat.Bp  = obj.interp_mtrx(obj.LAT.Bp, u_cmd, w_cmd);
        end

        function [Xlon, Xlat, Xlon_cmd, Xlat_cmd] = perturb_lin_ctrl(obj, state, ref, e_pos_body, e_chi, X0)
            ref_vel = ref.vel;
            ref_chi_dot = ref.chi_dot;

            omg_s = state(10:12);
            omg_0 = X0(4:6);
            omg_t = omg_s - omg_0;

            Vb_s = state(4:6);
            Vb_0 = X0(1:3);
            Vb_t = Vb_s - Vb_0;

            switch obj.w_ff_mode
                case 0,  w_ff = 0;
                case 1,  w_ff = ref_vel(1)*tan(X0(11)) - X0(3);
                case 2,  w_ff = state(4) *tan(state(8)) - X0(3);
            end
            Vb_cmd_t = ref_vel - Vb_0 + e_pos_body.*obj.k_pos + [0; 0; w_ff];
            
            chi_dot_cmd = ref_chi_dot + e_chi * 0.1;

            eta_s = state(7:9);
            eta_0 = X0(10:12);
            eta_t = eta_s - eta_0;

            Xlon = [Vb_t(1); Vb_t(3); omg_t(2); eta_t(2)];
            Xlat = [Vb_t(2); omg_t(1); omg_t(3); eta_t(1)];
            Xlon_cmd = [Vb_cmd_t(1); Vb_cmd_t(3); 0];
            Xlat_cmd = [Vb_cmd_t(2); chi_dot_cmd];
        end
    
        function mdes_lon = lon_ctrl(obj, ctrl_lon, Xlon_cmd, Xlon)
            obj.ctrl_lon_i = obj.ctrl_lon_i +...
             (ctrl_lon.F * Xlon_cmd - ctrl_lon.C * Xlon + ctrl_lon.Kv * (obj.theta_v - ctrl_lon.Cv * Xlon)) .* obj.dt;
            
            mdes_lon = ctrl_lon.G * Xlon_cmd + ctrl_lon.Ki * obj.ctrl_lon_i - ctrl_lon.Kx * Xlon;
        end
        
        function mdes_lat = lat_ctrl(obj, ctrl_lat, Xlat_cmd, Xlat)
            obj.ctrl_lat_i = obj.ctrl_lat_i +...
             (ctrl_lat.F * Xlat_cmd - ctrl_lat.C * Xlat + ctrl_lat.Kv * (obj.phi_v - ctrl_lat.Cv * Xlat)) .* obj.dt;
            
            mdes_lat = ctrl_lat.G * Xlat_cmd + ctrl_lat.Ki * obj.ctrl_lat_i - ctrl_lat.Kx * Xlat;
        end

        function eng_cmd = engine_cmd(~, u_lon, u_lat, U0)
            lift_cmd = u_lat(1:end-2) + u_lon(1:end-3);
            prop_cmd = u_lon(9);
            
            eng_cmd = [lift_cmd; prop_cmd] + U0(5:end);
        end
        
        function srf_cmd = srface_cmd(~, u_lon, u_lat, U0)
            ele_cmd = u_lon(10);
            flp_cmd = u_lon(11);
            ail_cmd = u_lat(9);
            rud_cmd = u_lat(10);
            srf_0 = [U0(1) - U0(2);...
                     U0(1) + U0(2);...
                     U0(3);...
                     U0(3);...
                     U0(4)];
                     
            srf_cmd = [flp_cmd - ail_cmd;...
                       flp_cmd + ail_cmd;...
                       ele_cmd;...
                       ele_cmd;...
                       rud_cmd] + srf_0;
        end

        function perturb_cmd = engine_surface_perturb(~, u_lon, u_lat)
            lift_cmd = u_lat(1:end-2) + u_lon(1:end-3);
            prop_cmd = u_lon(9);
            
            ele_cmd = u_lon(10);
            flp_cmd = u_lon(11);
            ail_cmd = u_lat(9);
            rud_cmd = u_lat(10);

            perturb_cmd = [lift_cmd; prop_cmd; flp_cmd; ail_cmd; ele_cmd; rud_cmd];
        end

        function [engine, surface] = total_cmd(~, perturb_cmd, U0)
            surface_0 = [U0(1) - U0(2);...
                         U0(1) + U0(2);...
                         U0(3);...
                         U0(3);...
                         U0(4)];
            engine_0 = U0(5:end);

            engine = perturb_cmd(1:9) + engine_0;
            flp_p = perturb_cmd(10);
            ail_p = perturb_cmd(11);
            ele_p = perturb_cmd(12);
            rud_p = perturb_cmd(13);

            surface = [flp_p - ail_p;...
                       flp_p + ail_p;...
                       ele_p;...
                       ele_p;...
                       rud_p] + surface_0;
        end

        function perturb_cmd = pseudo_alloc(obj, ctrl_lon, ctrl_lat, mdes_lon, mdes_lat)
            % Weighted pseudo-inverse allocation: M = W^-1 B' (B W^-1 B')^-1
            % Produces the nominal effector perturbation (pre-liveness-filter).
            M_lon = (ctrl_lon.W \ ctrl_lon.B') / (ctrl_lon.B * (ctrl_lon.W \ ctrl_lon.B'));
            M_lat = (ctrl_lat.W \ ctrl_lat.B') / (ctrl_lat.B * (ctrl_lat.W \ ctrl_lat.B'));

            act_lon = M_lon*mdes_lon;

            scale_factor = 1; % Initialize scale factor to 1

            % Check if any engines exceed max position limit
            % (matches vehicles/Lift+Cruise/Control/Limited_Long_Cont_Alloc.m)
            if any(act_lon(1:9) > obj.rslqrCfg.eng_max)
                for loop = 1:9
                    scale_temp = obj.rslqrCfg.eng_max(loop)/(M_lon(loop,:)*mdes_lon);
                    if scale_temp>0 && scale_temp<1 && scale_temp < scale_factor
                        scale_factor = scale_temp;
                    end
                end
            end
            % Check if elevator exceed limits
            if act_lon(10) > obj.rslqrCfg.ele_max
                scale_temp =  obj.rslqrCfg.ele_max/(M_lon(10,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            if act_lon(10) < obj.rslqrCfg.ele_min
                scale_temp = obj.rslqrCfg.ele_min/(M_lon(10,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            % Check if elevator or flaps exceed limits
            if act_lon(11) > obj.rslqrCfg.flp_max
                scale_temp = obj.rslqrCfg.flp_max/(M_lon(11,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            if act_lon(11) < obj.rslqrCfg.flp_min
                scale_temp = obj.rslqrCfg.flp_min/(M_lon(11,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end

            % Scale allocation if required
            if scale_factor ~= 1
                act_lon = M_lon*mdes_lon*scale_factor;
            end

            % Parse out the effector allocations. The virtual attitude
            % effectors (theta_v/phi_v) feed the next step's integral
            % compensator (lon_ctrl/lat_ctrl).
            u_lon       = act_lon(1:end-1);
            obj.theta_v = act_lon(end);

            act_lat = M_lat * mdes_lat;
            u_lat   = act_lat(1:end-1);
            obj.phi_v = act_lat(end);

            % Nominal effector perturbation (13x1):
            % [lift1..8; pusher; flp; ail; ele; rud]. The longitudinal liveness
            % filter and final command assembly are applied by Controller.
            perturb_cmd = obj.engine_surface_perturb(u_lon, u_lat);
        end

        function [perturb_cmd, U0] = control(obj, state, ref)
            % Baseline RSLQR pass: nominal effector perturbation + scheduled
            % trim input U0. The safety filter and absolute-command assembly
            % are handled by the owning Controller.
            % ref : struct with fields pos (3x1), vel (3x1, heading frame),
            %       chi (scalar), chi_dot (scalar)
            u_cmd = ref.vel(1);
            w_cmd = ref.vel(3);

            [e_pos_body, e_chi] = obj.guidance_error(state, ref);
            [X0, U0, ctrl_lon, ctrl_lat] = obj.scheduled_params(u_cmd, w_cmd);
            [Xlon, Xlat, Xlon_cmd, Xlat_cmd] = obj.perturb_lin_ctrl(state, ref, e_pos_body, e_chi, X0);
            mdes_lon = obj.lon_ctrl(ctrl_lon, Xlon_cmd, Xlon);
            mdes_lat = obj.lat_ctrl(ctrl_lat, Xlat_cmd, Xlat);

            perturb_cmd = obj.pseudo_alloc(ctrl_lon, ctrl_lat, mdes_lon, mdes_lat);
        end

        function reset(obj)
            obj.ctrl_lon_i = zeros(3, 1);
            obj.ctrl_lat_i = zeros(3, 1);
            obj.phi_v   = 0;
            obj.theta_v = 0;
        end

        function S = set_gains(obj)
            % Regenerate the gain-scheduled lon/lat controller tables offline
            % from the stored trim table obj.XU0 (25 x N_trim x M_trim).
            %
            % Port of ctrl_scheduler_GUAM.m (NASA GUAM, Jacob Cook / M. J. Acheson).
            cfg = obj.rslqrCfg;
            N   = cfg.N_trim;
            M   = cfg.M_trim;

            LON_ = [];
            LAT_ = [];
            % Fill order (ii fastest, then jj) must match the reshape below:
            % LON_ is 1 x (N*M) with the UH index varying fastest.
            for jj = 1 : M
                for ii = 1 : N
                    trim_point = obj.trim_xu_eq(obj.XU0(:, ii, jj));

                    [lon, lon_err] = obj.ctrl_lon(trim_point, cfg.Qlon0, cfg.Rlon0, cfg.Wlon0);
                    if lon_err
                        error('RSLQR:ctrl_lon', 'Error in ctrl_lon at UH=%g, WH=%g', trim_point(1), trim_point(2));
                    end
                    [lat, lat_err] = obj.ctrl_lat(trim_point, cfg.Qlat0, cfg.Rlat0, cfg.Wlat0);
                    if lat_err
                        error('RSLQR:ctrl_lat', 'Error in ctrl_lat at UH=%g, WH=%g', trim_point(1), trim_point(2));
                    end

                    LON_ = [LON_, lon]; %#ok<AGROW>
                    LAT_ = [LAT_, lat]; %#ok<AGROW>
                end
            end

            % Longitudinal System
            S.Ki_lon_interp = reshape([LON_.Ki], cfg.Ni_lon, cfg.Ni_lon, N, M);
            S.Kx_lon_interp = reshape([LON_.Kx], cfg.Ni_lon, cfg.Nx_lon, N, M);
            S.Kv_lon_interp = reshape([LON_.Kv], cfg.Ni_lon, cfg.Nv_lon, N, M);
            S.F_lon_interp  = reshape([LON_.F],  cfg.Ni_lon, cfg.Nr_lon, N, M);
            S.G_lon_interp  = reshape([LON_.G],  cfg.Ni_lon, cfg.Nr_lon, N, M);
            S.C_lon_interp  = reshape([LON_.C],  cfg.Ni_lon, cfg.Nx_lon, N, M);
            S.Cv_lon_interp = reshape([LON_.Cv], cfg.Nv_lon, cfg.Nx_lon, N, M);
            S.W_lon_interp  = reshape([LON_.W],  cfg.Nu_lon+cfg.Nv_lon, cfg.Nu_lon+cfg.Nv_lon, N, M);
            S.B_lon_interp  = reshape([LON_.B],  cfg.Ni_lon, cfg.Nu_lon+cfg.Nv_lon, N, M);
            S.Ap_lon_interp = reshape([LON_.Ap], cfg.Nx_lon, cfg.Nx_lon, N, M);
            S.Bp_lon_interp = reshape([LON_.Bp], cfg.Nx_lon, cfg.Nu_lon, N, M);

            % Lateral System
            S.Ki_lat_interp = reshape([LAT_.Ki], cfg.Ni_lat, cfg.Ni_lat, N, M);
            S.Kx_lat_interp = reshape([LAT_.Kx], cfg.Ni_lat, cfg.Nx_lat, N, M);
            S.Kv_lat_interp = reshape([LAT_.Kv], cfg.Ni_lat, cfg.Nv_lat, N, M);
            S.F_lat_interp  = reshape([LAT_.F],  cfg.Ni_lat, cfg.Nr_lat, N, M);
            S.G_lat_interp  = reshape([LAT_.G],  cfg.Ni_lat, cfg.Nr_lat, N, M);
            S.C_lat_interp  = reshape([LAT_.C],  cfg.Ni_lat, cfg.Nx_lat, N, M);
            S.Cv_lat_interp = reshape([LAT_.Cv], cfg.Nv_lat, cfg.Nx_lat, N, M);
            S.W_lat_interp  = reshape([LAT_.W],  cfg.Nu_lat+cfg.Nv_lat, cfg.Nu_lat+cfg.Nv_lat, N, M);
            S.B_lat_interp  = reshape([LAT_.B],  cfg.Ni_lat, cfg.Nu_lat+cfg.Nv_lat, N, M);
            S.Ap_lat_interp = reshape([LAT_.Ap], cfg.Nx_lat, cfg.Nx_lat, N, M);
            S.Bp_lat_interp = reshape([LAT_.Bp], cfg.Nx_lat, cfg.Nu_lat, N, M);
        end

        function [out, ctrl_error] = ctrl_lon(obj, xu_eq, q, r, wc)
            % Design a unified longitudinal controller at one trim point.
            % Port of ctrl_lon.m (NASA GUAM, Jacob Cook / M. J. Acheson).
            %   xu_eq : [ubar; wbar; Radius; th; phi; p; q; r; delf; dela; ...
            %            dele; delr; omp1..omp9]  (21x1, control-frame trim)
            %   q,r,wc: LQR state / acceleration / allocation weight vectors
            ctrl_error = 0;

            Q  = diag(q);
            R  = diag(r);
            Wc = diag(wc);

            xeq = xu_eq(1:8);    % trim states
            ueq = xu_eq(9:end);  % trim effectors

            NS = 4;  % aero control surfaces
            NP = 9;  % rotor/prop effectors

            [Alon, Blon, Clon, Dlon] = obj.get_long_dynamics_heading(xeq, ueq, NS, NP);

            if xeq(1) > obj.rslqrCfg.trans_end % cruise: zero the lifting rotors
                Blon(:,1:8) = zeros(4,8);
            end

            % size definitions
            Ni  = 3;  % integrator states
            Nr  = 3;  % reference
            Nu  = 11; % physical controls
            Nv  = 1;  % virtual controls
            Nxi = 3;  % general states

            % Performance design with general acceleration inputs
            Av = Alon([1 2 3], [1 2 3]);
            Cv = eye(Nxi);

            At = [ zeros(Ni,Ni)   Cv ;
                   zeros(Nxi,Ni)  Av ];
            Bt = [ zeros(Nxi); eye(Nxi) ];

            % LQR optimal feedback gains
            Kc  = lqr(At, Bt, Q, R);
            Ki0 = Kc(:,1:Ni);
            Kx0 = Kc(:,Ni+1:Ni+Nxi);

            % Control allocation design: Bu = [u u_v]
            Bu = [ Blon([1 2 3],:) Alon([1 2 3], 4) ];
            M  = Wc\Bu'/(Bu/(Wc)*Bu');

            A = Alon;
            B = Blon;
            C = Clon([1 2 3],:); % tracking states

            % add a column of zeros so theta comes back in as a state
            Kx = Kx0*[1 0 0 0; 0 1 0 0; 0 0 1 0];
            Ki = Ki0;

            % Break up the allocation matrix
            Mu = M(1:Nu,:);
            Mv = M(Nu+Nv,:);

            % Virtual control
            Cv = Clon(4,:);
            Kv = [0; 0; 1];

            F = eye(3);
            G = zeros(3);

            Acl = [ Kv*Mv*Ki    Kv*Mv*Kx+Kv*Cv+C ;
                    -B*Mu*Ki           A-B*Mu*Kx ];
            if sum(eig(Acl) > 0)
                fprintf('trim point has unstable poles\n');
                ctrl_error = 1;
            end

            % Assign outputs (only the fields consumed at runtime are kept)
            out.Ap = A;
            out.Bp = B;
            out.Ki = Ki;
            out.Kx = Kx;
            out.Kv = Kv;
            out.F  = F;
            out.G  = G;
            out.C  = C;
            out.Cv = Cv;
            out.W  = Wc;
            out.B  = Bu;
        end

        function [out, ctrl_error] = ctrl_lat(obj, xu_eq, q, r, wc)
            % Design a unified lateral controller at one trim point.
            % Port of ctrl_lat.m (NASA GUAM, Jacob Cook / M. J. Acheson).
            %   xu_eq : same 21x1 control-frame trim vector as ctrl_lon
            %   q,r,wc: LQR state / acceleration / allocation weight vectors
            ctrl_error = 0;

            Q  = diag(q);
            R  = diag(r);
            Wc = diag(wc);

            xeq = xu_eq(1:8);    % trim states
            ueq = xu_eq(9:end);  % trim effectors

            NS = 4;  % aero control surfaces
            NP = 9;  % rotor/prop effectors

            [Alat, Blat, Clat, Dlat, XU0] = obj.get_lat_dynamics_heading(xeq, ueq, NS, NP);

            if xeq(1) > obj.rslqrCfg.trans_end % cruise: zero the lifting rotors
                Blat(:,1:8) = zeros(4,8);
            end

            % size definitions
            Ni  = 3;  % integrator states
            Nu  = 10; % physical controls
            Nv  = 1;  % virtual controls
            Nxi = 3;  % general states

            % Performance design with general acceleration inputs
            Av = Alat([1 2 3], [1 2 3]);
            Cv = eye(Nxi);

            At = [ zeros(Ni,Ni)   Cv ;
                   zeros(Nxi,Ni)  Av ];
            Bt = [ zeros(Nxi); eye(Nxi) ];

            % LQR optimal feedback gains
            Kc  = lqr(At, Bt, Q, R);
            Ki0 = Kc(:,1:Ni);
            Kx0 = Kc(:,Ni+1:Ni+Nxi);

            % Control allocation design: Bu = [u u_v]
            Bu = [ Blat([1 2 3],:) Alat([1 2 3], 4) ];
            M  = Wc\Bu'/(Bu/(Wc)*Bu');

            A = Alat;
            B = Blat;
            C = Clat([1 2 3],:); % tracking states

            % add a column of zeros so phi comes back in as a state
            Kx = Kx0*[1 0 0 0; 0 1 0 0; 0 0 1 0];
            Ki = Ki0;

            % Break up the allocation matrix
            Mu = M(1:Nu,:);
            Mv = M(Nu+Nv,:);

            % Virtual control
            Cv = Clat(4,:);
            Kv = [0; 1; 0];

            th0  = XU0(11);
            phi0 = XU0(10);
            q0   = XU0(5);
            vc0  = XU0(1:3);
            vb0  = Rx(phi0)*Ry(th0)*vc0;
            ub0  = vb0(1);

            F = [1 0; 0 0; 0 1-phi0*q0];
            G = [0 ub0; 0 0; 0 0];

            Acl = [ Kv*Mv*Ki    Kv*Mv*Kx+Kv*Cv+C ;
                    -B*Mu*Ki           A-B*Mu*Kx ];
            if sum(eig(Acl) > 0)
                fprintf('trim point has unstable poles\n');
                ctrl_error = 1;
            end

            % Assign outputs (only the fields consumed at runtime are kept)
            out.Ap = A;
            out.Bp = B;
            out.Ki = Ki;
            out.Kx = Kx;
            out.Kv = Kv;
            out.F  = F;
            out.G  = G;
            out.C  = C;
            out.Cv = Cv;
            out.W  = Wc;
            out.B  = Bu;
        end

        function [Alon, Blon, Clon, Dlon, XU0] = get_long_dynamics_heading(obj, XEQ, UEQ, NS, NP)
            % Longitudinal slice of the full heading-frame linear dynamics.
            % Port of get_long_dynamics_heading.m (polynomial branch).
            %   x = [ ubar wbar q th ]',  u = [ omp1-9 dele delf ]'
            [A, B, C, D, XU0] = obj.get_lin_dynamics_heading(XEQ, UEQ, NS, NP);

            % get_lin_dynamics states/controls:
            %   x = [ px py pz phi th psi ubar vbar wbar p q r ]'
            %   u = [ delf dela dele delr omp1-9 ]'
            x_idx = [7 9 11 5];  % [ubar wbar q theta]
            u_idx = [5:13 3 1];  % [omp1-9 elev flap] (GUAM simulation order)
            y_idx = x_idx;

            Alon = A(x_idx,x_idx);
            Blon = B(x_idx,u_idx);
            Clon = C(y_idx,x_idx);
            Dlon = D(y_idx,u_idx);
        end

        function [Alat, Blat, Clat, Dlat, XU0] = get_lat_dynamics_heading(obj, XEQ, UEQ, NS, NP)
            % Lateral slice of the full heading-frame linear dynamics.
            % Port of get_lat_dynamics_heading.m (polynomial branch).
            %   x = [ vbar p r phi ]',  u = [ omp1-8 dela delr ]'
            [A, B, C, D, XU0] = obj.get_lin_dynamics_heading(XEQ, UEQ, NS, NP);

            x_idx = [8 10 12 4]; % [vbar p r phi]
            u_idx = [5:12 2 4];  % [omp1-8 ail rud] (GUAM simulation order)
            y_idx = x_idx;

            Alat = A(x_idx,x_idx);
            Blat = B(x_idx,u_idx);
            Clat = C(y_idx,x_idx);
            Dlat = D(y_idx,u_idx);
        end

        function [A, B, C, D, XU0] = get_lin_dynamics_heading(obj, XEQ, UEQ, NS, NP)
            % Full (12-state) linearized heading-frame dynamics, obtained by
            % numerically differentiating the polynomial aero/propulsive model
            % about a trim condition. Port of get_lin_dynamics_heading.m
            % (polynomial-database branch; strip-theory branch dropped).
            %
            %  dx = A*x + B*u,   y = C*x + D*u
            %   x = [ px py pz  phi th psi  ubar vbar wbar  p q r ]'
            %   u = [ delf dela dele delr  omp1..omp9 ]'
            %   y = x
            %
            %   XEQ : [ubar; wbar; Radius; th; phi; p; q; r]
            %   UEQ : [delf; dela; dele; delr; omp1..omp9]
            %   XU0 : trim vector [vbar; om; ab; eta; del; omp] (25x1)

            rho = obj.rslqrCfg.rho;
            g   = obj.rslqrCfg.grav;

            % unit vectors, rotations and skew-symmetric helper
            e1 = [1; 0; 0];
            e3 = [0; 0; 1];

            Rx  = @(x)  [1 0 0 ; 0 cos(x) sin(x); 0 -sin(x) cos(x)];
            Ry  = @(x)  [cos(x) 0 -sin(x); 0 1 0; sin(x) 0 cos(x)];
            Rz  = @(x)  [cos(x) sin(x) 0; -sin(x) cos(x) 0; 0 0 1];
            hat = @(x)  [ 0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

            % equilibrium state
            ubar   = XEQ(1);
            wbar   = XEQ(2);
            Radius = XEQ(3);
            th     = XEQ(4);
            phi    = XEQ(5);
            psi    = 0;
            p      = XEQ(6);
            q      = XEQ(7);
            r      = XEQ(8);

            % equilibrium controls
            del = UEQ(1:NS);
            omp = UEQ(NS+1:NS+NP);

            % convert Radius to turn rate (Inf => straight flight => dpsi = 0)
            dpsi = ubar/Radius;

            % body accelerations
            ab = [ -g*sin(th); g*cos(th)*sin(phi); g*cos(th)*cos(phi) ];

            eta  = [ phi; th; psi];
            vbar = [ubar; 0; wbar];
            vb   = Rx(phi)*Ry(th)*vbar;
            om   = [ p; q; r ];

            % save the trim state
            XU0 = [vbar; om; ab; eta; del; omp];

            % Interrogate the polynomial database for forces/moments and their
            % state/input derivatives.
            [F, ~, F_x, M_x, F_u, M_u] = poly_aero_wrapper_Mod_du([vb; om]', [del; omp]', rho);

            F_vb = F_x(1:3,1:3);
            M_vb = M_x(1:3,1:3);
            F_om = F_x(1:3,4:6);
            M_om = M_x(1:3,4:6);
            F_u  = F_u(1:3,:);
            M_u  = M_u(1:3,:);

            % mass and inertia
            m = obj.vehicleCfg.mass;
            J = obj.vehicleCfg.I;

            S     = [1   sin(phi)*tan(th)   cos(phi)*tan(th) ;
                     0   cos(phi)          -sin(phi)         ;
                     0   sin(phi)/cos(th)   cos(phi)/cos(th)];

            S_phi = [0   cos(phi)*tan(th)  -sin(phi)*tan(th) ;
                     0  -sin(phi)          -cos(phi)         ;
                     0   cos(phi)/cos(th)  -sin(phi)/cos(th)];

            S_th  = [0   sin(phi)*sec(th)^2        cos(phi)*sec(th)^2       ;
                     0   0                         0                        ;
                     0   sin(phi)*tan(th)/cos(th)  cos(phi)*tan(th)/cos(th)];

            Som_eta = [S_phi*om S_th*om zeros(3,1)];

            % intermediate derivatives
            Rbar     =  Rx(phi)*Ry(th);
            Rbar_phi = -Rx(phi)*hat(e1)*Ry(th);
            Rbar_th  = -Rx(phi)*Ry(th)*hat([0;1;0]);
            Rbar_psi =  zeros(3,3);

            vb_vbar = Rbar;
            vb_eta  = [Rbar_phi*vbar Rbar_th*vbar Rbar_psi*vbar];

            % force/moment derivatives w.r.t. states and inputs
            Fbar_vbar = Rbar'*F_vb*vb_vbar;
            Fbar_om   = Rbar'*F_om;
            Fbar_eta  = [Rbar_phi'*F Rbar_th'*F Rbar_psi'*F] + Rbar'*F_vb*vb_eta;
            Fbar_u    = Rbar'*F_u;

            M_vbar = M_vb*vb_vbar;
            M_eta  = M_vb*vb_eta;

            % state-space matrices
            A = [zeros(3,3)  Rz(psi)*hat(e3)  Rz(psi)                            zeros(3,3)                    ;
                 zeros(3,3)  Som_eta          zeros(3,3)                         S                             ;
                 zeros(3,3)  1/m*Fbar_eta    -hat([0; 0; dpsi])+1/m*Fbar_vbar    1/m*Fbar_om                   ;
                 zeros(3,3)  J\M_eta          J\M_vbar                           J\(-hat(om)*J+hat(J*om)+M_om) ];

            B = [ zeros(3,NS+NP) ;
                  zeros(3,NS+NP) ;
                  1/m*Fbar_u     ;
                  J\M_u          ];

            C = eye(12);
            D = zeros(12,NS+NP);
        end

    end % end of methods

    methods (Static)
        function xu_eq = trim_xu_eq(xu0)
            % Reconstruct the 21x1 control-frame trim vector expected by
            % ctrl_lon/ctrl_lat from a 25x1 XU0_interp column
            % [vbar(=[ubar;0;wbar]); om(=[p;q;r]); ab; eta(=[phi;th;psi]);
            %  del(=[flap;ail;ele;rud]); omp(1:9)].
            
            ubar = xu0(1);
            wbar = xu0(3);
            phi  = xu0(10);
            th   = xu0(11);
            p    = xu0(4);
            q    = xu0(5);
            r    = xu0(6);
            del  = xu0(13:16);
            omp  = xu0(17:25);

            xu_eq = [ubar; wbar; Inf; th; phi; p; q; r; del; omp];
        end

        function rotm = rotm_i2b(phi, theta, psi)
            cphi = cos(phi);
            sphi = sin(phi);
            rotm_phi = [1,     0,    0;...
                        0,  cphi, sphi;...
                        0, -sphi, cphi];

            ctheta = cos(theta);
            stheta = sin(theta);
            rotm_theta = [ctheta, 0, -stheta;...
                               0, 1,       0;...
                          stheta, 0,  ctheta];

            cpsi = cos(psi);
            spsi = sin(psi);
            rotm_psi = [cpsi, spsi, 0;...
                       -spsi, cpsi, 0;...
                           0,    0, 1];

            rotm = rotm_phi * rotm_theta * rotm_psi;
        end

        function rotm = rotm_b2i(phi, theta, psi)
            rotm = RSLQR.rotm_i2b(phi, theta, psi)';
        end
    end


end