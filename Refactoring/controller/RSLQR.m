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
    %   RSLQRConfig, get_long_dynamics_heading, get_lat_dynamics_heading, Rx, Ry

    properties
        cfg       % RSLQRConfig (weights and dimension constants)
        UH        % N_trim x 1  body-frame x-velocity (u) breakpoints [ft/s]
        WH        % M_trim x 1  body-frame z-velocity (w) breakpoints [ft/s]
                  % (verified: WH == XU0(3)=body w exactly; NOT heading/inertial.
                  %  At cruise pitch, body w>0 is ~level flight, not descent.)
        XU0
        LON       % 1x(N*M*L) struct array — longitudinal gains per schedule point
        LAT       % 1x(N*M*L) struct array — lateral gains per schedule point

        ctrl_lon_i  % Longitudinal integral term of controller
        ctrl_lat_i  % Lateral integral term of controller

        phi_v
        theta_v

        liveness_lon   % LivenessFilter (longitudinal) behind nominal allocation
        liveness_wh_anchor = []
    end

    methods

        function obj = RSLQR()
            % Load the precomputed gain-scheduled RSLQR design
            % (trim table + lon/lat gain tables over the UH x WH envelope).
            % trim_table_Poly_ConcatVer4p0.mat must be on the MATLAB path.

            obj.cfg = RSLQRConfig;
            S = load('tables/trim/trim_table_Poly_ConcatVer4p0.mat');
            obj.XU0 = S.XU0_interp;
            obj.UH  = S.UH;
            obj.WH  = S.WH;

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

            % Longitudinal liveness filter (production channel). Loads BRT
            % value functions from LivenessConfig.tables_dir_default; passes
            % through when tables are absent or (uh,wh) is outside coverage.
            % Default mode 'blend' (smooth blending). Callers may override
            % obj.liveness_lon.mode ('blend' | 'lr' | 'off').
            obj.liveness_lon = LivenessFilter('lon', 'blend', ...
                LivenessConfig.tables_dir_default, obj.UH, obj.WH);
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

        function [Xlon, Xlat, Xlon_cmd, Xlat_cmd] = perturb_lin_ctrl(~, state, ref, e_pos_body, e_chi, X0)
            ref_vel = ref.vel;
            ref_chi_dot = ref.chi_dot;

            omg_s = state(10:12);
            omg_0 = X0(4:6);
            omg_t = omg_s - omg_0;

            Vb_s = state(4:6);
            Vb_0 = X0(1:3);
            Vb_t = Vb_s - Vb_0;

            Vb_cmd_t = ref_vel - Vb_0 + e_pos_body.*0.1;
            
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
             (ctrl_lon.F * Xlon_cmd - ctrl_lon.C * Xlon + ctrl_lon.Kv * (obj.theta_v - ctrl_lon.Cv * Xlon)) .* obj.cfg.dt;
            
            mdes_lon = ctrl_lon.G * Xlon_cmd + ctrl_lon.Ki * obj.ctrl_lon_i - ctrl_lon.Kx * Xlon;
        end
        
        function mdes_lat = lat_ctrl(obj, ctrl_lat, Xlat_cmd, Xlat)
            obj.ctrl_lat_i = obj.ctrl_lat_i +...
             (ctrl_lat.F * Xlat_cmd - ctrl_lat.C * Xlat + ctrl_lat.Kv * (obj.phi_v - ctrl_lat.Cv * Xlat)) .* obj.cfg.dt;
            
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

        function [engine, surface] = pseudo_alloc(obj, U0, ctrl_lon, ctrl_lat, mdes_lon, mdes_lat, Xlon, uh, wh, state)
            %% Weighted pseudo-inverse allocation: M = W^-1 B' (B W^-1 B')^-1
            M_lon = (ctrl_lon.W \ ctrl_lon.B') / (ctrl_lon.B * (ctrl_lon.W \ ctrl_lon.B'));
            M_lat = (ctrl_lat.W \ ctrl_lat.B') / (ctrl_lat.B * (ctrl_lat.W \ ctrl_lat.B'));

            act_lon = M_lon*mdes_lon;

            scale_factor = 1; % Initialize scale factor to 1

            % Check if any engines exceed max position limit
            % (matches vehicles/Lift+Cruise/Control/Limited_Long_Cont_Alloc.m)
            if any(act_lon(1:9) > obj.cfg.eng_max)
                for loop = 1:9
                    scale_temp = obj.cfg.eng_max(loop)/(M_lon(loop,:)*mdes_lon);
                    if scale_temp>0 && scale_temp<1 && scale_temp < scale_factor
                        scale_factor = scale_temp;
                    end
                end
            end
            % Check if elevator exceed limits
            if act_lon(10) > obj.cfg.ele_max
                scale_temp =  obj.cfg.ele_max/(M_lon(10,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            if act_lon(10) < obj.cfg.ele_min
                scale_temp = obj.cfg.ele_min/(M_lon(10,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            % Check if elevator or flaps exceed limits
            if act_lon(11) > obj.cfg.flp_max
                scale_temp = obj.cfg.flp_max/(M_lon(11,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end
            if act_lon(11) < obj.cfg.flp_min
                scale_temp = obj.cfg.flp_min/(M_lon(11,:)*mdes_lon);
                if scale_temp>0 && scale_temp<1 && scale_temp<scale_factor
                    scale_factor = scale_temp;
                end
            end

            % Scale allocation if required
            if scale_factor ~= 1
                act_lon = M_lon*mdes_lon*scale_factor;
            end

            %% Longitudinal liveness filter (behind nominal allocation).
            u_lon_nom = act_lon(1:11);  % nominal allocation before liveness filter
            uhA = state(4);
            uhA = max(obj.UH(1), min(obj.UH(end), uhA));
            whA        = obj.liveness_wh_anchor;
            [X0a, U0a] = obj.interp_xu0(uhA, whA); % Anchor trim state and input
            Ap_a       = obj.interp_mtrx(obj.LON.Ap, uhA, whA); % Linear dynamics at anchor trim
            Bp_a       = obj.interp_mtrx(obj.LON.Bp, uhA, whA); % Linear dynamics at anchor trim
            % perturbation state relative to the anchor trim [u; w; q; theta]
            x_anchor   = [state(4)  - X0a(1); ...
                          state(6)  - X0a(3); ...
                          state(11) - X0a(5); ...
                          state(8)  - X0a(11)];
            idx        = obj.liveness_lon.spec.U0_idx;
            dtrim      = U0(idx) - U0a(idx);          % 11x1 effector trim offset
            u_anchor_nom  = u_lon_nom + dtrim;        % mission frame -> anchor frame
            u_anchor_f    = obj.liveness_lon.filter(x_anchor, u_anchor_nom, ...
                                Ap_a, Bp_a, U0a, uhA, whA, state);
            act_lon(1:11) = u_anchor_f - dtrim;       % anchor frame -> mission frame

            % Parse out the effector allocations
            u_lon       = act_lon(1:end-1);
            obj.theta_v = act_lon(end);

            act_lat = M_lat * mdes_lat;
            u_lat = act_lat(1:end-1);
            obj.phi_v = act_lat(end);

            engine  = obj.engine_cmd(u_lon, u_lat, U0);
            surface = obj.srface_cmd(u_lon, u_lat, U0);
        end

        function [engine, surface] = control(obj, state, ref)
            % ref : struct with fields pos (3x1), vel (3x1, heading frame),
            %       chi (scalar), chi_dot (scalar)
            u_cmd = ref.vel(1);
            w_cmd = ref.vel(3);

            [e_pos_body, e_chi] = obj.guidance_error(state, ref);
            [X0, U0, ctrl_lon, ctrl_lat] = obj.scheduled_params(u_cmd, w_cmd);
            [Xlon, Xlat, Xlon_cmd, Xlat_cmd] = obj.perturb_lin_ctrl(state, ref, e_pos_body, e_chi, X0);
            mdes_lon = obj.lon_ctrl(ctrl_lon, Xlon_cmd, Xlon);
            mdes_lat = obj.lat_ctrl(ctrl_lat, Xlat_cmd, Xlat);

            [engine, surface] = obj.pseudo_alloc(U0, ctrl_lon, ctrl_lat, mdes_lon, mdes_lat, Xlon, u_cmd, w_cmd, state);
        end

        function set_liveness_dynamics(obj, fdyn)
            if ~isempty(obj.liveness_lon)
                obj.liveness_lon.fdyn = fdyn;
            end
        end

        function reset(obj)
            obj.ctrl_lon_i = zeros(3, 1);
            obj.ctrl_lat_i = zeros(3, 1);
            obj.phi_v   = 0;
            obj.theta_v = 0;
            if ~isempty(obj.liveness_lon)
                obj.liveness_lon.reset_counters();
            end
        end
    end % methods

    methods (Static)
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
