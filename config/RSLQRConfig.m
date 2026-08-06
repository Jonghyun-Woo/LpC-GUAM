classdef RSLQRConfig
    
    properties (Constant)
        L_trim = 1;
        M_trim = 3;
        N_trim = 28;

        % Gain source flag. false => load the stored gains from
        % trim_table_Poly_ConcatVer4p0.mat (default). true => redesign the
        % LQR gains from the Qlon/Rlon/Wlon (and lat) weights below via
        % RSLQR.set_gains at construction.
        update_gains = false;

        rho  = 0.0023769;   % slugs/ft^3
        grav = 32.17405;    % ft/sec^2
        ft2kts = 1/(1852.0/0.3048/3600);

        % Transition-end body x-velocity [ft/s]. Above this speed the offline
        % gain design (ctrl_lon/ctrl_lat) zeros the 8 lift-rotor B-columns so
        % the cruise controller uses aero surfaces + pusher only.
        % Value matches Trim_Ver4p0 (trans_end = 94.8 kts) used to build
        % trim_table_Poly_ConcatVer4p0.mat.
        trans_end = 94.8 / (1/(1852.0/0.3048/3600));   % 94.8 kts -> ft/s

        Qlon0 = [ 0.01 0.01 1000 0 0 0]'; % original 
        %Qlon0 = [ 0.02 0.02 100 0 0 0]'; % original 
        
        % Control acceleration cost
        Rlon0 = [1 1 1]'; % original
        
        % Control allocation weighting
        %       [ omp1:omp9         dele dflap    th] 
        Wlon0 = [ 1 1 1 1 1 1 1 1 1 1000 10000000 0.1]'; % Modified effector output order to match the simulation allocation
        %Wlon0 = [ 1 1 1 1 1 1 1 1 1 1000 10000000 1]'; % Modified effector output order to match the simulation allocation
        %Wlon0 = [ 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 1 1000 10000000 0.01]'; % Modified effector output order to match the simulation allocation
        
        Qlon = repmat(RSLQRConfig.Qlon0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);
        Rlon = repmat(RSLQRConfig.Rlon0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);
        Wlon = repmat(RSLQRConfig.Wlon0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);
        
        % State Cost 
        %         vi   pi   ri   v p r 
        Qlat0 = [ 0.01 1000 1000 0 0 0]'; 
        
        % Control acceleration cost
        Rlat0 = [1 1 1]';
        
        % Control allocation weighting
        %       [ omp1:omp8      dela delr phi]   (pusher omp9 is lon-only; not in lat)
        Wlat0 = [1 1 1 1 1 1 1 1 1000 1000 1]'; % Modified effector output order to match the simulation allocation
        %Wlat0 = [1 1 1 1 1 1 1 1 1000 1000 0.2]'; % Modified effector output order to match the simulation allocation
        
        Qlat = repmat(RSLQRConfig.Qlat0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);
        Rlat = repmat(RSLQRConfig.Rlat0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);
        Wlat = repmat(RSLQRConfig.Wlat0, [1,RSLQRConfig.N_trim,RSLQRConfig.M_trim,RSLQRConfig.L_trim]);

        Nx_lon = 4;
        Ni_lon = 3;
        Nu_lon = 11;
        Nr_lon = 3;
        Nv_lon = 1;
        
        Nx_lat = 4;
        Ni_lat = 3;
        Nu_lat = 10;
        Nr_lat = 2;
        Nv_lat = 1;

        eng_max = [1600 *ones([8, 1]); 2000] .* (2 * pi / 60);
        ele_max = deg2rad(30);
        ele_min = -deg2rad(30);
        flp_max = deg2rad(30);
        flp_min = -deg2rad(30);

        % NOTE: the servo-compensator discretization step (dt) is NOT stored
        % here. It is owned by SimConfig and injected into RSLQR at
        % construction (RSLQR(rslqrCfg, dt)) so sim and controller
        % share a single dt. This class holds gains/limits only.
    end

end