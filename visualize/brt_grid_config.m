function S = brt_grid_config()
% Grid / target-set definition for the GUAM BRT tables (see
% reachable_data/GUAM_BRT_run_*/guam_analysis_config.yml). Grid vectors are in
% deviation-from-trim native units (ft, rad); scale converts to display units
% (m, deg). trim_rows map each state to its trim value in XU0_interp(:,idx,wh).
    ft2m = 0.3048;  r2d = 180/pi;

    S.lon.gmin = [-16.4042; -21.3255; -0.6109; -0.4363];
    S.lon.gmax = [ 16.4042;  21.3255;  0.6109;  0.4363];
    S.lon.gnum = [      33;       51;       49;       33];
    S.lon.tlb  = [    -3.0;     -3.0;    -0.10;    -0.10];
    S.lon.tub  = [     3.0;      3.0;     0.10;     0.10];
    S.lon.trim_rows = [1; 3; 5; 11];          % u, w, q, theta
    S.lon.scale  = [ft2m; ft2m; r2d; r2d];
    S.lon.labels = {'u [m/s]','w [m/s]','q [deg/s]','\theta [deg]'};

    S.lat.gmin = [-16.4042; -1.2217; -0.6109; -0.6981];
    S.lat.gmax = [ 16.4042;  1.2217;  0.6109;  0.6981];
    S.lat.gnum = [      33;      85;      49;      39];
    S.lat.tlb  = [    -3.0;   -0.10;   -0.10;   -0.10];
    S.lat.tub  = [     3.0;    0.10;    0.10;    0.10];
    S.lat.trim_rows = [2; 4; 6; 10];          % v, p, r, phi
    S.lat.scale  = [ft2m; r2d; r2d; r2d];
    S.lat.labels = {'v [m/s]','p [deg/s]','r [deg/s]','\phi [deg]'};

    for ax = ["lon","lat"]
        g = S.(ax);
        g.gv = arrayfun(@(a,b,n) linspace(a,b,n), g.gmin, g.gmax, g.gnum, ...
                        'UniformOutput', false);
        S.(ax) = g;
    end

    S.UH_idx = 1:20;
    S.WH_idx = 3;
end
