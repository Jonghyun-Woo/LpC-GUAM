function visualize_brt_transition_corridor_compare(idx_list)
% BRT transition corridors, with vs without disturbance. Companion to
% helperOC/visualize_transition_corridor.m: stacks the V=0 BRT boundary of
% every trim point along the trim-airspeed axis (u_trim), overlaying the
% with-disturbance run (red) and no-disturbance run (blue). Four figures:
% LON u-w, LON q-theta, LAT v-phi, LAT p-r. Display units m/s and deg.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here);

    S = brt_grid_config();
    if nargin < 1 || isempty(idx_list), idx_list = S.UH_idx; end

    distRoot   = fullfile(root, 'reachable_data', 'GUAM_BRT_run_0903');
    nodistRoot = fullfile(root, 'reachable_data', 'GUAM_BRT_run_0821_wo_disturbance');
    outDir     = fullfile(root, 'reachable_data', 'compare_figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    tt  = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp');
    XU0 = tt.XU0_interp;
    wh  = S.WH_idx;
    ft2m = 0.3048;

    cND = [0.00 0.45 0.74];
    cD  = [0.85 0.10 0.10];

    % {fig, axis, xdim, ydim, xlabel, ylabel, name, is2d}
    projs = {
        {1, 'lon', 1, 2, 'u [m/s]',   'w [m/s]',      'corridor_lon_uw',     true }
        {2, 'lon', 3, 4, 'q [deg/s]', '\theta [deg]', 'corridor_lon_qtheta', false}
        {3, 'lat', 1, 4, 'v [m/s]',   '\phi [deg]',   'corridor_lat_vphi',   false}
        {4, 'lat', 2, 3, 'p [deg/s]', 'r [deg/s]',    'corridor_lat_pr',     false} };

    for p = 1:numel(projs)
        figure(projs{p}{1}); clf; set(gcf, 'Color', 'w'); hold on; grid on; box on;
    end

    for idx = idx_list
        trim = XU0(:, idx, wh);
        zlev = trim(1)*ft2m;

        V.lon_d  = readNPY(fullfile(distRoot,   'LON_NPY', sprintf('GUAM_LON_BRT_UH%d_WH%d.npy', idx, wh)));
        V.lon_nd = readNPY(fullfile(nodistRoot, 'LON_NPY', sprintf('GUAM_LON_BRT_UH%d_WH%d.npy', idx, wh)));
        V.lat_d  = readNPY(fullfile(distRoot,   'LAT_NPY', sprintf('GUAM_LAT_BRT_UH%d_WH%d.npy', idx, wh)));
        V.lat_nd = readNPY(fullfile(nodistRoot, 'LAT_NPY', sprintf('GUAM_LAT_BRT_UH%d_WH%d.npy', idx, wh)));

        for p = 1:numel(projs)
            fnum = projs{p}{1}; ax = projs{p}{2}; xdim = projs{p}{3}; ydim = projs{p}{4};
            is2d = projs{p}{8};
            g   = S.(ax);
            Vnd = V.(sprintf('%s_nd', ax));
            Vd  = V.(sprintf('%s_d',  ax));

            [xv, yv, Znd] = brt_project(Vnd, g.gv, xdim, ydim);
            [~,  ~,  Zd ] = brt_project(Vd,  g.gv, xdim, ydim);
            tx = trim(g.trim_rows(xdim));  sx = g.scale(xdim);
            ty = trim(g.trim_rows(ydim));  sy = g.scale(ydim);

            figure(fnum);
            stack_segs(brt_zero_contour(xv, yv, Znd), tx, ty, sx, sy, zlev, cND, is2d);
            stack_segs(brt_zero_contour(xv, yv, Zd ), tx, ty, sx, sy, zlev, cD,  is2d);

            xlb = (g.tlb(xdim) + tx)*sx;  xub = (g.tub(xdim) + tx)*sx;
            ylb = (g.tlb(ydim) + ty)*sy;  yub = (g.tub(ydim) + ty)*sy;
            if is2d
                plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], 'g--', 'LineWidth', 1.0);
            else
                plot3([xlb xub xub xlb xlb], zlev*ones(1,5), [ylb ylb yub yub ylb], 'g--', 'LineWidth', 1.0);
            end
        end
    end

    for p = 1:numel(projs)
        figure(projs{p}{1});
        title(sprintf('BRT Transition Corridor: %s', projs{p}{7}), 'Interpreter', 'none');
        if projs{p}{8}
            xlabel(projs{p}{5}); ylabel(projs{p}{6});
            pbaspect([3, 1, 1]);
        else
            xlabel(projs{p}{5}); ylabel('u_{trim} [m/s]'); zlabel(projs{p}{6});
            view([60, 20]);
            pbaspect([1, 4, 1]);
        end
        exportgraphics(gcf, fullfile(outDir, [projs{p}{7} '.png']), 'Resolution', 150);
        fprintf('saved %s\n', fullfile(outDir, [projs{p}{7} '.png']));
    end
end

function stack_segs(segs, tx, ty, sx, sy, zlev, color, is2d)
    % for k = 1:numel(segs)
    s = segs{end};
    if is2d
        plot((s(1,:) + tx)*sx, (s(2,:) + ty)*sy, '-', 'Color', color, 'LineWidth', 1.3);
    else
        plot3((s(1,:) + tx)*sx, zlev*ones(1, size(s,2)), (s(2,:) + ty)*sy, '-', ...
                'Color', color, 'LineWidth', 1.3);
    end
    % end
end
