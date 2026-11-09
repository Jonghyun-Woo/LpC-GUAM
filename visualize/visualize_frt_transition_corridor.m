function visualize_frt_transition_corridor(idx_list)
% FRT transition corridor: stacks the V=0 FRT boundary of every trim point along
% the trim-airspeed axis (u_trim). Four figures: LON u-w, LON q-theta, LAT v-phi,
% LAT p-r. Display units m/s and deg. Source: guam_output (WH2, converged FRT).
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here);

    grids = struct('lon', axis_grid('lon'), 'lat', axis_grid('lat'));
    if nargin < 1 || isempty(idx_list), idx_list = 1:20; end

    frtRoot = fullfile(root, 'reachability_data', 'guam_output');
    outDir  = fullfile(root, 'reachability_data', 'corridor_figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    trimData = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp');
    XU0 = trimData.XU0_interp;
    wh_idx = 2;
    ft2m = 0.3048;

    colors = lines(2);
    cTube = colors(2, :);

    projs = {
        {1, 'lon', 1, 2, 'u [m/s]',   'w [m/s]',      'corridor_frt_lon_uw',     true }
        {2, 'lon', 3, 4, 'q [deg/s]', '\theta [deg]', 'corridor_frt_lon_qtheta', false}
        {3, 'lat', 1, 2, 'v [m/s]',   'p [deg/s]',    'corridor_frt_lat_vp',     false}
        {4, 'lat', 3, 4, 'r [deg/s]', '\phi [deg]',   'corridor_frt_lat_rphi',   false} };

    frt_line_width = 1.8;
    target_set_line_width = 1.5;
    label_font_size = 15;

    n_proj = numel(projs);
    for p = 1:n_proj
        figure(projs{p}{1}); clf; set(gcf, 'Color', 'w');
        hold on;
        grid on;
        box on;
    end

    for idx = idx_list
        trim = XU0(:, idx, wh_idx);
        zlev = trim(1)*ft2m;

        [V.lon, gv.lon] = load_tube(fullfile(frtRoot, 'LON_FRT', sprintf('GUAM_LON_FRT_UH%d_WH%d.mat', idx, wh_idx)));
        [V.lat, gv.lat] = load_tube(fullfile(frtRoot, 'LAT_FRT', sprintf('GUAM_LAT_FRT_UH%d_WH%d.mat', idx, wh_idx)));

        for p = 1:n_proj
            fnum = projs{p}{1}; axisName = projs{p}{2}; xdim = projs{p}{3}; ydim = projs{p}{4};
            is2d = projs{p}{8};
            axCfg = grids.(axisName);

            [xv, yv, Z] = brt_project(V.(axisName), gv.(axisName), xdim, ydim);
            tx = trim(axCfg.trim_rows(xdim));  sx = axCfg.scale(xdim);
            ty = trim(axCfg.trim_rows(ydim));  sy = axCfg.scale(ydim);

            figure(fnum);
            stack_segs(brt_zero_contour(xv, yv, Z), tx, ty, sx, sy, zlev, cTube, is2d, frt_line_width);

            xlb = (axCfg.tlb(xdim) + tx)*sx;  xub = (axCfg.tub(xdim) + tx)*sx;
            ylb = (axCfg.tlb(ydim) + ty)*sy;  yub = (axCfg.tub(ydim) + ty)*sy;
            if is2d
                plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], 'g--', 'LineWidth', target_set_line_width);
            else
                plot3([xlb xub xub xlb xlb], zlev*ones(1,5), [ylb ylb yub yub ylb], 'g--', 'LineWidth', target_set_line_width);
            end
        end
    end

    for p = 1:n_proj
        figure(projs{p}{1});
        xlabel(projs{p}{5}, 'FontSize', label_font_size);
        if projs{p}{8}
            ylabel(projs{p}{6}, 'FontSize', label_font_size);
            pbaspect([3, 1, 1]);
        else
            ylabel('u_{trim} [m/s]', 'FontSize', label_font_size);
            zlabel(projs{p}{6}, 'FontSize', label_font_size);
            view([60, 20]);
            pbaspect([1, 6, 1]);
        end
        exportgraphics(gcf, fullfile(outDir, [projs{p}{7} '.png']), 'Resolution', 150);
        fprintf('saved %s\n', fullfile(outDir, [projs{p}{7} '.png']));
    end
end

%%
function stack_segs(segs, tx, ty, sx, sy, zlev, color, is2d, line_width)
    if isempty(segs), return; end
    s = segs{end};
    if is2d
        plot((s(1,:) + tx)*sx, (s(2,:) + ty)*sy, '-', 'Color', color, 'LineWidth', line_width);
    else
        plot3((s(1,:) + tx)*sx, zlev*ones(1, size(s,2)), (s(2,:) + ty)*sy, '-', ...
              'Color', color, 'LineWidth', line_width);
    end
end

function [V, grid_vectors] = load_tube(fname)
    s = load(fname, 'values', 'grid_axes');
    V = double(s.values);
    grid_vectors = s.grid_axes;
end

function g = axis_grid(ax)
    spec = FilterConfig.axisSpec(ax);
    ft2m = 0.3048;  r2d = 180 / pi;
    g.trim_rows = spec.trim_rows(:);
    g.tlb       = -spec.target_ub(:);
    g.tub       =  spec.target_ub(:);
    if strcmpi(ax, 'lon')
        g.scale = [ft2m; ft2m; r2d; r2d];
    else
        g.scale = [ft2m; r2d; r2d; r2d];
    end
end
