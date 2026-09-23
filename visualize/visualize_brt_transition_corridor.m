function visualize_brt_transition_corridor(idx_list)
% BRT transition corridor from the CURRENT reachability tube (single run, no
% comparison). Stacks the V=0 BRT boundary of every trim point along the
% trim-airspeed axis (u_trim). Four figures: LON u-w, LON q-theta, LAT v-phi,
% LAT p-r. Display units m/s and deg. The tube is the full-horizon (tau=0) slice
% of each guam_timestack .mat stack.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here);

    grids = struct('lon', axis_grid('lon'), 'lat', axis_grid('lat'));
    if nargin < 1 || isempty(idx_list), idx_list = 1:20; end

    brtRoot = fullfile(root, 'reachable_data', 'guam_timestack', 'BRT');
    outDir  = fullfile(root, 'reachable_data', 'corridor_figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    trimData = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp');
    XU0 = trimData.XU0_interp;
    wh_idx = 2;
    ft2m = 0.3048;

    cTube = [0.00 0.45 0.74];

    % {fig, axis, xdim, ydim, xlabel, ylabel, name, is2d}
    projs = {
        {1, 'lon', 1, 2, 'u [m/s]',   'w [m/s]',      'corridor_lon_uw',     true }
        {2, 'lon', 3, 4, 'q [deg/s]', '\theta [deg]', 'corridor_lon_qtheta', false}
        {3, 'lat', 1, 2, 'v [m/s]',   'p [deg/s]',    'corridor_lat_vp',     false}
        {4, 'lat', 3, 4, 'r [deg/s]', '\phi [deg]',   'corridor_lat_rphi',   false} };

    brt_line_width = 1.8;
    label_font_size = 15;
    legend_font_size = 15;
    
    n_proj = numel(projs);
    for p = 1:n_proj
        figure(projs{p}{1}); clf; set(gcf, 'Color', 'w');
        hold on;
        grid on;
        box on;
    end

    % One representative handle per figure so the legend shows a single entry per
    % category (the loop draws one boundary + one target box per UH, which would
    % otherwise stack into ~20 duplicate legend rows).
    h_tube   = gobjects(1, n_proj);
    h_target = gobjects(1, n_proj);

    for idx = idx_list
        trim = XU0(:, idx, wh_idx);
        zlev = trim(1)*ft2m;

        V.lon = load_tube(fullfile(brtRoot, 'BRT_LON_MAT', sprintf('GUAM_LON_BRT_UH%d_WH%d_stack.mat', idx, wh_idx)));
        V.lat = load_tube(fullfile(brtRoot, 'BRT_LAT_MAT', sprintf('GUAM_LAT_BRT_UH%d_WH%d_stack.mat', idx, wh_idx)));

        for p = 1:n_proj
            fnum = projs{p}{1}; axisName = projs{p}{2}; xdim = projs{p}{3}; ydim = projs{p}{4};
            is2d = projs{p}{8};
            axCfg = grids.(axisName);

            [xv, yv, Z] = brt_project(V.(axisName), axCfg.gv, xdim, ydim);
            tx = trim(axCfg.trim_rows(xdim));  sx = axCfg.scale(xdim);
            ty = trim(axCfg.trim_rows(ydim));  sy = axCfg.scale(ydim);

            figure(fnum);
            h_boundary = stack_segs(brt_zero_contour(xv, yv, Z), tx, ty, sx, sy, zlev, cTube, is2d, brt_line_width);
            if ~isempty(h_boundary), h_tube(p) = h_boundary; end

            xlb = (axCfg.tlb(xdim) + tx)*sx;  xub = (axCfg.tub(xdim) + tx)*sx;
            ylb = (axCfg.tlb(ydim) + ty)*sy;  yub = (axCfg.tub(ydim) + ty)*sy;
            target_set_line_width = 1.5;
            if is2d
                h_target(p) = plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], ...
                                   'g--', 'LineWidth', target_set_line_width);
            else
                h_target(p) = plot3([xlb xub xub xlb xlb], zlev*ones(1,5), [ylb ylb yub yub ylb], ...
                                    'g--', 'LineWidth', target_set_line_width);
            end
        end
    end

    for p = 1:n_proj
        figure(projs{p}{1});
        if projs{p}{8}
            xlabel(projs{p}{5}, 'FontSize', label_font_size);
            ylabel(projs{p}{6}, 'FontSize', label_font_size);
            pbaspect([3, 1, 1]);
        else
            xlabel(projs{p}{5}, 'FontSize', label_font_size);
            ylabel('u_{trim} [m/s]', 'FontSize', label_font_size);
            zlabel(projs{p}{6}, 'FontSize', label_font_size);
            view([60, 20]);
            pbaspect([1, 6, 1]);
        end
        legend_handles = [h_tube(p), h_target(p)];
        legend_labels  = {'BRT boundary', 'target set'};
        shown = isgraphics(legend_handles);
        % legend(legend_handles(shown), legend_labels(shown), 'Location', 'best', 'FontSize', legend_font_size);
        exportgraphics(gcf, fullfile(outDir, [projs{p}{7} '.png']), 'Resolution', 150);
        fprintf('saved %s\n', fullfile(outDir, [projs{p}{7} '.png']));
    end
end

%%
function h = stack_segs(segs, tx, ty, sx, sy, zlev, color, is2d, line_width)
    h = gobjects(0);
    if isempty(segs), return; end
    s = segs{end};
    if is2d
        h = plot((s(1,:) + tx)*sx, (s(2,:) + ty)*sy, '-', 'Color', color, 'LineWidth', line_width);
    else
        h = plot3((s(1,:) + tx)*sx, zlev*ones(1, size(s,2)), (s(2,:) + ty)*sy, '-', ...
                'Color', color, 'LineWidth', line_width);
    end
end

function V = load_tube(fname)
    % Full-horizon BRT tube = first (tau=0) slice of the time-stack.
    s = load(fname, 'Vslices');
    V = s.Vslices{1};
end

function g = axis_grid(ax)
    % Native grid vectors + plot metadata, all derived from FilterConfig.axisSpec.
    spec = FilterConfig.axisSpec(ax);
    ft2m = 0.3048;  r2d = 180 / pi;
    g.gv        = arrayfun(@(a, b, n) linspace(a, b, n), ...
                           spec.grid_min, spec.grid_max, spec.grid_num, 'UniformOutput', false);
    g.trim_rows = spec.trim_rows(:);
    g.tlb       = -spec.target_ub(:);
    g.tub       =  spec.target_ub(:);
    if strcmpi(ax, 'lon')
        g.scale = [ft2m; ft2m; r2d; r2d];
    else
        g.scale = [ft2m; r2d; r2d; r2d];
    end
end
