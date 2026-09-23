function visualize_reachable_tube_pair()
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here);

    dataRoot = fullfile(root, 'reachable_data', 'guam_output');
    outDir   = fullfile(root, 'reachable_data', 'tube_pair_figures');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    grids    = struct('lon', axis_grid('lon'), 'lat', axis_grid('lat'));
    trimData = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp');
    XU0      = trimData.XU0_interp;
    wh_idx   = 2;
    uh_all   = 1:20;

    col_brt    = [0.00 0.45 0.74];
    col_frt    = [0.85 0.33 0.10];
    col_target = [0.10 0.60 0.10];
    tube_line_width   = 1.8;
    target_line_width = 1.5;
    label_font_size   = 14;
    legend_font_size  = 12;

    projs = struct( ...
        'axis', {'lon',       'lon',          'lat',       'lat'         }, ...
        'xdim', {1,           3,              1,           3             }, ...
        'ydim', {2,           4,              2,           4             }, ...
        'xlab', {'u [m/s]',   'q [deg/s]',    'v [m/s]',   'r [deg/s]'   }, ...
        'ylab', {'w [m/s]',   '\theta [deg]', 'p [deg/s]', '\phi [deg]'  }, ...
        'tag',  {'lon_uw',    'lon_qtheta',   'lat_vp',    'lat_rphi'    });

    for k = 1:numel(uh_all) - 1
        uh_first  = uh_all(k);
        uh_second = uh_all(k + 1);

        val_first.lon.brt  = load_value(dataRoot, 'lon', 'BRT', uh_first,  wh_idx);
        val_first.lon.frt  = load_value(dataRoot, 'lon', 'FRT', uh_first,  wh_idx);
        val_first.lat.brt  = load_value(dataRoot, 'lat', 'BRT', uh_first,  wh_idx);
        val_first.lat.frt  = load_value(dataRoot, 'lat', 'FRT', uh_first,  wh_idx);
        val_second.lon.brt = load_value(dataRoot, 'lon', 'BRT', uh_second, wh_idx);
        val_second.lat.brt = load_value(dataRoot, 'lat', 'BRT', uh_second, wh_idx);

        trim_first  = XU0(:, uh_first,  wh_idx);
        trim_second = XU0(:, uh_second, wh_idx);

        for pj = 1:numel(projs)
            axisName = projs(pj).axis;  xdim = projs(pj).xdim;  ydim = projs(pj).ydim;
            axCfg    = grids.(axisName);

            figure; clf; set(gcf, 'Color', 'w'); hold on; grid on; box on;

            h_brt_first  = draw_tube_boundary(val_first.(axisName).brt, axCfg, xdim, ydim, ...
                                              trim_first,  '-',  col_brt,    tube_line_width);
            h_frt_first  = draw_tube_boundary(val_first.(axisName).frt, axCfg, xdim, ydim, ...
                                              trim_first,  '-',  col_frt,    tube_line_width);
            h_tgt_first  = draw_target_box(axCfg, xdim, ydim, ...
                                           trim_first,  '-',  col_target, target_line_width);
            h_brt_second = draw_tube_boundary(val_second.(axisName).brt, axCfg, xdim, ydim, ...
                                              trim_second, '--', col_brt,    tube_line_width);
            h_tgt_second = draw_target_box(axCfg, xdim, ydim, ...
                                           trim_second, '--', col_target, target_line_width);

            xlabel(projs(pj).xlab, 'FontSize', label_font_size);
            ylabel(projs(pj).ylab, 'FontSize', label_font_size);

            handles = [h_brt_first, h_frt_first, h_tgt_first, h_brt_second, h_tgt_second];
            labels  = {sprintf('BRT UH%d', uh_first),    sprintf('FRT UH%d', uh_first), ...
                       sprintf('target UH%d', uh_first), sprintf('BRT UH%d', uh_second), ...
                       sprintf('target UH%d', uh_second)};
            shown = isgraphics(handles);
            legend(handles(shown), labels(shown), 'Location', 'northeast', 'FontSize', legend_font_size);

            fname = fullfile(outDir, sprintf('tubepair_UH%d_%d_%s.png', ...
                             uh_first, uh_second, projs(pj).tag));
            exportgraphics(gcf, fname, 'Resolution', 150);
            fprintf('saved %s\n', fname);
            close(gcf);
        end
    end
end

function V = load_value(dataRoot, axisName, mode, uh, wh)
    AX = upper(axisName);  MODE = upper(mode);
    fname = fullfile(dataRoot, sprintf('%s_%s', AX, MODE), ...
                     sprintf('GUAM_%s_%s_UH%d_WH%d.mat', AX, MODE, uh, wh));
    s = load(fname, 'values');
    V = s.values;
end

function h = draw_tube_boundary(V, axCfg, xdim, ydim, trim, style, color, line_width)
    h = gobjects(1);
    [xv, yv, Z] = brt_project(V, axCfg.gv, xdim, ydim);
    segs = brt_zero_contour(xv, yv, Z);
    if isempty(segs), return; end

    tx = trim(axCfg.trim_rows(xdim));  sx = axCfg.scale(xdim);
    ty = trim(axCfg.trim_rows(ydim));  sy = axCfg.scale(ydim);

    X = [];  Y = [];
    for k = 1:numel(segs)
        s = segs{k};
        X = [X, (s(1,:) + tx)*sx, NaN];  %#ok<AGROW>
        Y = [Y, (s(2,:) + ty)*sy, NaN];  %#ok<AGROW>
    end
    h = plot(X, Y, style, 'Color', color, 'LineWidth', line_width);
end

function h = draw_target_box(axCfg, xdim, ydim, trim, style, color, line_width)
    tx = trim(axCfg.trim_rows(xdim));  sx = axCfg.scale(xdim);
    ty = trim(axCfg.trim_rows(ydim));  sy = axCfg.scale(ydim);
    xlb = (axCfg.tlb(xdim) + tx)*sx;  xub = (axCfg.tub(xdim) + tx)*sx;
    ylb = (axCfg.tlb(ydim) + ty)*sy;  yub = (axCfg.tub(ydim) + ty)*sy;
    h = plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], style, ...
             'Color', color, 'LineWidth', line_width);
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
