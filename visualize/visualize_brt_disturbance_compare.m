here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here);

S = brt_grid_config();
idx_list = S.UH_idx;

distRoot   = fullfile(root, 'reachable_data', 'GUAM_BRT_run_0903');
nodistRoot = fullfile(root, 'reachable_data', 'GUAM_BRT_run_0821_wo_disturbance');
outDir     = fullfile(root, 'reachable_data', 'compare_figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

tt  = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp');
XU0 = tt.XU0_interp;
wh  = S.WH_idx;
keepOpen = numel(idx_list) <= 6;

projs = {{'lon',1,2}, {'lon',3,4}, {'lat',1,4}, {'lat',2,3}};
cND = [0.00 0.45 0.74];
cD  = [0.85 0.10 0.10];

for idx = idx_list
    trim = XU0(:, idx, wh);

    Vlon_d  = readNPY(fullfile(distRoot,   'LON_NPY', sprintf('GUAM_LON_BRT_UH%d_WH%d.npy', idx, wh)));
    Vlon_nd = readNPY(fullfile(nodistRoot, 'LON_NPY', sprintf('GUAM_LON_BRT_UH%d_WH%d.npy', idx, wh)));
    Vlat_d  = readNPY(fullfile(distRoot,   'LAT_NPY', sprintf('GUAM_LAT_BRT_UH%d_WH%d.npy', idx, wh)));
    Vlat_nd = readNPY(fullfile(nodistRoot, 'LAT_NPY', sprintf('GUAM_LAT_BRT_UH%d_WH%d.npy', idx, wh)));
    
    f = figure('Name', sprintf('BRT compare UH%d WH%d', idx, wh), ...
                'Color', 'w', 'Position', [80 80 1000 820]);

    for p = 1:4
        axis_name = projs{p}{1};
        xdim = projs{p}{2};  ydim = projs{p}{3};
        g   = S.(axis_name);
        if strcmp(axis_name, 'lon'), Vnd = Vlon_nd; Vd = Vlon_d;
        else,                        Vnd = Vlat_nd; Vd = Vlat_d;  end

        subplot(2, 2, p);
        plot_projection(g, Vnd, Vd, xdim, ydim, trim, cND, cD);
    end

    sgtitle(sprintf('GUAM BRT  UH%d (u_{trim}=%.1f m/s, w_{trim}=%.1f m/s)', ...
                    idx, trim(1)*0.3048, trim(3)*0.3048), 'FontWeight', 'bold');

    pngPath = fullfile(outDir, sprintf('BRT_compare_UH%d_WH%d.png', idx, wh));
    exportgraphics(f, pngPath, 'Resolution', 150);
    % if ~keepOpen, close(f); end
    fprintf('saved %s\n', pngPath);
end

function plot_projection(g, Vnd, Vd, xdim, ydim, trim, cND, cD)
    [xv, yv, Znd] = brt_project(Vnd, g.gv, xdim, ydim);
    [~,  ~,  Zd ] = brt_project(Vd,  g.gv, xdim, ydim);

    tx = trim(g.trim_rows(xdim));  sx = g.scale(xdim);
    ty = trim(g.trim_rows(ydim));  sy = g.scale(ydim);

    hold on; grid on; box on;
    plot_segs(brt_zero_contour(xv, yv, Znd), tx, ty, sx, sy, cND);
    plot_segs(brt_zero_contour(xv, yv, Zd ), tx, ty, sx, sy, cD);

    xlb = (g.tlb(xdim) + tx)*sx;  xub = (g.tub(xdim) + tx)*sx;
    ylb = (g.tlb(ydim) + ty)*sy;  yub = (g.tub(ydim) + ty)*sy;
    plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], 'g--', 'LineWidth', 1.4);
    plot(tx*sx, ty*sy, 'k+', 'MarkerSize', 9, 'LineWidth', 1.2);

    xlabel(g.labels{xdim});  ylabel(g.labels{ydim});
end

function plot_segs(segs, tx, ty, sx, sy, color)
    % for k = 1:numel(segs)
    %     s = segs{k};
    %     plot((s(1,:) + tx)*sx, (s(2,:) + ty)*sy, '-', 'Color', color, 'LineWidth', 1.6);
    % end
    s = segs{end};
    plot((s(1,:) + tx)*sx, (s(2,:) + ty)*sy, '-', 'Color', color, 'LineWidth', 1.6);
end
