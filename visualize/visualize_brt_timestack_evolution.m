function visualize_brt_timestack_evolution()
% Draw how the longitudinal/lateral BRT shrinks over the finite horizon, using
% the time-stacked value functions in reachable_data/guam_timestack/BRT. One
% square figure per requested time slice; axis limits are held fixed across the
% slices (sized to the largest, i.e. longest-horizon, tube) so the temporal
% shrinkage is visible directly. Coordinates are deviation-from-trim in display
% units (m/s, deg), centred at trim. The target set is drawn as a green dashed
% box in every figure.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here);

    dataRoot = fullfile(root, 'reachable_data', 'guam_timestack', 'BRT');
    outDir   = fullfile(root, 'reachable_data', 'brt_timestack_evolution');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    uh = 5;  wh = 2;
    req_times = [0.6, 0.7, 0.8, 0.9, 1.0];    % s; snapped to nearest tau slice

    ft2m = 0.3048;  r2d = 180/pi;
    col_brt    = [0.00 0.45 0.74];
    col_target = [0.10 0.60 0.10];

    axes_cfg = struct( ...
        'axis',    {'LON',            'LAT'          }, ...
        'sub',     {'BRT_LON_MAT',    'BRT_LAT_MAT'  }, ...
        'xdim',    {3,                2              }, ...
        'ydim',    {4,                4              }, ...
        'scale',   {[ft2m ft2m],      [ft2m r2d]     }, ...
        'tlb',     {[-0.10 -0.10],    [-0.10 -0.10]  }, ...
        'tub',     {[ 0.10  0.10],    [ 0.10  0.10]  }, ...
        'xlab',    {' ',  ' '  }, ...
        'ylab',    {' ',  ' '}, ...
        'tag',     {' ',             ' '           });

    for ac = axes_cfg
        fname = fullfile(dataRoot, ac.sub, ...
                         sprintf('GUAM_%s_BRT_UH%d_WH%d_stack.mat', ac.axis, uh, wh));
        s  = load(fname, 'Vslices', 'taus', 'grid_axes');
        gv = s.grid_axes;

        % snap requested times to the nearest available tau slice
        tidx = arrayfun(@(t) nearest_idx(s.taus, t), req_times);

        % fixed square limits from the largest tube across all shown slices
        [Lx, Ly] = tube_extent(s.Vslices(tidx), gv, ac);

        for j = 1:numel(tidx)
            V   = s.Vslices{tidx(j)};
            tau = s.taus(tidx(j));

            [xv, yv, Z] = brt_project(V, gv, ac.xdim, ac.ydim);
            segs = brt_zero_contour(xv, yv, Z);

            figure; clf; set(gcf, 'Color', 'w'); hold on; grid on; box on;

            draw_target_box(ac);
            for k = 1:numel(segs)
                seg = segs{k};
                plot(seg(1,:)*ac.scale(1), seg(2,:)*ac.scale(2), '-', ...
                     'Color', col_brt, 'LineWidth', 2.0);
            end

            xlim([-Lx Lx]);  ylim([-Ly Ly]);  pbaspect([1 1 1]);
            xlabel(ac.xlab, 'FontSize', 14);
            ylabel(ac.ylab, 'FontSize', 14);
            title(sprintf('%s BRT  UH%d  \\tau = %.2f s', ac.axis, uh, tau), ...
                  'FontSize', 14);

            outName = fullfile(outDir, sprintf('brt_evo_%s_%s_UH%d_t%03.0fms.png', ...
                               ac.axis, ac.tag, uh, tau*1000));
            exportgraphics(gcf, outName, 'Resolution', 150);
            fprintf('saved %s\n', outName);
            close(gcf);
        end
    end

    function draw_target_box(ac)
        xlb = ac.tlb(1)*ac.scale(1);  xub = ac.tub(1)*ac.scale(1);
        ylb = ac.tlb(2)*ac.scale(2);  yub = ac.tub(2)*ac.scale(2);
        plot([xlb xub xub xlb xlb], [ylb ylb yub yub ylb], '--', ...
             'Color', col_target, 'LineWidth', 1.6);
    end
end

function i0 = nearest_idx(v, t)
    [~, i0] = min(abs(v - t));
end

function [Lx, Ly] = tube_extent(Vcells, gv, ac)
% Half-extent covering every shown tube's zero contour in display units, with
% margin. Clamped to at least the target box so it is always visible.
    xmax = ac.tub(1);  ymax = ac.tub(2);
    for i = 1:numel(Vcells)
        [xv, yv, Z] = brt_project(Vcells{i}, gv, ac.xdim, ac.ydim);
        segs = brt_zero_contour(xv, yv, Z);
        for k = 1:numel(segs)
            seg = segs{k};
            xmax = max(xmax, max(abs(seg(1,:))));
            ymax = max(ymax, max(abs(seg(2,:))));
        end
    end
    Lx = 1.1 * xmax * ac.scale(1);
    Ly = 1.1 * ymax * ac.scale(2);
end
