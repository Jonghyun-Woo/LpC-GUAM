% Visualizes the Monte-Carlo reach test written by tests/test_mc_verify.m.
% Per case: trajectories projected onto the two sampled states, one tile per UH,
% truncated at first target entry, framed on the BRT tube boundary.
% Figures go to reachability_data/mc_verify_timestack/figures/.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
dataDir = fullfile(root, 'reachability_data', 'mc_verify_timestack');
figDir  = fullfile(dataDir, 'figures');
CASES   = {'uw', 'qtheta', 'vp', 'rphi'};
MARGIN  = 0.08;
% -------------------------------------------------------------------------

if ~exist(figDir, 'dir'), mkdir(figDir); end
ft2m = 0.3048;  r2d = 180 / pi;

view = struct( ...
    'axis',  {'lon', 'lat'}, ...
    'prow',  {[4; 6; 11; 8], [5; 10; 12; 7]}, ...
    'trow',  {[1; 3; 5; 11], [2; 4; 6; 10]}, ...
    'scale', {[ft2m; ft2m; r2d; r2d], [ft2m; r2d; r2d; r2d]}, ...
    'label', {{'\Deltau [m/s]', '\Deltaw [m/s]', '\Deltaq [deg/s]', '\Delta\theta [deg]'}, ...
              {'\Deltav [m/s]', '\Deltap [deg/s]', '\Deltar [deg/s]', '\Delta\phi [deg]'}});

cReach = [0.00 0.45 0.74];
cFail  = [0.85 0.10 0.10];
rate   = struct('case', {}, 'uh', {}, 'pct', {});

for ic = 1:numel(CASES)
    files = dir(fullfile(dataDir, sprintf('mc_%s_UH*.mat', CASES{ic})));
    if isempty(files)
        warning('viz_mc_verify:missing', 'no files for case %s', CASES{ic});
        continue;
    end
    uh_all = arrayfun(@(f) sscanf(f.name, ['mc_' CASES{ic} '_UH%d.mat']), files);
    [~, ord] = sort(uh_all);
    files = files(ord);

    nt  = numel(files);
    fig = figure('Name', sprintf('MC reach: %s', CASES{ic}), ...
                 'Position', [40 40 1400 400 * ceil(nt / 3)]);
    tl = tiledlayout(fig, ceil(nt / 3), min(3, nt), 'Padding', 'compact', ...
                     'TileSpacing', 'compact');

    for k = 1:nt
        S  = load(fullfile(dataDir, files(k).name));
        vw = view(1 + strcmp(S.axis, 'lat'));
        d1 = S.free_dims(1);  d2 = S.free_dims(2);
        s1 = vw.scale(d1);    s2 = vw.scale(d2);
        n  = numel(S.reach);

        dev1 = (S.X(:, :, vw.prow(d1)) - S.trim_state(vw.trow(d1))) * s1;
        dev2 = (S.X(:, :, vw.prow(d2)) - S.trim_state(vw.trow(d2))) * s2;
        stop = size(S.X, 2) * ones(n, 1);
        stop(S.reach) = round(S.t_reach(S.reach) / S.dt) + 1;

        nexttile(tl);  hold on; grid on; box on;
        for s = 1:n
            if S.reach(s), col = cReach; else, col = cFail; end
            plot(dev1(s, 1:stop(s)), dev2(s, 1:stop(s)), '-', ...
                 'Color', [col 0.20], 'LineWidth', 0.5);
        end
        scatter(S.x0_dev(~S.reach, d1) * s1, S.x0_dev(~S.reach, d2) * s2, 10, cFail,  'filled');
        scatter(S.x0_dev( S.reach, d1) * s1, S.x0_dev( S.reach, d2) * s2, 10, cReach, 'filled');

        e1 = arrayfun(@(s) dev1(s, stop(s)), (1:n)');
        e2 = arrayfun(@(s) dev2(s, stop(s)), (1:n)');
        plot(e1, e2, 'k.', 'MarkerSize', 7);

        b1 = S.target_ub(d1) * s1;  b2 = S.target_ub(d2) * s2;
        plot([-b1 b1 b1 -b1 -b1], [-b2 -b2 b2 b2 -b2], 'g-', 'LineWidth', 1.8);

        pts = [S.brt_contour{:}];
        if isempty(pts), pts = [S.x0_dev(:, d1)'; S.x0_dev(:, d2)']; end
        for q = 1:numel(S.brt_contour)
            seg = S.brt_contour{q};
            plot(seg(1, :) * s1, seg(2, :) * s2, 'k-', 'LineWidth', 1.0);
        end
        xlim(pad([min(pts(1, :)) max(pts(1, :))] * s1, MARGIN));
        ylim(pad([min(pts(2, :)) max(pts(2, :))] * s2, MARGIN));

        xlabel(vw.label{d1});  ylabel(vw.label{d2});
        title(sprintf('UH%d (u_{trim} = %.1f m/s)', ...
                      S.uh_idx, S.uh_vel * ft2m));

        rate(end + 1) = struct('case', CASES{ic}, 'uh', S.uh_idx, ...
                               'pct', 100 * mean(S.reach)); %#ok<AGROW>
    end

    title(tl, sprintf(['BRT reach test, case %s   ' ...
                       '(dashed = BRT boundary, green = target set, black = final state)'], CASES{ic}));
    exportgraphics(fig, fullfile(figDir, sprintf('mc_reach_%s.png', CASES{ic})), 'Resolution', 150);
end

fig = figure('Name', 'MC reach summary', 'Color', 'w', 'Position', [60 60 720 460]);
hold on; grid on; box on;
mk = {'-o', '-s', '-^', '-d'};
for ic = 1:numel(CASES)
    sel = rate(strcmp({rate.case}, CASES{ic}));
    if isempty(sel), continue; end
    plot([sel.uh], [sel.pct], mk{ic}, 'LineWidth', 1.6, 'MarkerFaceColor', 'w');
end
ylim([0 105]);  xlabel('UH index');  ylabel('reach rate [%]');
legend(CASES, 'Location', 'southwest');
title('Target-set reach rate on the nonlinear plant');
exportgraphics(fig, fullfile(figDir, 'mc_reach_summary.png'), 'Resolution', 150);

fprintf('figures saved to %s\n', figDir);


function r = pad(r, frac)
    w = diff(r);
    if w <= 0, w = max(abs(r(1)), 1); end
    r = r + frac * w * [-1 1];
end
