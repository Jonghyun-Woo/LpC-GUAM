function fig = plot_trimref_handover(traj, ev, brt, gv, k, guam)
% PLOT_TRIMREF_HANDOVER (u,theta) handover slices for ONE segment.
%
%   fig = plot_trimref_handover(traj, ev, brt, gv, k, guam)
%
% Inputs
%   traj : TrimScheduleTrajectory.build() output
%   ev   : struct with V_cur, V_next (1xN) from check_trimref_liveness
%   brt  : 1 x n_trim cell of 4D BRT value arrays (u,w,q,theta)
%   gv   : 1x4 cell of BRT grid vectors
%   k    : segment index in [1, n_trim-1] (trim k -> k+1)
%   guam : optional flown GUAM trace overlay (red), struct with
%            .time (1xM) [s], .lon (4xM) [u; w; q; theta]
%
% Six subplots at fixed progress fractions of segment k: BRT_k (blue) vs
% BRT_{k+1} (red) sliced at the reference's own (w,q); black star = ref
% point, red line/dot = flown GUAM path over the segment / at that time.
if nargin < 6, guam = []; end

t   = traj.time;  lon = traj.lon;
kn  = k + 1;
fracs = [0 0.25 0.4 0.5 0.6 0.75];

fig = figure('Name', sprintf('trimref vs BRT: segment %d handover', k), ...
             'Position', [100 100 1400 800]);
m  = find(traj.seg == k & traj.progress < kn);
fr = traj.progress(m) - k;
for p = 1:numel(fracs)
    axp = subplot(2, 3, p); hold(axp, 'on'); grid(axp, 'on');
    [~, j] = min(abs(fr - fracs(p)));
    i = m(j);
    BRTSlice.draw(axp, brt{k},  gv, lon(:, i), traj.trim_lon(:, k),  [0 0.447 0.741]);
    BRTSlice.draw(axp, brt{kn}, gv, lon(:, i), traj.trim_lon(:, kn), [0.85 0.1 0.1]);
    plot(axp, lon(1, m), rad2deg(lon(4, m)), 'k-', 'LineWidth', 0.8);
    if ~isempty(guam)
        mg = guam.time >= traj.t_node(k) & guam.time <= traj.t_node(kn);
        plot(axp, guam.lon(1, mg), rad2deg(guam.lon(4, mg)), '-', ...
             'Color', [1 0 0], 'LineWidth', 1.8);
        [~, jg] = min(abs(guam.time - t(i)));
        plot(axp, guam.lon(1, jg), rad2deg(guam.lon(4, jg)), 'o', ...
             'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 7);
    end
    plot(axp, lon(1, i), rad2deg(lon(4, i)), 'kp', 'MarkerFaceColor', 'k', ...
         'MarkerSize', 14);
    title(axp, sprintf('progress %.0f%%  (V_{cur}=%+.3f, V_{next}=%+.3f)', ...
          100 * fracs(p), ev.V_cur(i), ev.V_next(i)));
    xlabel(axp, 'u [ft/s]'); ylabel(axp, '\theta [deg]');

    if p == 1
        % Legend via proxy handles (drawn objects carry no DisplayName).
        hp = gobjects(1, 0);
        hp(end+1) = patch(axp, 'XData', nan, 'YData', nan, 'FaceColor', [0 0.447 0.741], ...
            'FaceAlpha', 0.15, 'EdgeColor', [0 0.447 0.741], 'DisplayName', sprintf('BRT_{%d} (current)', k)); %#ok<AGROW>
        hp(end+1) = patch(axp, 'XData', nan, 'YData', nan, 'FaceColor', [0.85 0.1 0.1], ...
            'FaceAlpha', 0.15, 'EdgeColor', [0.85 0.1 0.1], 'DisplayName', sprintf('BRT_{%d} (next)', kn)); %#ok<AGROW>
        hp(end+1) = plot(axp, nan, nan, '^', 'Color', [0 0.447 0.741], 'MarkerFaceColor', ...
            [0 0.447 0.741], 'DisplayName', 'trim k / k+1 (BRT target origins)'); %#ok<AGROW>
        hp(end+1) = plot(axp, nan, nan, 'k-', 'DisplayName', 'reference path (this segment)'); %#ok<AGROW>
        hp(end+1) = plot(axp, nan, nan, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12, ...
            'DisplayName', 'reference point @ this progress'); %#ok<AGROW>
        if ~isempty(guam)
            hp(end+1) = plot(axp, nan, nan, '-', 'Color', [1 0 0], 'LineWidth', 1.8, ...
                'DisplayName', 'GUAM flown path (this segment)'); %#ok<AGROW>
            hp(end+1) = plot(axp, nan, nan, 'o', 'Color', [1 0 0], 'MarkerFaceColor', ...
                [1 0 0], 'DisplayName', 'GUAM position @ this progress'); %#ok<AGROW>
        end
        legend(axp, hp, 'Location', 'southeast', 'FontSize', 7);
    end
end
sgtitle(fig, sprintf(['Segment %d handover: BRT_{%d} (blue) vs BRT_{%d} (red), ' ...
        'sliced at trajectory (w,q); star = ref, red = flown GUAM'], k, k, kn));
end
