function figs = plot_trimref_brt(traj, ev, brt, gv, opts)
% PLOT_TRIMREF_BRT visualize a trim-schedule reference trajectory vs LON BRTs.
%
%   figs = plot_trimref_brt(traj, ev, brt, gv, opts)
%
% Inputs
%   traj : TrimScheduleTrajectory.build() output
%   ev   : struct with V_cur, V_next (1xN) from check_trimref_liveness
%   brt  : 1 x n_trim cell of 4D BRT value arrays (u,w,q,theta)
%   gv   : 1x4 cell of BRT grid vectors
%   opts : optional struct
%            .hand_seg  : segment index for the handover figure  (default 1)
%            .zoom_seg  : segment index for the zoomed tube      (default 1)
%            .zoom_dt   : slice spacing inside the zoomed tube [s] (default 0.2)
%            .guam      : flown GUAM trace to overlay in red on the tube and
%                         handover figures: struct with
%                           .time (1xM) [s], .lon (4xM) [u; w; q; theta]
%
% Figures
%   figs.timeseries : V_cur / V_next along the trajectory + handover windows
%   figs.tube3d     : (u, t, theta) zero-contour tube, one slice per trim
%                     node (= segment start/end), sliced at the trajectory's
%                     own (w,q) perturbation
%   figs.tube3d_zoom: same tube zoomed into segment zoom_seg with fine
%                     time slices (current + next anchor)
%   figs.handover   : (u,theta) slices across one segment's handover
%
% NOTE (dimension squish): the BRT lives in (u,w,q,theta). Every 2D/3D view
% here slices the 4D value function AT THE TRAJECTORY'S OWN (w,q) point of
% the respective anchor, so inside/outside read directly off the plot. For
% the pure trim-schedule reference this slice is exact ((w,q) pert = 0);
% for a flown GUAM trace pass the actual states through the same code path.
if nargin < 5, opts = struct(); end
hand_seg = getfield_default(opts, 'hand_seg', 1);
zoom_seg = getfield_default(opts, 'zoom_seg', 1);
zoom_dt  = getfield_default(opts, 'zoom_dt', 0.2);
guam     = getfield_default(opts, 'guam', []);

t   = traj.time;   lon = traj.lon;   n_trim = traj.n_trim;
figs = struct();

% ---------------------------------------------------------------- fig 1
f1 = figure('Name', 'trimref vs BRT: V timeseries', 'Position', [60 60 1250 700]);
ax1 = subplot(4, 1, 1:3); hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, t, ev.V_cur,  'Color', [0 0.447 0.741], 'LineWidth', 1.6, ...
     'DisplayName', 'V (current-anchor BRT_k)');
plot(ax1, t, ev.V_next, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.6, ...
     'DisplayName', 'V (next-anchor BRT_{k+1})');
yline(ax1, 0, 'k-', 'HandleVisibility', 'off');
for k = 1:n_trim - 1
    m = traj.seg == k & traj.progress < n_trim;
    hand = m & ev.V_cur < 0 & ev.V_next < 0;
    if any(hand)
        th = t(hand);
        patch(ax1, [th(1) th(end) th(end) th(1)], ...
              [min(ylim(ax1))*[1 1] max(ylim(ax1))*[1 1]], ...
              [0 0.6 0], 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
    end
    xline(ax1, traj.t_node(k), 'Color', [0.6 0.6 0.6], 'Alpha', 0.5, ...
          'HandleVisibility', 'off');
end
xline(ax1, traj.t_node(end), 'Color', [0.6 0.6 0.6], 'Alpha', 0.5, ...
      'HandleVisibility', 'off');
ylabel(ax1, 'BRT value V   (V<0 = inside)');
title(ax1, ['Reference trajectory vs LON BRT value functions ' ...
            '(green = handover window: inside both BRT_k and BRT_{k+1})']);
legend(ax1, 'Location', 'northeast');

ax2 = subplot(4, 1, 4); grid(ax2, 'on');
plot(ax2, t, lon(1, :), 'k-', 'LineWidth', 1.4);
xlabel(ax2, 'time [s]'); ylabel(ax2, 'ref u [ft/s]');
linkaxes([ax1 ax2], 'x');
figs.timeseries = f1;

% ---------------------------------------------------------------- fig 2
% One slice per trim node: at t_node(k) the crossed anchor is BRT_k, which
% is simultaneously segment (k-1)'s end tube and segment k's start tube.
f2 = figure('Name', 'trimref vs BRT: (u,t,theta) tube', 'Position', [80 80 1250 800]);
ax = axes(f2); hold(ax, 'on'); grid(ax, 'on');
cmap = parula(256);
for k = 1:n_trim
    ts = traj.t_node(k);
    [~, i] = min(abs(t - ts));
    col = cmap(max(1, round(ts / t(end) * 255) + 1), :);
    C = BRTSlice.zero_contour(brt{k}, gv, lon(:, i), traj.trim_lon(:, k));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
plot3(ax, lon(1, :), t, rad2deg(lon(4, :)), 'k-', 'LineWidth', 2.5, ...
      'DisplayName', 'reference trajectory');
plot3(ax, traj.trim_lon(1, :), traj.t_node, rad2deg(traj.trim_lon(4, :)), ...
      'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, ...
      'DisplayName', 'trim points (BRT target origins)');
if ~isempty(guam)
    plot3(ax, guam.lon(1, :), guam.time, rad2deg(guam.lon(4, :)), '-', ...
          'Color', [1 0 0], 'LineWidth', 2.0, ...
          'DisplayName', 'GUAM flown trajectory');
    legend(ax, 'Location', 'northeast');
end
xlabel(ax, 'u [ft/s]'); ylabel(ax, 'time [s]'); zlabel(ax, '\theta [deg]');
title(ax, ['LON BRT zero-level tube, one slice per trim node ' ...
           '(= segment start/end), sliced at the trajectory''s (w,q)']);
view(ax, -60, 22);
figs.tube3d = f2;

% ---------------------------------------------------------------- fig 2b
figs.tube3d_zoom = plot_trimref_tube_zoom(traj, brt, gv, zoom_seg, ...
                       struct('guam', guam, 'zoom_dt', zoom_dt));

% ---------------------------------------------------------------- fig 3
figs.handover = plot_trimref_handover(traj, ev, brt, gv, hand_seg, guam);
end
