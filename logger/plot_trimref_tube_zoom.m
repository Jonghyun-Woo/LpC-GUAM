function fig = plot_trimref_tube_zoom(traj, brt, gv, k, opts)
% PLOT_TRIMREF_TUBE_ZOOM (u, t, theta) BRT tube for ONE segment.
%
%   fig = plot_trimref_tube_zoom(traj, brt, gv, k, opts)
%
% Inputs
%   traj : TrimScheduleTrajectory.build() output
%   brt  : 1 x n_trim cell of 4D BRT value arrays (u,w,q,theta)
%   gv   : 1x4 cell of BRT grid vectors
%   k    : segment index in [1, n_trim-1] (trim k -> k+1)
%   opts : optional struct
%            .zoom_dt : slice spacing [s]                     (default 0.2)
%            .guam    : flown GUAM trace overlay (red line), struct with
%                         .time (1xM) [s], .lon (4xM) [u; w; q; theta]
%
% Fine time slices over the segment's scheduled window [t_node(k),
% t_node(k+1)]: current anchor BRT_k colored blue->yellow by time, next
% anchor BRT_{k+1} in red; black = reference trajectory, red = flown GUAM.
% Slices are taken at the REFERENCE's own (w,q) perturbation per anchor.
if nargin < 5, opts = struct(); end
zoom_dt = getfield_default(opts, 'zoom_dt', 0.2);
guam    = getfield_default(opts, 'guam', []);

t   = traj.time;  lon = traj.lon;
kn  = k + 1;
t0  = traj.t_node(k);  t1 = traj.t_node(kn);
cmap = parula(256);

fig = figure('Name', sprintf('trimref vs BRT: segment %d tube zoom', k), ...
             'Position', [90 90 1250 800]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on');

t_sl = t0:zoom_dt:t1;
if abs(t_sl(end) - t1) > 1e-9, t_sl(end + 1) = t1; end
for ts = t_sl
    [~, i] = min(abs(t - ts));
    fr  = (ts - t0) / (t1 - t0);
    col = cmap(max(1, round(fr * 255) + 1), :);
    C = BRTSlice.zero_contour(brt{k}, gv, lon(:, i), traj.trim_lon(:, k));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', col, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    C = BRTSlice.zero_contour(brt{kn}, gv, lon(:, i), traj.trim_lon(:, kn));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', [0.85 0.1 0.1 0.45], 'LineWidth', 0.9, ...
              'HandleVisibility', 'off');
    end
end

m = t >= t0 & t <= t1;
plot3(ax, lon(1, m), t(m), rad2deg(lon(4, m)), 'k-', 'LineWidth', 2.5, ...
      'DisplayName', 'reference trajectory');
plot3(ax, traj.trim_lon(1, k:kn), [t0 t1], rad2deg(traj.trim_lon(4, k:kn)), ...
      'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
      'DisplayName', 'trim points (BRT target origins)');
if ~isempty(guam)
    mg = guam.time >= t0 & guam.time <= t1;
    plot3(ax, guam.lon(1, mg), guam.time(mg), rad2deg(guam.lon(4, mg)), ...
          '-', 'Color', [1 0 0], 'LineWidth', 2.0, ...
          'DisplayName', 'GUAM flown trajectory');
end
legend(ax, 'Location', 'northeast');
xlabel(ax, 'u [ft/s]'); ylabel(ax, 'time [s]'); zlabel(ax, '\theta [deg]');
title(ax, sprintf(['Segment %d (UH_{%d} \\rightarrow UH_{%d}) tube, ' ...
      '\\Deltat = %.2g s — colored: BRT_{%d} (current), red: BRT_{%d} (next)'], ...
      k, k, kn, zoom_dt, k, kn));
view(ax, -60, 22);
end
