function fig = plot_governor_run(traj, brt, gv, G, k, opts)
% PLOT_GOVERNOR_RUN governed-run figure: segment tube + command timeline.
%
%   fig = plot_governor_run(traj, brt, gv, G, k, opts)
%
% Inputs
%   traj : TrimScheduleTrajectory.build() output (NOMINAL schedule)
%   brt  : 1 x n_trim cell of 4D BRT value arrays
%   gv   : 1x4 cell of BRT grid vectors
%   G    : governed run log (run_governor_option1) with fields
%            .t (1xN), .st (12xN), .kappa (1xN), .progress (1xN), .u_ref (1xN)
%   k    : segment index in [1, n_trim-1] for the top tube view
%   opts : optional struct
%            .n_slices : tube slices over the segment window (default 12)
%
% Top   : (u, t, theta) BRT_k / BRT_{k+1} tube over the GOVERNED window of
%         segment k, with the governed command (black, on the trim line)
%         and the flown GUAM state (red).
% Middle: governed command u_ref vs time; kappa = 0 (governor intervening)
%         shaded translucent green; nominal T_seg node times as vertical
%         dashed lines - the horizontal offset of the governed schedule
%         from those lines is exactly the accumulated governor hold time.
% Bottom: the PREDICTION (drawn only when G.V_hat is logged). V_next is
%         what the CBF sees now; V_hat is the forecast at the handover
%         deadline. V_hat crossing zero BEFORE V_next does is the lead
%         time the predictive gate buys.
if nargin < 6, opts = struct(); end
n_slices = getfield_default(opts, 'n_slices', 12);

n_tr = traj.n_trim;
kn   = k + 1;
dt   = G.t(2) - G.t(1);

% governed command point on the trim line (u_ref, w_trim, q_trim, th_ref)
ks    = min(floor(G.progress), n_tr - 1);
frc   = G.progress - ks;
th_ref = traj.trim_lon(4, ks) + frc .* (traj.trim_lon(4, ks + 1) - traj.trim_lon(4, ks));
w_trim = traj.trim_lon(2, 1);          % identical for the whole WH family
q_trim = traj.trim_lon(3, 1);

% governed time window of segment k: node-k crossing -> node-(k+1) crossing
if k == 1
    i0 = find(G.t >= traj.t_node(1), 1, 'first');
else
    i0 = find(G.progress >= k, 1, 'first');
end
i1 = find(G.progress >= kn, 1, 'first');
if isempty(i1), i1 = numel(G.t); end
t0 = G.t(i0);  t1 = G.t(i1);

has_pred = isfield(G, 'V_hat') && any(G.V_hat ~= 0);
nrow = 2 + has_pred;

fig = figure('Name', sprintf('governor run: segment %d', k), ...
             'Position', [70 70 1200 300 * nrow]);

% ------------------------------------------------------------- top: tube
ax1 = subplot(nrow, 1, 1); hold(ax1, 'on'); grid(ax1, 'on');
cmap = parula(256);
for ts = linspace(t0, t1, n_slices)
    [~, i] = min(abs(G.t - ts));
    cmd_pt = [G.u_ref(i); w_trim; q_trim; th_ref(i)];
    col = cmap(max(1, round((ts - t0) / max(t1 - t0, eps) * 255) + 1), :);
    C = BRTSlice.zero_contour(brt{k}, gv, cmd_pt, traj.trim_lon(:, k));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax1, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', col, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
    C = BRTSlice.zero_contour(brt{kn}, gv, cmd_pt, traj.trim_lon(:, kn));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax1, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', [0.85 0.1 0.1 0.45], 'LineWidth', 0.9, ...
              'HandleVisibility', 'off');
    end
end
m = i0:i1;
plot3(ax1, G.u_ref(m), G.t(m), rad2deg(th_ref(m)), 'k-', 'LineWidth', 2.5, ...
      'DisplayName', 'governed command (trim line)');
plot3(ax1, G.st(4, m), G.t(m), rad2deg(G.st(8, m)), '-', 'Color', [1 0 0], ...
      'LineWidth', 2.0, 'DisplayName', 'GUAM flown trajectory');
plot3(ax1, traj.trim_lon(1, k:kn), [t0 t1], rad2deg(traj.trim_lon(4, k:kn)), ...
      'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
      'DisplayName', 'trim points (BRT target origins)');
legend(ax1, 'Location', 'northeast');
xlabel(ax1, 'u [ft/s]'); ylabel(ax1, 'time [s]'); zlabel(ax1, '\theta [deg]');
title(ax1, sprintf(['Governed segment %d (UH_{%d} \\rightarrow UH_{%d}): ' ...
      'colored BRT_{%d}, red BRT_{%d}, window %.2f-%.2f s'], ...
      k, k, kn, k, kn, t0, t1));
view(ax1, -60, 22);

% ------------------------------------------------- middle: command timeline
ax2 = subplot(nrow, 1, 2); hold(ax2, 'on'); grid(ax2, 'on');
yl = [-5, max(G.u_ref) + 15];

% translucent green background where the governor intervened (kappa = 0)
hold_mask = (G.kappa == 0) & (G.t >= traj.t_node(1));
d = diff([0, hold_mask, 0]);
sIdx = find(d == 1);  eIdx = find(d == -1) - 1;
for j = 1:numel(sIdx)
    patch(ax2, G.t([sIdx(j) eIdx(j) eIdx(j) sIdx(j)]), ...
          [yl(1) yl(1) yl(2) yl(2)], [0 0.7 0], 'FaceAlpha', 0.15, ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% nominal T_seg boundaries (schedule without governor)
for j = 1:n_tr
    xline(ax2, traj.t_node(j), 'k--', 'Alpha', 0.4, 'HandleVisibility', 'off');
end

plot(ax2, traj.time, traj.lon(1, :), 'k--', 'LineWidth', 1.0, ...
     'DisplayName', 'nominal schedule u_{ref} (T_{seg} grid, dashed lines)');
plot(ax2, G.t, G.u_ref, 'b-', 'LineWidth', 1.6, ...
     'DisplayName', 'governed command u_{ref}');
plot(ax2, G.t, G.st(4, :), '-', 'Color', [1 0 0 0.6], 'LineWidth', 1.2, ...
     'DisplayName', 'flown u');
patch(ax2, nan, nan, [0 0.7 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
      'DisplayName', 'governor active (\kappa = 0)');

t_hold = nnz(hold_mask) * dt;
i_fin  = find(G.progress >= n_tr, 1, 'first');
if isempty(i_fin), t_fin = G.t(end); else, t_fin = G.t(i_fin); end
ylim(ax2, yl);
ylabel(ax2, 'u_{ref} [ft/s]');
title(ax2, sprintf(['Command timeline - final trim reached %.2f s vs ' ...
      'nominal %.2f s: +%.2f s = governor hold %.2f s'], ...
      t_fin, traj.t_node(end), t_fin - traj.t_node(end), t_hold));
legend(ax2, 'Location', 'southeast');
if ~has_pred
    xlabel(ax2, 'time [s]');
    return;
end

% ------------------------------------------------------ bottom: prediction
ax3 = subplot(nrow, 1, 3); hold(ax3, 'on'); grid(ax3, 'on');
ylp = [-0.5, 0.5];
for j = 1:numel(sIdx)
    patch(ax3, G.t([sIdx(j) eIdx(j) eIdx(j) sIdx(j)]), ...
          [ylp(1) ylp(1) ylp(2) ylp(2)], [0 0.7 0], 'FaceAlpha', 0.15, ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
end
for j = 1:n_tr
    xline(ax3, traj.t_node(j), 'k--', 'Alpha', 0.4, 'HandleVisibility', 'off');
end
plot(ax3, G.t, G.V_next, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5, ...
     'DisplayName', 'V_{next} (measured now)');
plot(ax3, G.t, G.V_hat, '-', 'Color', [0.2 0.2 0.85], 'LineWidth', 1.3, ...
     'DisplayName', 'V_{hat} (forecast at deadline)');
yline(ax3, 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylim(ax3, ylp);
xlabel(ax3, 'time [s]'); ylabel(ax3, 'BRT value');
title(ax3, ['Prediction: V_{hat} = V_{next} + dV/dt\cdot(t_{rem} \mp t_{mar}). ' ...
            'V_{hat} \geq 0 while V_{next} < 0 is the early warning']);
legend(ax3, 'Location', 'southeast');
linkaxes([ax2 ax3], 'x');
end
