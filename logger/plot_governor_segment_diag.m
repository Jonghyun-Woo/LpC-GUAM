function fig = plot_governor_segment_diag(traj, brt, gv, G, k, opts)
% PLOT_GOVERNOR_SEGMENT_DIAG why the governor acted in one governed segment.
%
%   fig = plot_governor_segment_diag(traj, brt, gv, G, k, opts)
%
% Inputs
%   traj : TrimScheduleTrajectory.build() output (nominal schedule)
%   brt  : 1 x n_trim cell of 4D BRT value arrays
%   gv   : 1x4 cell of BRT grid vectors
%   G    : governed run log; needs .t .st(12xN) .kappa .progress
%   k    : segment index (1..n_trim-1); window = where floor(progress) == k
%   opts : .n_slices (tube slices, default 12)
%
% Panels
%   1) (u, t, theta) BRT tube over the GOVERNED window. The flown trajectory
%      is RED where the governor let the schedule run and GREEN where it
%      froze the command (kappa = 0).
%   2) The states that carry the momentum: u, w (linear, momentum = m*u,
%      m*w) and theta, q (attitude / pitch rate, angular momentum = Iyy*q).
%      Green bands mark the frozen spans.
%   3) ATTRIBUTION. dV/dt decomposed per state via the chain rule,
%
%          dV/dt = (dV/du)*udot + (dV/dw)*wdot + (dV/dq)*qdot + (dV/dth)*thdot
%
%      against the NEXT anchor BRT_{k+1}. Whichever curve is most positive
%      is the state whose motion is pushing the vehicle out of the BRT -
%      this is the panel that answers "which angle moved how".
if nargin < 6, opts = struct(); end
n_slices = getfield_default(opts, 'n_slices', 12);

kn = k + 1;
m  = find(floor(G.progress) == k);
if isempty(m)
    error('plot_governor_segment_diag:noSegment', ...
          'Segment %d never active in this run.', k);
end
t   = G.t(m);
lon = G.st([4 6 11 8], m);        % [u; w; q; theta]
frz = G.kappa(m) == 0;
dt  = G.t(2) - G.t(1);
t0  = t(1);  t1 = t(end);

fig = figure('Name', sprintf('governor diagnosis: segment %d', k), ...
             'Position', [50 30 1250 950]);

% ------------------------------------------------------------ 1) tube
ax1 = subplot(3, 1, 1); hold(ax1, 'on'); grid(ax1, 'on');
cmap = parula(256);
for ts = linspace(t0, t1, n_slices)
    [~, j] = min(abs(t - ts));
    col = cmap(max(1, round((ts - t0) / max(t1 - t0, eps) * 255) + 1), :);
    C = BRTSlice.zero_contour(brt{k}, gv, lon(:, j), traj.trim_lon(:, k));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax1, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    C = BRTSlice.zero_contour(brt{kn}, gv, lon(:, j), traj.trim_lon(:, kn));
    for c = 1:numel(C)
        v = C{c};
        plot3(ax1, v(1, :), ts * ones(1, size(v, 2)), v(2, :), ...
              'Color', [0.85 0.1 0.1 0.35], 'LineWidth', 0.8, ...
              'HandleVisibility', 'off');
    end
end
plot_runs(ax1, ~frz, lon(1, :), t, rad2deg(lon(4, :)), [1 0 0], ...
          'flown: governor idle');
plot_runs(ax1,  frz, lon(1, :), t, rad2deg(lon(4, :)), [0 0.65 0], ...
          'flown: governor FROZEN (\kappa=0)');
plot3(ax1, traj.trim_lon(1, k:kn), [t0 t1], rad2deg(traj.trim_lon(4, k:kn)), ...
      'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
      'DisplayName', 'trim points');
xlabel(ax1, 'u [ft/s]'); ylabel(ax1, 'time [s]'); zlabel(ax1, '\theta [deg]');
title(ax1, sprintf(['Segment %d (UH%d\\rightarrowUH%d) governed window ' ...
      '%.2f-%.2f s   |   colored BRT_%d, red BRT_%d'], k, k, kn, t0, t1, k, kn));
legend(ax1, 'Location', 'northeast');
view(ax1, -60, 22);

% ------------------------------------------------- 2) momentum-carrying states
% Data first, then shade using the resulting limits: a patch drawn before
% the plots would stretch the axis to the patch's own extent.
ax2 = subplot(3, 1, 2); hold(ax2, 'on'); grid(ax2, 'on');
yyaxis(ax2, 'left');
plot(ax2, t, lon(1, :), '-',  'LineWidth', 1.6, 'DisplayName', 'u [ft/s]');
plot(ax2, t, lon(2, :), '--', 'LineWidth', 1.4, 'DisplayName', 'w [ft/s]');
ylabel(ax2, 'u, w [ft/s]   (linear momentum = m\cdotu, m\cdotw)');
ylL = pad_lim([lon(1, :), lon(2, :)]);
ylim(ax2, ylL);
shade(ax2, t, frz, ylL);
yyaxis(ax2, 'right');
plot(ax2, t, rad2deg(lon(4, :)), '-',  'LineWidth', 1.6, 'DisplayName', '\theta [deg]');
plot(ax2, t, rad2deg(lon(3, :)), '--', 'LineWidth', 1.4, 'DisplayName', 'q [deg/s]');
ylabel(ax2, '\theta [deg],  q [deg/s]   (ang. momentum = I_{yy}\cdotq)');
ylim(ax2, pad_lim(rad2deg([lon(4, :), lon(3, :)])));
title(ax2, 'States carrying the momentum (green = governor frozen)');
legend(ax2, 'Location', 'best', 'NumColumns', 2);

% ---------------------------------------------------------- 3) attribution
% dV/dt = grad(V)' * xdot, split per state, against the NEXT anchor.
xdot = zeros(4, numel(m));
for d = 1:4
    xdot(d, :) = gradient(lon(d, :), dt);
end
gradV = zeros(4, numel(m));
for j = 1:numel(m)
    gradV(:, j) = brt_gradient(brt{kn}, gv, lon(:, j) - traj.trim_lon(:, kn));
end
contrib = gradV .* xdot;                       % 4 x N
contrib = movmean(contrib, round(0.15 / dt), 2);   % de-noise the derivative

ax3 = subplot(3, 1, 3); hold(ax3, 'on'); grid(ax3, 'on');
lbl = {'u', 'w', 'q', '\theta'};
col = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.93 0.69 0.13; 0.85 0.10 0.10];
for d = 1:4
    plot(ax3, t, contrib(d, :), '-', 'Color', col(d, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s contribution', lbl{d}));
end
plot(ax3, t, sum(contrib, 1), 'k-', 'LineWidth', 2.0, ...
     'DisplayName', 'total dV_{next}/dt');
yline(ax3, 0, 'k-', 'HandleVisibility', 'off');
% Scale to the contributions themselves, not to a single spike.
yl = max(0.05, 1.2 * prctile(abs(contrib(:)), 98));
ylim(ax3, [-yl, yl]);
shade(ax3, t, frz, [-yl, yl]);
xlabel(ax3, 'time [s]'); ylabel(ax3, 'dV_{next}/dt  [1/s]');
title(ax3, ['Attribution: which state''s motion pushes V toward the BRT ' ...
            'boundary (positive = outward)']);
legend(ax3, 'Location', 'best', 'NumColumns', 2);
linkaxes([ax2 ax3], 'x');  xlim(ax2, [t0 t1]);

% ---- worst-offender summary in the command window ----------------------
if any(frz)
    jf = find(frz, 1, 'first');
    [~, dmax] = max(contrib(:, jf));
    fprintf(['\nsegment %d: governor first froze at t = %.2f s; ' ...
             'largest outward push was from %s ' ...
             '(%+.3f /s of a total %+.3f /s)\n'], ...
            k, t(jf), lbl{dmax}, contrib(dmax, jf), sum(contrib(:, jf)));
    fprintf('   state then: u=%.1f ft/s  w=%.1f ft/s  q=%+.2f deg/s  theta=%+.2f deg\n', ...
            lon(1, jf), lon(2, jf), rad2deg(lon(3, jf)), rad2deg(lon(4, jf)));
end
end

% -------------------------------------------------------------------------
function plot_runs(ax, mask, x, y, z, col, name)
% plot3 contiguous runs of mask in one color, single legend entry
d = diff([false, mask, false]);
s = find(d == 1);  e = find(d == -1) - 1;
first = true;
for j = 1:numel(s)
    idx = s(j):e(j);
    if numel(idx) < 2, continue; end
    if first
        plot3(ax, x(idx), y(idx), z(idx), '-', 'Color', col, ...
              'LineWidth', 2.2, 'DisplayName', name);
        first = false;
    else
        plot3(ax, x(idx), y(idx), z(idx), '-', 'Color', col, ...
              'LineWidth', 2.2, 'HandleVisibility', 'off');
    end
end
if first   % mask never true: emit an invisible proxy so the legend is stable
    plot3(ax, nan, nan, nan, '-', 'Color', col, 'LineWidth', 2.2, ...
          'DisplayName', name);
end
end

function yl = pad_lim(v)
lo = min(v);  hi = max(v);
if hi - lo < 1e-9, hi = lo + 1; end
pad = 0.08 * (hi - lo);
yl = [lo - pad, hi + pad];
end

function shade(ax, t, mask, yl)
d = diff([false, mask, false]);
s = find(d == 1);  e = find(d == -1) - 1;
for j = 1:numel(s)
    patch(ax, t([s(j) e(j) e(j) s(j)]), [yl(1) yl(1) yl(2) yl(2)], ...
          [0 0.7 0], 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
end
end

function g = brt_gradient(data, gv, x_pert)
% Central-difference gradient of the interpolated value function at x_pert,
% using half the grid spacing as the step.
g = zeros(4, 1);
xc = x_pert(:)';
for d = 1:4
    xc(d) = min(max(xc(d), gv{d}(1)), gv{d}(end));
end
for d = 1:4
    h  = 0.5 * (gv{d}(2) - gv{d}(1));
    xp = xc;  xm = xc;
    xp(d) = min(xc(d) + h, gv{d}(end));
    xm(d) = max(xc(d) - h, gv{d}(1));
    den = xp(d) - xm(d);
    if den <= 0, g(d) = 0; continue; end
    Vp = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xp(1), xp(2), xp(3), xp(4), 'linear');
    Vm = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xm(1), xm(2), xm(3), xm(4), 'linear');
    g(d) = (Vp - Vm) / den;
end
end
