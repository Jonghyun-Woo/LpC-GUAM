function figs = plot_rg_run(traj, brtV, L, opts)
% PLOT_RG_RUN figures for a reference-governor run (TrimRefGovernor).
%
%   figs = plot_rg_run(traj, brtV, L, opts)
%
% opts.seg      : transition segment to zoom into, 1..n_trim-1. Segment k runs
%                 from trim point k to trim point k+1, timed on the FLOWN
%                 command (when v actually reached each trim speed), not on
%                 the plan. Takes priority over win_t/win_half.
% opts.seg_pad  : margin added either side of a segment [s] (default 0.5)
% opts.win_half : half-width of the zoom window [s]        (default 2)
% opts.win_t    : centre of the zoom window [s]            (default: the
%                 most-constrained solve, i.e. max V_hold)
% opts.n_slices : tube slices in the zoom figure           (default 12)
%
% The whole-mission figure takes no slice count: it draws one BRT per trim
% point, at the instant the applied command first reaches that trim speed.
%
% Red always means "the schedule was not followed":
%   red area between r and v : how much of the plan is outstanding
%   red ticks along the top  : instants at which r was judged inadmissible
%   red shading (zoom only)  : spans where kappa < 1
%
% Note on the ticks: with a 10 Hz governor and a plan the vehicle cannot keep
% up with, r is rejected at almost every solve. Drawing one full-height line
% per rejection paints the whole axis, so figure 1 shows them as a density
% rug at the top edge and reserves full lines for the zoom figures, where
% they are countable.
%
% Figures
%   figs.main     speeds / liveness / commands + kappa
%   figs.mission  (u,t,theta) over the whole run: plan, command, flown state
%   figs.zoom     (u,t,theta) tube + command timeline + predictions
%   figs.decision the line search itself: V versus kappa, solve by solve
if nargin < 4, opts = struct(); end
n_slices = getfield_default(opts, 'n_slices', 12);
win_half = getfield_default(opts, 'win_half', 2.0);

t   = L.t;
follow = L.kappa >= 1 - 1e-9;
[sIdx, eIdx] = runs_of(~follow);
t_rej = L.dec_t(L.dec_reject);

% When the applied command first reached each trim speed. Segment boundaries
% are taken from this, so the window follows the vehicle rather than the plan.
[t_cross, i_cross] = trim_cross_times(brtV, L);

% Zoom window. Priority: opts.seg > opts.win_t > the most-constrained solve.
if isfield(opts, 'seg') && ~isempty(opts.seg)
    kseg = opts.seg;
    n_seg = numel(t_cross) - 1;
    assert(kseg >= 1 && kseg <= n_seg, ...
           'opts.seg must be 1..%d (there are %d transition segments).', ...
           n_seg, n_seg);
    assert(~isnan(t_cross(kseg+1)), ...
           ['segment %d was never completed in this run - the command only ' ...
            'reached trim point %d.'], kseg, find(~isnan(t_cross), 1, 'last'));
    pad = getfield_default(opts, 'seg_pad', 0.5);
    tw  = [max(t(1),   t_cross(kseg)   - pad), ...
           min(t(end), t_cross(kseg+1) + pad)];
    win_label = sprintf('segment %d  (trim %d \\rightarrow %d,  %.1f \\rightarrow %.1f ft/s)', ...
                        kseg, kseg, kseg+1, ...
                        brtV.u_anchor(kseg), brtV.u_anchor(kseg+1));
else
    if isfield(opts, 'win_t') && ~isempty(opts.win_t)
        tc = opts.win_t;
    else
        vh = L.dec_Vh;  vh(isnan(vh)) = -inf;
        [~, jc] = max(vh);
        tc = L.dec_t(jc);
    end
    tw = [max(t(1), tc - win_half), min(t(end), tc + win_half)];
    win_label = '';
end
fprintf('zoom window: %.2f - %.2f s', tw(1), tw(2));
if ~isempty(win_label)
    fprintf('   %s', strrep(win_label, '\rightarrow', '->'));
end
fprintf('\n');

figs = struct();

%% ------------------------------------------------------------ Figure 1
f1 = figure('Name', 'RG governor: overview', 'Position', [60 40 1250 880]);
ax = gobjects(1,3);

ax(1) = subplot(3,1,1); hold on; grid on;
plot(ax(1), traj.time, traj.lon(1,:), 'k--', 'LineWidth', 1.0, ...
     'DisplayName', 'scheduled command r(t)');
plot(ax(1), t, L.v, 'b-', 'LineWidth', 1.6, 'DisplayName', 'applied command v(t)');
plot(ax(1), t, L.st(4,:), '-', 'Color', [0.85 0.2 0.2 0.6], 'LineWidth', 1.2, ...
     'DisplayName', 'flown u');
ylabel(ax(1), 'speed [ft/s]'); legend(ax(1), 'Location', 'southeast');
title(ax(1), sprintf(['Reference governor: final command at %.2f s ' ...
      '(the plan wanted %.2f s)'], L.t_done, traj.t_node(end)));

ax(2) = subplot(3,1,2); hold on; grid on;
plot(ax(2), t, L.Vnow, 'Color', [0 0.35 0.75], 'LineWidth', 1.4, ...
     'DisplayName', 'V(x, v) at the applied command');
yline(ax(2), 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylabel(ax(2), 'V'); ylim(ax(2), [-0.5 0.3]);
legend(ax(2), 'Location', 'southeast');
title(ax(2), 'Liveness actually achieved (V \geq 0 would mean the certificate is lost)');

% --- bottom: both commands, the outstanding plan, and kappa -------------
ax(3) = subplot(3,1,3); hold on; grid on;
yyaxis(ax(3), 'left');
% the gap between plan and applied command, as an area
fill(ax(3), [t, fliplr(t)], [L.r, fliplr(L.v)], [0.85 0.1 0.1], ...
     'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
     'DisplayName', 'plan not yet followed (r - v)');
plot(ax(3), t, L.r, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0, ...
     'DisplayName', 'scheduled r(t)');
plot(ax(3), t, L.v, '-', 'Color', [0 0.2 0.8], 'LineWidth', 1.5, ...
     'DisplayName', 'applied v(t)');
ylabel(ax(3), 'command [ft/s]');
ylL = [-8, max(L.r)*1.25];  ylim(ax(3), ylL);
% rejection instants as a density rug along the top edge
y0 = ylL(2)*0.955;  y1 = ylL(2)*0.995;
for j = 1:numel(t_rej)
    plot(ax(3), [t_rej(j) t_rej(j)], [y0 y1], '-', ...
         'Color', [0.85 0.1 0.1 0.55], 'LineWidth', 0.6, ...
         'HandleVisibility', 'off');
end
plot(ax(3), nan, nan, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.0, ...
     'DisplayName', 'r judged inadmissible (tick)');

yyaxis(ax(3), 'right');
stairs(ax(3), t, L.kappa, '-', 'Color', [0 0.55 0.2], 'LineWidth', 1.3, ...
       'DisplayName', '\kappa');
ylabel(ax(3), '\kappa'); ylim(ax(3), [-0.05 1.6]);
xlabel(ax(3), 'time [s]');
legend(ax(3), 'Location', 'northwest', 'NumColumns', 2, 'FontSize', 8);
title(ax(3), sprintf(['Governor action: r rejected at %d of %d solves; ' ...
      'schedule not followed for %.1f %% of the mission'], ...
      nnz(L.dec_reject), numel(L.dec_t), 100*nnz(~follow)/numel(t)));
linkaxes(ax, 'x');  xlim(ax(1), [0 t(end)]);
figs.main = f1;

%% ------------------------------------------------------------ Figure 2
% Same (u, t, theta) space as the zoom figure, but the whole mission: what
% the plan asked for, what the governor let through, and where the vehicle
% actually went. The three curves live in the same axes so the gap between
% them IS the governor's action.
f2 = figure('Name', 'RG governor: full mission', 'Position', [70 50 1150 800]);
mx = axes(f2); hold(mx, 'on'); grid(mx, 'on');

% One slice per trim point, drawn at the instant the applied command first
% reaches that trim speed. At those instants v sits exactly on a trim point,
% so each contour is the zero level of ONE stored BRT table rather than a
% blend of two - the 20 curves correspond one-to-one with the 20 tables.
% (The w,q the table is sliced at still comes from the flown state.)
keep = ~isnan(t_cross);
t_slice = t_cross(keep);  i_slice = i_cross(keep);

cmapM = parula(256);
for k = 1:numel(t_slice)
    colM = cmapM(round((k-1)/max(numel(t_slice)-1,1)*255)+1, :);
    C = brt_zero_contour(brtV, L.v(i_slice(k)), L.st([4 6 11 8], i_slice(k)));
    for c = 1:numel(C)
        p = C{c};
        plot3(mx, p(1,:), t_slice(k)*ones(1,size(p,2)), p(2,:), ...
              'Color', [colM 0.45], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
end
plot3(mx, nan, nan, nan, '-', 'Color', [0.35 0.55 0.55], 'LineWidth', 1.0, ...
      'DisplayName', sprintf('BRT of each trim point (%d of %d reached)', ...
                             numel(t_slice), numel(t_cross)));

% the three curves, all taken through the same trim map x_e(.)
xr = zeros(2, numel(t));  xv = zeros(2, numel(t));
for i = 1:numel(t)
    er = brtV.trim_at(L.r(i));  xr(:,i) = [er(1); rad2deg(er(4))];
    ev = brtV.trim_at(L.v(i));  xv(:,i) = [ev(1); rad2deg(ev(4))];
end
plot3(mx, xr(1,:), t, xr(2,:), '--', 'Color', [0.15 0.15 0.15], ...
      'LineWidth', 1.8, 'DisplayName', 'scheduled trajectory  x_e(r)');
plot3(mx, xv(1,:), t, xv(2,:), '-', 'Color', [0 0.2 0.85], ...
      'LineWidth', 2.4, 'DisplayName', 'governed command  x_e(v)');
plot3(mx, L.st(4,:), t, rad2deg(L.st(8,:)), '-', 'Color', [0.85 0.2 0.2], ...
      'LineWidth', 2.0, 'DisplayName', 'flown state  (u, \theta)');

xlabel(mx, 'u [ft/s]'); ylabel(mx, 'time [s]'); zlabel(mx, '\theta [deg]');
title(mx, sprintf(['Whole mission: the plan reached the final command at %.2f s, ' ...
      'the governor at %.2f s'], traj.t_node(end), L.t_done));
legend(mx, 'Location', 'northeast'); view(mx, -60, 22);
ylim(mx, [t(1) t(end)]);
figs.mission = f2;

%% ------------------------------------------------------------ Figure 3
w1 = find(t >= tw(1), 1, 'first');
w2 = find(t <= tw(2), 1, 'last');

f3 = figure('Name', 'RG governor: zoom', 'Position', [80 30 1250 900]);
bx = gobjects(1,3);

bx(1) = subplot(3,1,1); hold on; grid on;
cmap = parula(256);
for ts = linspace(tw(1), tw(2), n_slices)
    [~, i] = min(abs(t - ts));
    col = cmap(max(1, round((ts-tw(1))/max(diff(tw),eps)*255)+1), :);
    C = brt_zero_contour(brtV, L.v(i), L.st([4 6 11 8], i));
    for c = 1:numel(C)
        p = C{c};
        plot3(bx(1), p(1,:), ts*ones(1,size(p,2)), p(2,:), 'Color', col, ...
              'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
end
% Only two curves here, in the same colours as the whole-mission figure:
% blue = what the governor let through, red = flown.
%
% The plan x_e(r) is deliberately NOT drawn in this panel. The tube shown is
% the BRT of the APPLIED command, centred on x_e(v); x_e(r) is a different
% trim's equilibrium, so "is it inside this tube" is not the safety question
% and reads as an alarm when it is not. Worse, segment windows are timed on
% the flown command, and the plan runs ahead of it, so within one window the
% plan sits at a quite different speed and just stretches the u axis. The
% question it does answer - "would the schedule itself be admissible" - is
% V_sched, already plotted as red stairs in the bottom panel.
xv_w = zeros(2, w2-w1+1);  q = 0;
for i = w1:w2
    q = q + 1;
    ev = brtV.trim_at(L.v(i));  xv_w(:,q) = [ev(1); rad2deg(ev(4))];
end
plot3(bx(1), xv_w(1,:), t(w1:w2), xv_w(2,:), '-', 'Color', [0 0.2 0.85], ...
      'LineWidth', 2.4, 'DisplayName', 'governed command  x_e(v)');
plot3(bx(1), L.st(4,w1:w2), t(w1:w2), rad2deg(L.st(8,w1:w2)), '-', ...
      'Color', [0.85 0.2 0.2], 'LineWidth', 2.0, ...
      'DisplayName', 'flown state  (u, \theta)');
xlabel(bx(1), 'u [ft/s]'); ylabel(bx(1), 'time [s]'); zlabel(bx(1), '\theta [deg]');
if isempty(win_label)
    title(bx(1), sprintf('BRT of the commanded trim, %.1f-%.1f s', tw(1), tw(2)));
else
    title(bx(1), sprintf('BRT of the commanded trim - %s,  %.1f-%.1f s', ...
                         win_label, tw(1), tw(2)));
end
legend(bx(1), 'Location', 'northeast'); view(bx(1), -60, 22);

bx(2) = subplot(3,1,2); hold on; grid on;
ylC = [min(L.v(w1:w2))-3, max(L.r(w1:w2))+3];
shade_spans(bx(2), t, sIdx, eIdx, ylC, [0.85 0.1 0.1], 0.10);
for j = 1:numel(t_rej)
    if t_rej(j) >= tw(1) && t_rej(j) <= tw(2)
        xline(bx(2), t_rej(j), '-', 'Color', [0.85 0.1 0.1 0.55], ...
              'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
plot(bx(2), t, L.r, 'k--', 'LineWidth', 1.2, 'DisplayName', 'scheduled r(t)');
plot(bx(2), t, L.v, 'b-',  'LineWidth', 1.8, 'DisplayName', 'applied v(t)');
plot(bx(2), t, L.st(4,:), '-', 'Color', [0.85 0.2 0.2 0.7], 'LineWidth', 1.3, ...
     'DisplayName', 'flown u');
ylim(bx(2), ylC); ylabel(bx(2), 'speed [ft/s]');
legend(bx(2), 'Location', 'southeast');
title(bx(2), 'Command timeline (red line = r rejected, red band = \kappa < 1)');

bx(3) = subplot(3,1,3); hold on; grid on;
inw  = L.dec_t >= tw(1) & L.dec_t <= tw(2);
vtop = max([0.3, 1.15*max(L.dec_Vs(inw & ~isnan(L.dec_Vs)))]);
ylP  = [min(-0.4, 1.15*min(L.Vnow(w1:w2))), vtop];
shade_spans(bx(3), t, sIdx, eIdx, ylP, [0.85 0.1 0.1], 0.10);
plot(bx(3), t, L.Vnow, '-', 'Color', [0 0.35 0.75], 'LineWidth', 1.6, ...
     'DisplayName', 'V achieved by the applied command');
stairs(bx(3), L.dec_t, L.dec_Vs, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5, ...
       'DisplayName', 'predicted V if r were applied');
stairs(bx(3), L.dec_t, L.dec_Vh, '-', 'Color', [0 0.55 0.2], 'LineWidth', 1.5, ...
       'DisplayName', 'predicted V if v(t-1) were held');
yline(bx(3), 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylim(bx(3), ylP); xlabel(bx(3), 'time [s]'); ylabel(bx(3), 'V');
legend(bx(3), 'Location', 'northeast');
title(bx(3), 'The predictions behind each decision (above 0 = rejected)');
linkaxes(bx(2:3), 'x');  xlim(bx(2), tw);
figs.zoom = f3;

%% ------------------------------------------------------------ Figure 4
% The line search itself: for each solve in the window, V against kappa.
sel = find(L.dec_t >= tw(1) & L.dec_t <= tw(2));
if numel(sel) > 8, sel = sel(round(linspace(1, numel(sel), 8))); end

f4 = figure('Name', 'RG governor: the line search', 'Position', [100 60 1200 760]);
cx = gobjects(1,2);

cx(1) = subplot(2,1,1); hold on; grid on;
col = turbo(max(numel(sel),1));
for a = 1:numel(sel)
    j = sel(a);  C = L.dec_cands{j};
    if isempty(C), continue; end
    [ks, ord] = sort(C(1,:));  Vs = C(2, ord);
    plot(cx(1), ks, Vs, '-o', 'Color', col(a,:), 'LineWidth', 1.4, ...
         'MarkerSize', 4, 'DisplayName', sprintf('t = %.2f s', L.dec_t(j)));
    [~, ib] = min(abs(ks - L.dec_k(j)));
    plot(cx(1), ks(ib), Vs(ib), 'p', 'Color', col(a,:), ...
         'MarkerFaceColor', col(a,:), 'MarkerSize', 13, 'HandleVisibility', 'off');
end
yline(cx(1), 0, 'k-', 'LineWidth', 1.4, 'HandleVisibility', 'off');
xlabel(cx(1), '\kappa candidate'); ylabel(cx(1), 'predicted V');
allV = cell2mat(cellfun(@(C) C(2,:), L.dec_cands(sel), 'UniformOutput', false));
ylim(cx(1), [min(-0.4, 1.15*min(allV)), max(0.3, 1.05*max(allV))]);
xlim(cx(1), [-0.05 1.05]);
legend(cx(1), 'Location', 'northwest', 'FontSize', 8);
title(cx(1), ['Line search: predicted V against \kappa. Star = the value applied. ' ...
              'Only points below the black line are admissible.']);

cx(2) = subplot(2,1,2); hold on; grid on;
yyaxis(cx(2), 'left');
plot(cx(2), t, L.r - L.v, '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5, ...
     'DisplayName', 'r - v  (plan outstanding)');
ylabel(cx(2), 'r - v [ft/s]');
yyaxis(cx(2), 'right');
stairs(cx(2), L.dec_t, L.dec_k, '-', 'Color', [0 0.35 0.75], 'LineWidth', 1.6, ...
       'DisplayName', 'applied \kappa');
ylabel(cx(2), '\kappa'); ylim(cx(2), [-0.05 1.05]);
xlim(cx(2), tw); xlabel(cx(2), 'time [s]');
legend(cx(2), 'Location', 'northeast');
title(cx(2), 'What the search produced over the same window');
figs.decision = f4;
end

% =========================================================================
function [tc, ic] = trim_cross_times(brtV, L)
% Time (and sample index) at which the applied command first reached each trim
% speed. NaN for trim points the run never got to. With the governor active
% kappa >= 0, so v is non-decreasing and each speed is crossed exactly once.
u_a = brtV.u_anchor;
tc = nan(1, numel(u_a));  ic = nan(1, numel(u_a));
for k = 1:numel(u_a)
    i = find(L.v >= u_a(k) - 1e-9, 1, 'first');
    if ~isempty(i), tc(k) = L.t(i);  ic(k) = i; end
end
end

function [s, e] = runs_of(mask)
d = diff([false, mask(:)', false]);
s = find(d == 1);  e = find(d == -1) - 1;
end

function shade_spans(ax, t, sIdx, eIdx, yl, col, alpha)
for j = 1:numel(sIdx)
    patch(ax, t([sIdx(j) eIdx(j) eIdx(j) sIdx(j)]), ...
          [yl(1) yl(1) yl(2) yl(2)], col, 'FaceAlpha', alpha, ...
          'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end

% brt_zero_contour now lives in its own file so the option-1 (dV/dt) plots can
% use the same tube slicing.
