function fig = plot_cbf_detection(traj, brt, gv, guam_trace, opts)
% PLOT_CBF_DETECTION per-segment CBF risk detection on the flown trajectory.
%
%   fig = plot_cbf_detection(traj, brt, gv, guam_trace, opts)
%
% Inputs
%   traj       : TrimScheduleTrajectory.build() output
%   brt        : 1 x n_trim cell of 4D BRT value arrays
%   gv         : 1x4 cell of BRT grid vectors
%   guam_trace : flown trace, struct with .time (1xM), .lon (4xM)
%   opts       : optional struct
%                  .mark : segment indices to highlight (default [])
%
% Every sample of the flown trajectory is classified against the two BRTs
% scheduled for its segment, on a per-segment progress axis (0-100 %), so
% each segment is one legible row instead of one tooth of a sawtooth.
%
% Top    : status band per segment. Colors answer, at every instant,
%          "which BRT is the vehicle in?"
%            green  - inside BOTH (handover certified, healthy)
%            blue   - inside CURRENT only (has not reached the next BRT yet)
%            orange - inside NEXT only (already left the current one)
%            red    - inside NEITHER (no liveness certificate at all)
%          The marker at the right edge (100 % = handover instant) is the
%          CBF transition-ready verdict: filled circle = V_next < 0 (ok),
%          cross = V_next >= 0 (NOT ready, governor must intervene).
% Bottom : V_next per segment on the same progress axis, one arc per
%          segment, y clipped so the +-0.1 decisions are visible. The
%          zero crossing of each arc is the moment the vehicle enters the
%          next BRT; an arc rising back above 0 before 100 % is a re-exit.
if nargin < 5, opts = struct(); end
mark = getfield_default(opts, 'mark', []);

nSeg = traj.n_trim - 1;
nB   = 201;
frac = linspace(0, 1, nB);

Vc = nan(nSeg, nB);
Vn = nan(nSeg, nB);
for k = 1:nSeg
    tq = traj.t_node(k) + frac * (traj.t_node(k + 1) - traj.t_node(k));
    for j = 1:nB
        [~, i] = min(abs(guam_trace.time - tq(j)));
        x = guam_trace.lon(:, i);
        Vc(k, j) = brt_value(brt{k},     gv, x - traj.trim_lon(:, k));
        Vn(k, j) = brt_value(brt{k + 1}, gv, x - traj.trim_lon(:, k + 1));
    end
end

% status code: 1 neither, 2 current only, 3 next only, 4 both
S = ones(nSeg, nB);
S(Vc <  0 & Vn >= 0) = 2;
S(Vc >= 0 & Vn <  0) = 3;
S(Vc <  0 & Vn <  0) = 4;

cmapS = [0.85 0.20 0.20;      % 1 neither  - red
         0.30 0.50 0.85;      % 2 cur only - blue
         0.95 0.70 0.20;      % 3 next only- orange
         0.30 0.75 0.40];     % 4 both     - green

fig = figure('Name', 'CBF detection per segment', 'Position', [60 60 1250 900]);

% ------------------------------------------------------------ status bands
ax1 = subplot(2, 1, 1); hold(ax1, 'on');
imagesc(ax1, frac * 100, 1:nSeg, S);
colormap(ax1, cmapS);
clim(ax1, [0.5, 4.5]);
set(ax1, 'YDir', 'reverse');
ylim(ax1, [0.5, nSeg + 0.5]);
xlim(ax1, [0, 100]);

% handover verdict markers just outside the right edge
for k = 1:nSeg
    if Vn(k, end) < 0
        plot(ax1, 103, k, 'o', 'Color', [0.1 0.55 0.2], ...
             'MarkerFaceColor', [0.3 0.75 0.4], 'MarkerSize', 8, ...
             'Clipping', 'off');
    else
        plot(ax1, 103, k, 'x', 'Color', [0.75 0 0], 'LineWidth', 2.5, ...
             'MarkerSize', 11, 'Clipping', 'off');
    end
end
for k = mark(:)'
    plot(ax1, [0 100 100 0 0], k + [-.5 -.5 .5 .5 -.5], 'k-', 'LineWidth', 2);
end

ylabel(ax1, 'segment  (UH_k \rightarrow UH_{k+1})');
yticks(ax1, 1:nSeg);
xlabel(ax1, 'progress within segment [%]        (x = 103: handover verdict)');
title(ax1, 'Where is the vehicle at each instant of each segment?');
% legend proxies
hp = [patch(ax1, 'XData', nan, 'YData', nan, 'FaceColor', cmapS(4, :), 'DisplayName', 'inside BOTH (healthy)'), ...
      patch(ax1, 'XData', nan, 'YData', nan, 'FaceColor', cmapS(2, :), 'DisplayName', 'inside CURRENT only (not reached next yet)'), ...
      patch(ax1, 'XData', nan, 'YData', nan, 'FaceColor', cmapS(3, :), 'DisplayName', 'inside NEXT only (already left current)'), ...
      patch(ax1, 'XData', nan, 'YData', nan, 'FaceColor', cmapS(1, :), 'DisplayName', 'inside NEITHER (no certificate)')];
legend(ax1, hp, 'Location', 'southoutside', 'Orientation', 'horizontal');

% ------------------------------------------------------------ V_next arcs
% Healthy segments are drawn thin/grey; segments that are ever uncertified
% (not ready at handover, or outside both at some point) are colored and
% labelled, so the eye goes straight to the ones that matter.
ax2 = subplot(2, 1, 2); hold(ax2, 'on'); grid(ax2, 'on');
bad = false(1, nSeg);
for k = 1:nSeg
    bad(k) = Vn(k, end) >= 0 || any(Vc(k, :) >= 0 & Vn(k, :) >= 0);
end
cmapL = lines(max(nnz(bad), 1));
ib = 0;
for k = 1:nSeg
    if bad(k)
        ib = ib + 1;
        col = cmapL(ib, :);  lw = 2.2;
        nm  = sprintf('seg %d', k);
        vis = 'on';
    else
        col = [0.6 0.6 0.6];  lw = 0.9;
        nm  = '';  vis = 'off';
    end
    plot(ax2, frac * 100, Vn(k, :), '-', 'Color', col, 'LineWidth', lw, ...
         'DisplayName', nm, 'HandleVisibility', vis);
    if Vn(k, end) < 0
        plot(ax2, 100, Vn(k, end), 'o', 'Color', col, 'MarkerFaceColor', col, ...
             'MarkerSize', 6, 'HandleVisibility', 'off');
    else
        plot(ax2, 100, Vn(k, end), 'x', 'Color', [0.75 0 0], ...
             'LineWidth', 2.5, 'MarkerSize', 12, 'HandleVisibility', 'off');
        text(ax2, 101, Vn(k, end), sprintf(' %d', k), 'Color', [0.75 0 0], ...
             'FontWeight', 'bold', 'FontSize', 9);
    end
end
plot(ax2, nan, nan, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9, ...
     'DisplayName', 'always-certified segments');
yline(ax2, 0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
ylim(ax2, [-0.4, 0.6]);
xlim(ax2, [0, 106]);
xlabel(ax2, 'progress within segment [%]');
ylabel(ax2, 'V_{next}  (V<0 = inside next BRT)');
title(ax2, ['V_{next} at the flown state: crossing 0 downward = enters next BRT, ' ...
            'rising back above 0 = re-exit,  x at 100 % = NOT transition-ready']);
legend(ax2, 'Location', 'eastoutside', 'FontSize', 8);
end

% -------------------------------------------------------------------------
function V = brt_value(data, gv, x_pert)
xc = x_pert(:)';
for d = 1:4
    xc(d) = min(max(xc(d), gv{d}(1)), gv{d}(end));
end
V = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xc(1), xc(2), xc(3), xc(4), 'linear');
end
