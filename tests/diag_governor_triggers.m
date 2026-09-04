% diag_governor_triggers - which gate fired, when, and why, per segment.
%
% The governed run is a DIFFERENT trajectory from the ungoverned baseline:
% once the schedule is delayed, every later segment executes from a
% different vehicle state. So "segment k was clean in the baseline" does
% not imply it stays clean under the governor. This script reports the
% actual trigger history so that surprises can be traced instead of guessed.
%
% Outputs
%   1) per-segment table: freeze time broken down by trigger reason
%   2) detailed time history for one chosen segment (DETAIL_SEG)
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

DETAIL_SEG = 13;        % <- segment to dump in detail

dt_sim = 0.01;
traj = TrimScheduleTrajectory.build(struct('dt', dt_sim));
spec = FilterConfig.channelSpec('lon');
gv = cell(1, 4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
brt = cell(1, traj.n_trim);
for k = 1:traj.n_trim
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, traj.wh_idx)), 'data');
    brt{k} = S.data;
end
cfg = Config('trim_schedule', struct('filter_mode', 'off'));

%% Governed run with full trigger logging
T_max = 80;  Ng = round(T_max / dt_sim);
guam  = LpC_GUAM(cfg);
gov   = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05);
guam.reset();  gov.reset();

L = struct('t', zeros(1, Ng), 'seg', zeros(1, Ng), 'kappa', zeros(1, Ng), ...
           'V_cur', zeros(1, Ng), 'V_next', zeros(1, Ng), 'V_hat', zeros(1, Ng), ...
           'dVdt', zeros(1, Ng), 'u', zeros(1, Ng), 'th', zeros(1, Ng), ...
           'u_ref', zeros(1, Ng), 'progress', zeros(1, Ng));
gate = cell(1, Ng);
for i = 1:Ng
    tt = (i - 1) * dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    L.t(i) = tt;              L.seg(i) = gi.seg;
    L.kappa(i) = gi.kappa;    L.V_cur(i) = gi.V_cur;
    L.V_next(i) = gi.V_next;  L.V_hat(i) = gi.V_hat;
    L.dVdt(i) = gi.dVdt;      L.u(i) = guam.state(4);
    L.th(i) = guam.state(8);  L.u_ref(i) = gi.u_ref;
    L.progress(i) = gi.progress;
    gate{i} = gi.gate;
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 2
        fn = fieldnames(L);
        for ii = 1:numel(fn), L.(fn{ii}) = L.(fn{ii})(1:i); end
        gate = gate(1:i);
        break;
    end
end

%% 1) Per-segment trigger breakdown
reasons = {'outside-both-hold', 'predict-reexit-freeze', 'node-gate-not-ready'};
fprintf('\n k  |   window [s]    | freeze | outside-both | reexit-pred | node-gate\n');
fprintf('%s\n', repmat('-', 1, 78));
for k = 1:traj.n_trim - 1
    m = find(L.seg == k);
    if isempty(m), continue; end
    tf = zeros(1, numel(reasons));
    for r = 1:numel(reasons)
        tf(r) = nnz(strcmp(gate(m), reasons{r})) * dt_sim;
    end
    tot = nnz(L.kappa(m) == 0) * dt_sim;
    flag = '';
    if tot > 0.05, flag = '  <-- governor acted'; end
    fprintf('%2d  | %6.2f - %6.2f |  %5.2f |    %5.2f     |    %5.2f    |   %5.2f%s\n', ...
            k, L.t(m(1)), L.t(m(end)), tot, tf(1), tf(2), tf(3), flag);
end

%% 2) Detail dump for one segment
k = DETAIL_SEG;
m = find(L.seg == k);
fprintf('\n=== segment %d detail (UH%d -> UH%d) ===\n', k, k, k + 1);
fprintf('   t     u      th    | V_cur   V_next   V_hat   dVdt  | kappa  gate\n');
step = max(1, round(numel(m) / 25));
for j = 1:step:numel(m)
    i = m(j);
    fprintf('%6.2f %6.1f %+6.2f | %+6.3f %+7.3f %+7.3f %+6.2f |  %3.1f   %s\n', ...
        L.t(i), L.u(i), rad2deg(L.th(i)), L.V_cur(i), L.V_next(i), ...
        L.V_hat(i), L.dVdt(i), L.kappa(i), gate{i});
end

%% 3) Figure: the chosen segment in time
f = figure('Name', sprintf('segment %d trigger detail', k), ...
           'Position', [60 60 1200 800]);
tw = [L.t(m(1)), L.t(m(end))];
axd(1) = subplot(3, 1, 1); hold on; grid on;
plot(L.t(m), L.u(m), 'r-', 'LineWidth', 1.6, 'DisplayName', 'flown u');
plot(L.t(m), L.u_ref(m), 'b-', 'LineWidth', 1.4, 'DisplayName', 'commanded u_{ref}');
ylabel('u [ft/s]'); legend('Location', 'best');
title(sprintf('Segment %d: what the governor saw', k));

axd(2) = subplot(3, 1, 2); hold on; grid on;
plot(L.t(m), rad2deg(L.th(m)), 'Color', [0.8 0.4 0], 'LineWidth', 1.6);
ylabel('\theta [deg]');

axd(3) = subplot(3, 1, 3); hold on; grid on;
plot(L.t(m), L.V_cur(m),  'Color', [0 0.447 0.741], 'LineWidth', 1.5, ...
     'DisplayName', 'V_{cur}');
plot(L.t(m), L.V_next(m), 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5, ...
     'DisplayName', 'V_{next}');
plot(L.t(m), L.V_hat(m),  'Color', [0.2 0.2 0.85], 'LineWidth', 1.2, ...
     'DisplayName', 'V_{hat} (forecast)');
yline(0, 'k-', 'HandleVisibility', 'off');
frz = m(L.kappa(m) == 0);
if ~isempty(frz)
    plot(L.t(frz), zeros(1, numel(frz)), 'g.', 'MarkerSize', 10, ...
         'DisplayName', 'governor frozen');
end
ylim([-0.5 0.5]); ylabel('BRT value'); xlabel('time [s]');
legend('Location', 'best');
linkaxes(axd, 'x'); xlim(axd(1), tw);
