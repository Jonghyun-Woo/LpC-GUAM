% run_governor_option1 - closed-loop test of the option-1 scalar reference
% governor (TrimProgressGovernor) on the trim-schedule transition mission,
% against the ungoverned baseline.
%
% Governor rule (see TrimProgressGovernor): the progress along the trim
% schedule advances at the nominal rate, but (a) may cross a trim node only
% when the vehicle is already inside the next BRT, and (b) freezes while
% the vehicle is outside both scheduled BRTs.
%
% Prints per-node handover times/values and draws a comparison figure.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 1) Schedule, BRT tables, plant
dt_sim = 0.01;
params = struct('filter_mode', 'off');
params.T_seg = 2.0;      % <<< per-segment time [s]. Single knob: it feeds
                         % both the governor's schedule and the plant-side
                         % reference, so they cannot drift apart.

traj = TrimScheduleTrajectory.build( ...
           ReferenceTrajectory.trimOpts(dt_sim, params));   % nominal schedule
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

cfg  = Config('trim_schedule', params);

%% 2) Governed closed-loop run
T_max  = 80;                       % allow for governor-induced delay [s]
Ng     = round(T_max / dt_sim);
guam   = LpC_GUAM(cfg);
gov    = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05, ...
                              struct('use_predictor', true));
guam.reset();
gov.reset();

G = struct('t', zeros(1, Ng), 'st', zeros(12, Ng), 'kappa', zeros(1, Ng), ...
           'progress', zeros(1, Ng), 'V_cur', zeros(1, Ng), ...
           'V_next', zeros(1, Ng), 'u_ref', zeros(1, Ng), 'V_hat', zeros(1, Ng));
for i = 1:Ng
    tt = (i - 1) * dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    G.t(i) = tt;            G.st(:, i) = guam.state;
    G.kappa(i) = gi.kappa;  G.progress(i) = gi.progress;
    G.V_cur(i) = gi.V_cur;  G.V_next(i) = gi.V_next;
    G.u_ref(i) = gi.u_ref;  G.V_hat(i) = gi.V_hat;
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 5
        G = truncate(G, i);  break;
    end
end
i_fin = find(G.progress >= traj.n_trim, 1, 'first');
fprintf(['option-1 run: freeze %.2f s, slow %.2f s, ' ...
         'final trim reached %.2f s (nominal %.2f s)\n'], ...
        gov.n_hold * dt_sim, gov.n_slow * dt_sim, ...
        G.t(i_fin), traj.t_node(end));

%% 3) Ungoverned baseline (nominal schedule fed open-loop)
guam2 = LpC_GUAM(cfg);
rt = guam2.refTraj;
M  = size(rt.pos, 2);
guam2.reset();
st0 = zeros(12, M);
for i = 1:M
    ref = struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                 'chi', rt.chi(i), 'chi_dot', rt.chidot(i));
    st0(:, i) = guam2.state;
    guam2.step(ref);
end
base = struct('t', rt.time, 'st', st0);

%% 4) Node-crossing report
% governed crossing time of node k+1 = first t with progress >= k+1
fprintf('\n k+1   nominal t_node   governed cross t   delay [s]   V_next@cross\n');
for k = 1:traj.n_trim - 1
    ic = find(G.progress >= k + 1, 1, 'first');
    if isempty(ic)
        fprintf('%3d   %10.2f        (never)\n', k + 1, traj.t_node(k + 1));
        continue;
    end
    x_lon = G.st([4, 6, 11, 8], ic);
    Vn = brt_value(brt{k + 1}, gv, x_lon - traj.trim_lon(:, k + 1));
    fprintf('%3d   %10.2f        %10.2f      %+7.2f      %+9.3f\n', ...
            k + 1, traj.t_node(k + 1), G.t(ic), ...
            G.t(ic) - traj.t_node(k + 1), Vn);
end

% outside-both statistics (governed anchors follow governed progress)
gapG = G.V_cur >= 0 & G.V_next >= 0;
fprintf('\ngoverned run samples outside BOTH scheduled BRTs: %d / %d\n', ...
        nnz(gapG), numel(G.t));

%% 5) Comparison figure
f = figure('Name', 'option-1 reference governor vs baseline', ...
           'Position', [60 60 1250 850]);
axA = subplot(3, 1, 1); hold(axA, 'on'); grid(axA, 'on');
plot(axA, traj.time, traj.lon(1, :), 'k--', 'LineWidth', 1.0, ...
     'DisplayName', 'nominal schedule u_{ref}');
plot(axA, G.t, G.u_ref, 'b-', 'LineWidth', 1.4, 'DisplayName', 'governed u_{ref}');
plot(axA, base.t, base.st(4, :), '-', 'Color', [1 0 0 0.45], 'LineWidth', 1.2, ...
     'DisplayName', 'flown u (baseline)');
plot(axA, G.t, G.st(4, :), '-', 'Color', [0.55 0 0.8], 'LineWidth', 1.4, ...
     'DisplayName', 'flown u (governed)');
ylabel(axA, 'u [ft/s]'); legend(axA, 'Location', 'southeast');
title(axA, 'Option-1 scalar reference governor on the trim schedule');

axB = subplot(3, 1, 2); hold(axB, 'on'); grid(axB, 'on');
plot(axB, G.t, G.V_cur, 'Color', [0 0.447 0.741], 'LineWidth', 1.3, ...
     'DisplayName', 'V_{cur} (governed run)');
plot(axB, G.t, G.V_next, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.3, ...
     'DisplayName', 'V_{next} (governed run)');
yline(axB, 0, 'k-', 'HandleVisibility', 'off');
ylabel(axB, 'BRT value V'); legend(axB, 'Location', 'northeast');

axC = subplot(3, 1, 3); hold(axC, 'on'); grid(axC, 'on');
stairs(axC, G.t, G.kappa, 'b-', 'LineWidth', 1.0, 'DisplayName', '\kappa');
plot(axC, G.t, (G.progress - 1) / (traj.n_trim - 1), 'Color', [0.55 0 0.8], ...
     'LineWidth', 1.3, 'DisplayName', 'progress (normalized)');
ylim(axC, [-0.05 1.05]);
xlabel(axC, 'time [s]'); legend(axC, 'Location', 'east');
linkaxes([axA axB axC], 'x');

%% 6) Segment tube + command timeline (edit gov_seg, run this section)
% Top: governed window of one segment in (u, t, theta) with BRT tube,
% governed command (black) and flown GUAM (red). Bottom: command vs time,
% governor-active spans shaded green, nominal T_seg grid as dashed lines.
gov_seg = 5;                      % <- segment to inspect (1..19)
plot_governor_run(traj, brt, gv, G, gov_seg);

%% 7) WHY the governor acted in that segment (edit gov_seg, run this section)
% Tube with the flown trajectory RED while the schedule ran and GREEN while
% the command was frozen, the momentum-carrying states, and dV/dt split per
% state so the responsible motion is identifiable.
plot_governor_segment_diag(traj, brt, gv, G, gov_seg);

% -------------------------------------------------------------------------
function G = truncate(G, n)
fn = fieldnames(G);
for ii = 1:numel(fn)
    G.(fn{ii}) = G.(fn{ii})(:, 1:n);
end
end

function V = brt_value(data, gv, x_pert)
xc = x_pert(:)';
for d = 1:4
    xc(d) = min(max(xc(d), gv{d}(1)), gv{d}(end));
end
V = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xc(1), xc(2), xc(3), xc(4), 'linear');
end
