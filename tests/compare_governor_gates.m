% compare_governor_gates - which predictive mode actually helps?
%
% The predictive gate has two failure modes, and they are NOT symmetric:
%
%   LATE ENTRY : acting on it means slowing the command down. But the
%                vehicle only enters BRT_{k+1} because the command is
%                marching toward trim k+1 (the center of that BRT), so
%                slowing the command also slows the entry. The response
%                may be self-defeating.
%   RE-EXIT    : acting on it means freezing before the excursion develops,
%                which genuinely prevents the departure.
%
% This script runs the option-1 governor with each combination and reports
% mission time, certification margins and time spent uncertified, so the
% choice is made on evidence rather than on the elegance of the formula.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

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

CASES = { ...
    'present-tense only', struct('use_predictor', false); ...
    'predict re-exit'   , struct('use_predictor', true,  'predict_entry', false, 'predict_reexit', true); ...
    'predict entry'     , struct('use_predictor', true,  'predict_entry', true,  'predict_reexit', false); ...
    'predict both'      , struct('use_predictor', true,  'predict_entry', true,  'predict_reexit', true)};

fprintf('\n%-20s | mission | freeze | slow  | worst   | uncert.\n', 'gate configuration');
fprintf('%-20s |   [s]   |  [s]   | [s]   | margin  |  [%%]\n', '');
fprintf('%s\n', repmat('-', 1, 74));

R = struct();
for c = 1:size(CASES, 1)
    [G, gov] = run_case(traj, brt, gv, cfg, dt_sim, CASES{c, 2});
    i_fin = find(G.progress >= traj.n_trim, 1, 'first');
    if isempty(i_fin), t_fin = NaN; else, t_fin = G.t(i_fin); end

    % worst certification margin over the 19 node crossings
    worst = -Inf;
    for k = 1:traj.n_trim - 1
        ic = find(G.progress >= k + 1, 1, 'first');
        if isempty(ic), worst = NaN; break; end
        Vn = brt_value(brt{k + 1}, gv, ...
                       G.st([4 6 11 8], ic) - traj.trim_lon(:, k + 1));
        worst = max(worst, Vn);
    end
    uncert = 100 * nnz(G.V_cur >= 0 & G.V_next >= 0) / numel(G.t);

    fprintf('%-20s | %6.2f  | %5.2f  | %5.2f | %+7.3f | %5.1f\n', ...
            CASES{c, 1}, t_fin, gov.n_hold * dt_sim, gov.n_slow * dt_sim, ...
            worst, uncert);
    R(c).name = CASES{c, 1};  R(c).G = G;  R(c).t_fin = t_fin;
end

fprintf(['\nworst margin = the least-negative V_next at any of the 19 node\n' ...
         'crossings (must stay below -eps_ready = -0.05); uncert. = %% of\n' ...
         'samples outside BOTH scheduled BRTs.\n']);

%% Command overlay of all four gate configurations
f = figure('Name', 'gate configuration comparison', 'Position', [60 60 1250 620]);
ax = axes(f); hold(ax, 'on'); grid(ax, 'on');
plot(ax, traj.time, traj.lon(1, :), 'k--', 'LineWidth', 1.0, ...
     'DisplayName', 'nominal schedule');
col = lines(numel(R));
for c = 1:numel(R)
    plot(ax, R(c).G.t, R(c).G.u_ref, '-', 'Color', col(c, :), 'LineWidth', 1.6, ...
         'DisplayName', sprintf('%s (%.1f s)', R(c).name, R(c).t_fin));
end
xlabel(ax, 'time [s]'); ylabel(ax, 'governed u_{ref} [ft/s]');
title(ax, 'Option-1 governor: effect of each predictive mode on the command');
legend(ax, 'Location', 'southeast');

% -------------------------------------------------------------------------
function [G, gov] = run_case(traj, brt, gv, cfg, dt_sim, opts)
T_max = 90;
Ng    = round(T_max / dt_sim);
guam  = LpC_GUAM(cfg);
gov   = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05, opts);
guam.reset();  gov.reset();
G = struct('t', zeros(1, Ng), 'st', zeros(12, Ng), 'kappa', zeros(1, Ng), ...
           'progress', zeros(1, Ng), 'V_cur', zeros(1, Ng), ...
           'V_next', zeros(1, Ng), 'u_ref', zeros(1, Ng));
for i = 1:Ng
    tt = (i - 1) * dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    G.t(i) = tt;            G.st(:, i) = guam.state;
    G.kappa(i) = gi.kappa;  G.progress(i) = gi.progress;
    G.V_cur(i) = gi.V_cur;  G.V_next(i) = gi.V_next;
    G.u_ref(i) = gi.u_ref;
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 3
        fn = fieldnames(G);
        for ii = 1:numel(fn), G.(fn{ii}) = G.(fn{ii})(:, 1:i); end
        break;
    end
end
end

function V = brt_value(data, gv, x_pert)
xc = x_pert(:)';
for d = 1:4
    xc(d) = min(max(xc(d), gv{d}(1)), gv{d}(end));
end
V = interpn(gv{1}, gv{2}, gv{3}, gv{4}, data, xc(1), xc(2), xc(3), xc(4), 'linear');
end
