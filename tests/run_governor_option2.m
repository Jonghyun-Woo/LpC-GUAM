% run_governor_option2 - closed-loop test of the option-2 trajectory
% time-scaling governor (TrimPointGovernor) on the trim-schedule mission.
%
% Option 2 (Di Cairano, Goldsmith et al. 2015 style): the trim-line path is
% fixed; each period the governor picks how many points to advance
% (0 = hold, 1 = nominal, up to n_max = catch-up), subject to
%   - proximity |u_ref - u_vehicle| <= s_f  (command never runs away)
%   - BRT node gate V_next < -eps_ready     (certified anchor handover)
%   - outside-both hold                     (recursive feasibility)
%
% Prints the same node-crossing report as option 1 and draws the
% plot_governor_run figure. Compare 'final trim reached' against option 1.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 1) Schedule, BRT tables, plant
dt_sim = 0.01;
params = struct('filter_mode', 'off');
params.T_seg = 2.0;      % <<< per-segment time [s]; feeds both the
                         % governor schedule and the plant-side reference
traj = TrimScheduleTrajectory.build( ...
           ReferenceTrajectory.trimOpts(dt_sim, params));
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

%% 2) Governed closed-loop run (option 2)
s_f       = 6.0;    % proximity cap [ft/s] (> natural ~4 ft/s tracking lag)
n_max     = 1.0;    % catch-up DISABLED: n_max >= 1.5 ends in an altitude
                    % divergence the (u,w,q,theta) BRT cannot detect
eps_ready = 0.05;   % node-gate margin

T_max  = 80;
Ng     = round(T_max / dt_sim);
guam   = LpC_GUAM(cfg);
gov    = TrimPointGovernor(traj, brt, gv, dt_sim, s_f, n_max, eps_ready, ...
                           struct('use_predictor', true));
guam.reset();
gov.reset();

G = struct('t', zeros(1, Ng), 'st', zeros(12, Ng), 'kappa', zeros(1, Ng), ...
           'n', zeros(1, Ng), 'progress', zeros(1, Ng), 'V_hat', zeros(1, Ng), ...
           'V_cur', zeros(1, Ng), 'V_next', zeros(1, Ng), 'u_ref', zeros(1, Ng));
for i = 1:Ng
    tt = (i - 1) * dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    G.t(i) = tt;            G.st(:, i) = guam.state;
    G.kappa(i) = gi.kappa;  G.n(i) = gi.n;
    G.progress(i) = gi.progress;
    G.V_cur(i) = gi.V_cur;  G.V_next(i) = gi.V_next;
    G.u_ref(i) = gi.u_ref;  G.V_hat(i) = gi.V_hat;
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 5
        G = truncate(G, i);  break;
    end
end
i_fin = find(G.progress >= traj.n_trim, 1, 'first');
fprintf(['option-2 run: hold %.2f s, catch-up %.2f s, ' ...
         'final trim reached %.2f s (nominal %.2f s, option-1 was 49.65 s)\n'], ...
        gov.n_hold * dt_sim, gov.n_catchup * dt_sim, ...
        G.t(i_fin), traj.t_node(end));

%% 3) Node-crossing report
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
gapG = G.V_cur >= 0 & G.V_next >= 0;
fprintf('\nsamples outside BOTH scheduled BRTs: %d / %d\n', nnz(gapG), numel(G.t));

%% 4) Segment tube + command timeline (edit gov_seg, run this section)
gov_seg = 1;                      % <- segment to inspect (1..19)
plot_governor_run(traj, brt, gv, G, gov_seg);

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
