% run_governor_rg_linear - the same governed mission as run_governor_rg, but
% predicting with the linearised closed loop instead of stepping the plant.
%
% Everything else is identical: same schedule, same BRT, same Algorithm 2,
% same delta/eps/M. Only TrimRefGovernor.predict changes, so any difference
% in the result is attributable to the prediction model alone.
%
% Reference numbers from the nonlinear rollout (T_seg = 2.0, M = 40):
%   final command 51.70 s | V >= 0 in 0.0 % of the mission | worst V -0.000
%   altitude deviation 29.2 ft | wall clock 1285 s
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
GOV_HZ = 10;  N_dec = round(1/(GOV_HZ*dt_sim));
T_HOR  = 6.0;                 % keep the same horizon for a like-for-like run

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));

S = load(fullfile(root, 'logger', 'closed_loop_model.mat'), 'M');
lin = S.M;
fprintf('loaded closed-loop model: %d trim points, max|eig| = %.6f\n', ...
        numel(lin.Phi), max(cellfun(@(P) max(abs(eig(P))), lin.Phi)));

gov = TrimRefGovernor(traj, brtV, guam, dt_sim, ...
        struct('T', T_HOR, 'delta', 0.3, 'eps', 0.02, 'M', 40, 'lin', lin));
guam.reset();  gov.reset();

%% Run
SIM_MAX = 90;  Ng = round(SIM_MAX/dt_sim);
L = struct('t', zeros(1,Ng), 'st', zeros(12,Ng), 'v', zeros(1,Ng), ...
           'r', zeros(1,Ng), 'kappa', zeros(1,Ng), 'Vnow', zeros(1,Ng));
gate = cell(1,Ng);  last = [];
v_final = traj.trim_lon(1,end);  t_done = NaN;
tic;
for i = 1:Ng
    tt = (i-1)*dt_sim;
    if mod(i-1, N_dec) == 0, [ref, gi] = gov.step(tt);  last = gi;
    else,                    ref = gov.hold();  end
    L.t(i) = tt;                L.st(:,i) = guam.state;
    L.v(i) = gov.v;             L.r(i) = last.r;
    L.kappa(i) = last.kappa;    gate{i} = last.gate;
    L.Vnow(i) = brtV.value(guam.state([4 6 11 8]), gov.v);
    if isnan(t_done) && gov.v >= v_final - 1e-6, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done + 3, break; end
end
n = i;  fn = fieldnames(L);
for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:,1:n); end
gate = gate(1:n);  wall = toc;

%% Report, alongside the nonlinear-rollout numbers
fprintf('\n=== linear-prediction governor, T_seg = %.2f s, T = %.1f s ===\n', ...
        params.T_seg, gov.T);
fprintf('%-24s %12s %12s\n', '', 'LINEAR', 'nonlinear');
fprintf('%s\n', repmat('-', 1, 50));
fprintf('%-24s %12.2f %12.2f\n', 'final command [s]', t_done, 51.70);
fprintf('%-24s %12.1f %12.1f\n', 'V >= 0 [% of mission]', ...
        100*nnz(L.Vnow >= 0)/n, 0.0);
fprintf('%-24s %+12.3f %+12.3f\n', 'worst V', max(L.Vnow), -0.000);
fprintf('%-24s %12.1f %12.1f\n', 'altitude deviation [ft]', ...
        max(abs(-L.st(3,:) - 80)), 29.2);
fprintf('%-24s %12.1f %12.1f\n', 'wall clock [s]', wall, 1285);
fprintf('%-24s %12.0f %12.0f\n', 'speedup', 1285/max(wall,eps), 1);
fprintf('%-24s %12d %12d\n', 'rollouts', gov.n_rollout, 18356);

frac = @(s) 100*nnz(strcmp(gate, s))/n;
fprintf('\ngate breakdown [%% of time]\n');
for g = unique(gate), fprintf('   %-18s %5.1f\n', g{1}, frac(g{1})); end
fprintf('\nschedule not followed : %.1f %%\n', 100*nnz(L.kappa < 1-1e-9)/n);

save(fullfile(root,'logger','rg_linear_run.mat'), 'L', 'gate', 'wall', '-v7.3');
