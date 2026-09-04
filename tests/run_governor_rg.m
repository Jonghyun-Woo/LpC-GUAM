% run_governor_rg - main entry point for the reference-governor transition.
%
% The governor sits between the mission plan and the autopilot. Before passing on
% the speed the plan asks for, it simulates holding that speed for T seconds and
% keeps it only if the liveness certificate survives; otherwise it backs off along
% v = v(t-1) + kappa*(r - v(t-1)) to the largest surviving kappa, or holds v(t-1).
%
% Implementation notes: V comes from rolling the real closed loop forward
% (G&K 2002, Eq. 23) rather than reading the tube directly; V is continuous in
% the command (BRTValue); kappa is searched on a grid because V is not monotone
% in kappa (tests/check_kappa_monotone.m).
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% Setup
dt_sim = 0.01;
params = struct('filter_mode', 'off');
params.T_seg = 2.0;              % how fast the PLAN wants to move [s/segment]

GOV_HZ = 10;                     % governor solve rate
N_dec  = round(1/(GOV_HZ*dt_sim));

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
cfg  = Config('trim_schedule', params);
guam = LpC_GUAM(cfg);

gov = TrimRefGovernor(traj, brtV, guam, dt_sim, ...
        struct('T', 6.0, 'delta', 0.3, 'eps', 0.02, 'M', 40));

guam.reset();  gov.reset();

%% Governed run
SIM_MAX = 90;  Ng = round(SIM_MAX/dt_sim);
L = struct('t', zeros(1,Ng), 'st', zeros(12,Ng), 'v', zeros(1,Ng), ...
           'r', zeros(1,Ng), 'kappa', zeros(1,Ng), 'Vnow', zeros(1,Ng));
gate = cell(1,Ng);
dec  = struct('t', [], 'Vs', [], 'Vh', [], 'k', [], 'cands', {{}}, 'reject', []);
v_final = traj.trim_lon(1, end);
t_done  = NaN;
tic;
for i = 1:Ng
    tt = (i-1)*dt_sim;
    if mod(i-1, N_dec) == 0
        [ref, gi] = gov.step(tt);
        last_gi = gi;
        dec.t(end+1)      = tt;
        dec.Vs(end+1)     = gi.V_sched;
        dec.Vh(end+1)     = gi.V_hold;
        dec.k(end+1)      = gi.kappa;
        dec.cands{end+1}  = gi.cands;
        dec.reject(end+1) = ~isnan(gi.V_sched) && gi.V_sched > 0;
    else
        ref = gov.hold();
    end
    L.t(i) = tt;                L.st(:,i) = guam.state;
    L.v(i) = gov.v;             L.r(i) = last_gi.r;
    L.kappa(i) = last_gi.kappa; gate{i} = last_gi.gate;
    L.Vnow(i) = brtV.value(guam.state([4 6 11 8]), gov.v);

    if isnan(t_done) && gov.v >= v_final - 1e-6, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done + 3, break; end
end
n = i;  fn = fieldnames(L);
for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:,1:n); end
gate = gate(1:n);
wall = toc;

L.dec_t = dec.t;  L.dec_Vs = dec.Vs;  L.dec_Vh = dec.Vh;
L.dec_k = dec.k;  L.dec_cands = dec.cands;
L.dec_reject = logical(dec.reject);
L.t_done = t_done;

%% Report
fprintf('\n=== reference governor, T_seg = %.2f s, T = %.1f s ===\n', ...
        params.T_seg, gov.T);
fprintf('final command reached : %.2f s   (plan wanted %.2f s, delay %+.2f s)\n', ...
        t_done, traj.t_node(end), t_done - traj.t_node(end));
fprintf('rollouts              : %d   (%.1f per governor solve)\n', ...
        gov.n_rollout, gov.n_rollout / max(gov.n_call,1));
fprintf('wall clock            : %.1f s\n', wall);

frac = @(s) 100*nnz(strcmp(gate, s))/n;
fprintf('\ngate breakdown [%% of time]\n');
for g = unique(gate), fprintf('   %-18s %5.1f\n', g{1}, frac(g{1})); end

fprintf('\nschedule not followed  : %.1f %% of the mission\n', ...
        100*nnz(L.kappa < 1-1e-9)/n);
fprintf('r judged inadmissible  : %d of %d governor solves\n', ...
        nnz(L.dec_reject), numel(L.dec_t));
fprintf('samples with V >= 0    : %d / %d (%.1f %%)\n', ...
        nnz(L.Vnow >= 0), n, 100*nnz(L.Vnow >= 0)/n);
fprintf('worst V                : %+.3f\n', max(L.Vnow));
fprintf('theta range            : %+.2f .. %+.2f deg\n', ...
        rad2deg(min(L.st(8,:))), rad2deg(max(L.st(8,:))));
fprintf('altitude deviation     : %.1f ft\n', max(abs(-L.st(3,:) - 80)));

%% Figures
% Figures 3 and 4 zoom into one window. Set seg to the transition segment you
% want to look at: seg = k covers trim point k -> k+1, so with 20 trim points
% seg runs 1..19. Leave seg empty to fall back on win_t/win_half (or, with
% win_t empty too, on whichever solve was most constrained).
%
% Re-run THIS SECTION ONLY (Ctrl+Enter) to change the window - L is already in
% the workspace, so the simulation does not have to run again.
popts = struct('seg', 6, 'seg_pad', 0.5, 'win_t', [], 'win_half', 2.0);

figs = plot_rg_run(traj, brtV, L, popts);   %#ok<NASGU>