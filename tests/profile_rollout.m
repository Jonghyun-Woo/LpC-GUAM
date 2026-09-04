% profile_rollout - where does a governor rollout actually spend its time?
%
% One TrimRefGovernor.predict() call rolls the closed loop forward T seconds.
% At 33.5 rollouts per solve and 0.107 s per rollout the governor runs ~38x
% slower than real time, so before porting anything to C it is worth knowing
% which part of the loop that time belongs to: the aero polynomial, the
% controller, the rigid-body integration, or MATLAB overhead.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
params = struct('filter_mode', 'off');
params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));
gov  = TrimRefGovernor(traj, brtV, guam, dt_sim, ...
        struct('T', 6.0, 'delta', 0.3, 'eps', 0.02, 'M', 40));
guam.reset();  gov.reset();

%% Fly to a mid-mission state so the rollout is representative
N_dec = 10;
for i = 1:round(12/dt_sim)
    tt = (i-1)*dt_sim;
    if mod(i-1, N_dec) == 0, ref = gov.step(tt); else, ref = gov.hold(); end
    guam.step(ref);
end
fprintf('warm-up done: t = 12 s, v = %.1f ft/s\n', gov.v);

%% Time a rollout without the profiler (the profiler itself is not free)
n_rep = 10;
v_test = gov.v + 3;
tic;
for i = 1:n_rep, gov.predict(v_test, false); end       % full horizon, no early exit
t_raw = toc / n_rep;
N = gov.N;
fprintf('\nclean timing (no profiler)\n');
fprintf('  full-horizon rollout : %.4f s   (%d steps of %.3f ms)\n', ...
        t_raw, N, 1e3*t_raw/N);
fprintf('  10 Hz budget         : %.4f s   -> %.1fx too slow per rollout\n', ...
        1/10, t_raw/(1/10));

%% Profile
profile clear; profile on;
for i = 1:n_rep, gov.predict(v_test, false); end
profile off;
P = profile('info');

[~, ord] = sort([P.FunctionTable.TotalTime], 'descend');
fprintf('\ntop functions by total time (%d rollouts, %d steps each)\n', n_rep, N);
fprintf('%-42s %8s %8s %10s\n', 'function', 'total[s]', 'self[s]', 'calls');
fprintf('%s\n', repmat('-', 1, 72));
tot = 0;
for a = 1:min(18, numel(ord))
    f = P.FunctionTable(ord(a));
    nm = f.FunctionName;
    if numel(nm) > 42, nm = ['...' nm(end-38:end)]; end
    fprintf('%-42s %8.3f %8.3f %10d\n', ...
            nm, f.TotalTime, f.TotalTime - sum([f.Children.TotalTime]), f.NumCalls);
    tot = tot + f.TotalTime;
end

%% Group the self time into the parts that matter for a C port
grp = struct('aero', 0, 'control', 0, 'rbd', 0, 'brt', 0, 'other', 0);
for a = 1:numel(P.FunctionTable)
    f = P.FunctionTable(a);
    self = f.TotalTime - sum([f.Children.TotalTime]);
    nm = lower(f.FunctionName);
    if contains(nm, 'aero') || contains(nm, 'prop') || contains(nm, 'rotor') ...
            || contains(nm, 'wing') || contains(nm, 'fuse') || contains(nm, 'tail')
        grp.aero = grp.aero + self;
    elseif contains(nm, 'rslqr') || contains(nm, 'alloc') || contains(nm, 'ctrl')
        grp.control = grp.control + self;
    elseif contains(nm, 'rigidbody') || contains(nm, 'dynamics') || contains(nm, 'rotm')
        grp.rbd = grp.rbd + self;
    elseif contains(nm, 'brtvalue') || contains(nm, 'interpn')
        grp.brt = grp.brt + self;
    else
        grp.other = grp.other + self;
    end
end
tt = grp.aero + grp.control + grp.rbd + grp.brt + grp.other;
fprintf('\nself time by group\n');
fn = fieldnames(grp);
for a = 1:numel(fn)
    fprintf('  %-9s %7.3f s   %5.1f %%\n', fn{a}, grp.(fn{a}), 100*grp.(fn{a})/tt);
end
fprintf('  %-9s %7.3f s\n', 'TOTAL', tt);
