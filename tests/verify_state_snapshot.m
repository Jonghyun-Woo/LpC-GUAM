% verify_state_snapshot - is LpC_GUAM.saveState/restoreState complete?
%
% The reference governor predicts by rolling the closed loop forward from a
% snapshot and then rewinding. If the snapshot misses ANY memory (servo
% positions, RSLQR integrators, ...), the prediction silently disagrees with
% what the vehicle would really do.
%
% Test: run N steps from a snapshot, rewind, run the SAME N steps again.
% The two traces must be bit-identical. Then repeat with a deliberately
% incomplete snapshot to show the test actually has teeth.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
cfg  = Config('trim_schedule', params);
guam = LpC_GUAM(cfg);
rt   = guam.refTraj;
guam.reset();

% Fly to a non-trivial point (mid-transient) so integrators are wound up.
i0 = round(16 / dt_sim);
for i = 1:i0
    guam.step(struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                     'chi', rt.chi(i), 'chi_dot', rt.chidot(i)));
end
fprintf('snapshot taken at t = %.2f s, u = %.2f ft/s\n', guam.time, guam.state(4));

N = 300;                       % 3 s of rollout
snap = guam.saveState();

trace1 = run_n(guam, rt, i0, N);
guam.restoreState(snap);
trace2 = run_n(guam, rt, i0, N);

err = max(abs(trace1(:) - trace2(:)));
fprintf('\nfull snapshot   : max |trace1 - trace2| = %.3e', err);
if err == 0
    fprintf('   -> EXACT\n');
else
    fprintf('   -> MISMATCH (snapshot incomplete!)\n');
end

% --- negative control: drop the RSLQR integrators from the snapshot ------
guam.restoreState(snap);
bad = snap;  bad.ctrl_lon_i = zeros(3, 1);  bad.ctrl_lat_i = zeros(3, 1);
guam.restoreState(bad);
trace3 = run_n(guam, rt, i0, N);
err_bad = max(abs(trace1(:) - trace3(:)));
fprintf('without integrators: max |trace1 - trace3| = %.3e', err_bad);
if err_bad > 0
    fprintf('   -> differs, as expected (integrators matter)\n');
else
    fprintf('   -> no difference?! test has no teeth here\n');
end

guam.restoreState(snap);   % leave the object where we found it

% -------------------------------------------------------------------------
function tr = run_n(guam, rt, i0, N)
tr = zeros(4, N);
for j = 1:N
    i = i0 + j;
    guam.step(struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                     'chi', rt.chi(i), 'chi_dot', rt.chidot(i)));
    tr(:, j) = guam.state([4 6 11 8]);
end
end
