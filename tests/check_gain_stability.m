% check_gain_stability - how far can the scheduled gains be scaled before the
% held-command closed loop stops being stable?
%
% The ungoverned mission at gain_scale = 1.5 looked fine (BRT violation fell
% from 25.7 % to 7.3 %), but that run sweeps through the speed range and
% never holds a command. The governor's whole question is "what happens if I
% HOLD this command", so what matters is whether each trim point is stable
% when held - i.e. whether Phi is Schur. Scaling an LQR gain is not a safe
% way to raise bandwidth: the design is optimal at scale 1, and going above
% it can and does destabilise.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));

h0 = zeros(ClosedLoopModel.N_XI,1);
h0(1:3)=0.5; h0(4:6)=0.1; h0(7:9)=1e-3; h0(10:12)=1e-3;
h0(13:21)=0.5; h0(22:26)=1e-3; h0(27:32)=1e-2; h0(33:34)=1e-3;

GS_LIST = [1.0 1.1 1.2 1.3 1.4 1.5];
LAM = zeros(numel(GS_LIST), traj.n_trim);

for a = 1:numel(GS_LIST)
    guam = LpC_GUAM(Config('trim_schedule', params));
    guam.controller.gain_scale = GS_LIST(a);
    M = ClosedLoopModel(dt);
    M.build(guam, traj.trim_lon(1,:), h0, 15.0, false);
    LAM(a,:) = cellfun(@(P) max(abs(eig(P))), M.Phi);
    fprintf('gain_scale %.2f | max|eig| = %.6f | unstable trims: %s\n', ...
            GS_LIST(a), max(LAM(a,:)), mat2str(find(LAM(a,:) >= 1)));
end

fprintf('\nmax|eig| per trim (>= 1 means that trim cannot be held)\n');
fprintf('%6s', 'trim');  fprintf('%9.2f', GS_LIST);  fprintf('\n');
fprintf('%s\n', repmat('-', 1, 6 + 9*numel(GS_LIST)));
for k = 1:traj.n_trim
    fprintf('%3d %2.0f', k, traj.trim_lon(1,k));
    for a = 1:numel(GS_LIST)
        if LAM(a,k) >= 1, fprintf('  *%6.4f', LAM(a,k));
        else,             fprintf('   %6.4f', LAM(a,k));  end
    end
    fprintf('\n');
end

ok = GS_LIST(all(LAM < 1, 2));
if isempty(ok)
    fprintf('\nno scale keeps every trim stable\n');
else
    fprintf('\nlargest scale with every trim stable: %.2f\n', max(ok));
end
