% verify_two_line_change - confirm the two edited lines reach the running
% controller through the ordinary construction path.
%
% Everything so far was measured with the gains and the mode forced onto a
% local copy from inside a test script. This builds the vehicle the normal
% way and touches nothing, so it also checks that update_gains = true really
% picks up the new Qlon0.
%
% Expected, from the earlier sweeps:
%   Ki(1,1) = sqrt(q1/r1) = sqrt(0.08) = 0.2828
%   w_ff_mode = 2
%   peak speed error 4.45 ft/s, peak position error 9.6 ft, 2.24 % outside
clear;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

guam = LpC_GUAM(Config('trim_schedule', params));   % nothing forced
bc   = guam.controller.baseline_controller;
cfg  = RSLQRConfig;

fprintf('\nQlon0        = [%s]\n', strtrim(sprintf('%g ', cfg.Qlon0)));
fprintf('update_gains = %d\n', cfg.update_gains);
fprintf('w_ff_mode    = %d\n', bc.w_ff_mode);
fprintf('k_pos        = %g\n\n', bc.k_pos);

fprintf('%-6s %10s %10s %10s\n', 'trim', 'Ki(1,1)', 'sqrt(q1)', 'Kx(1,1)');
fprintf('%s\n', repmat('-',1,40));
for i = [1 6 12 20]
    fprintf('%-6d %10.4f %10.4f %10.4f\n', i, ...
            bc.LON.Ki(1,1,i,1), sqrt(cfg.Qlon0(1)/cfg.Rlon0(1)), bc.LON.Kx(1,1,i,1));
end

guam.reset();
rt = guam.refTraj;  n = size(rt.pos,2);
u=zeros(1,n); z=zeros(1,n); ep=zeros(1,n); V=zeros(1,n);
for i = 1:n
    s = guam.state;
    R = RSLQR.rotm_i2b(s(7),s(8),s(9));
    e = R*(rt.pos(:,i) - s(1:3));
    u(i)=s(4); z(i)=s(3); ep(i)=e(1);
    V(i) = brtV.value(s([4 6 11 8]), rt.vel(1,i));
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
k = rt.time <= traj.t_node(end) + 3;
alt = -z(k);  t = rt.time(k);  late = t >= 15;

fprintf('\nmission, built the ordinary way:\n');
fprintf('  peak speed error      %8.2f ft/s   (expected 4.45)\n', ...
        max(abs(u(k) - rt.vel(1,k))));
fprintf('  peak position error   %8.1f ft     (expected  9.6)\n', max(abs(ep(k))));
fprintf('  peak altitude dev     %8.1f ft\n', max(abs(alt - 80)));
fprintf('  altitude swing        %8.1f ft\n', max(alt(late)) - min(alt(late)));
fprintf('  outside the tube      %8.2f %%     (expected 2.24)\n', ...
        100*nnz(V(k) >= 0)/nnz(k));
