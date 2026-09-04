% diag_wff - does supplying the level-flight vertical speed remove the
% standing altitude offset?
%
% Two tests per mode:
%   (a) free response - hold a fixed speed and see where the altitude settles.
%       This is where the offset was measured (151 ft at 100 ft/s) and where
%       the feedforward should show up most clearly.
%   (b) the mission - the number that actually matters.
%
% mode 0 is the shipped behaviour and is the default in RSLQR, so nothing
% changes unless a caller sets it.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
MODE = [0 1 2];
NM   = {'0  off (shipped)', '1  ff from command', '2  ff from flown theta'};
V_TEST = 100;

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
g0   = LpC_GUAM(Config('trim_schedule', params));
bc0  = g0.controller.baseline_controller;
cfg  = RSLQRConfig;

Q = [0.08 0.01 1000 2.0 0 0]';
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

fprintf('\nq1 = 0.08, q4 = 2.0, k_pos = 0.1\n');
fprintf('\n(a) hold u = %d ft/s for 40 s, where does the altitude settle?\n\n', V_TEST);
fprintf('%-24s %14s\n', 'w feedforward', 'altitude [ft]');
fprintf('%s\n', repmat('-', 1, 40));
for a = 1:numel(MODE)
    z = settle_alt(params, Ki, Kx, MODE(a), V_TEST, dt);
    fprintf('%-24s %14.1f\n', NM{a}, z);
end
fprintf('%-24s %14.1f\n', 'reference', 80.0);

fprintf('\n(b) full mission\n\n');
fprintf('%-24s %9s %9s %9s %10s %8s\n', 'w feedforward', ...
        'max|u-r|','e_pos','alt dev','alt swing','BRT>=0');
fprintf('%s\n', repmat('-', 1, 74));
S = cell(1,numel(MODE));
for a = 1:numel(MODE)
    S{a} = mission(params, traj, brtV, Ki, Kx, MODE(a));
    if S{a}.blown
        fprintf('%-24s %s\n', NM{a}, 'DIVERGED');
    else
        fprintf('%-24s %9.2f %9.1f %9.1f %10.1f %8.1f\n', NM{a}, ...
                S{a}.eu, S{a}.ep, S{a}.ad, S{a}.swing, S{a}.viol);
    end
end

%% figure
f = figure('Name','level-flight w feedforward','Position',[60 60 1150 780]);
ax = gobjects(1,3);
co = {[.85 .2 .2], [0 .5 .2], [0 .35 .75]};
ax(1) = subplot(3,1,1); hold on; grid on;
ax(2) = subplot(3,1,2); hold on; grid on;
ax(3) = subplot(3,1,3); hold on; grid on;
for a = 1:numel(MODE)
    if S{a}.blown, continue; end
    plot(ax(1), S{a}.t, S{a}.u - S{a}.r, '-','Color',co{a},'LineWidth',1.4, ...
         'DisplayName', NM{a});
    plot(ax(2), S{a}.t, S{a}.alt, '-','Color',co{a},'LineWidth',1.5);
    plot(ax(3), S{a}.t, S{a}.w, '-','Color',co{a},'LineWidth',1.4);
end
yline(ax(1),0,'k-','HandleVisibility','off');
yline(ax(2),80,'k--');  yline(ax(3),0,'k-');
ylabel(ax(1),'u error [ft/s]');  legend(ax(1),'Location','southwest');
ylabel(ax(2),'altitude [ft]');
ylabel(ax(3),'body w [ft/s]');  xlabel(ax(3),'time [s]');
title(ax(1),'level-flight vertical-speed feedforward');
linkaxes(ax,'x');  xlim(ax(1),[0 S{1}.t(end)]);
exportgraphics(f, fullfile(root,'logger','wff.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','wff.fig'));
fprintf('\nsaved logger/wff(.fig/.png)\n');

% =========================================================================
function z = settle_alt(params, Ki, Kx, mode, v, dt)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.w_ff_mode = mode;
rt = guam.refTraj;  m = size(rt.pos,2);
rt.time=(0:m-1)*dt; rt.vel=repmat([v;0;0],1,m);
rt.pos=[(0:m-1)*dt*v; zeros(1,m); -80*ones(1,m)];
rt.chi=zeros(1,m); rt.chidot=zeros(1,m);
guam.refTraj = rt;  guam.reset();
p = 0;
for i = 1:round(40/dt)
    p = p + v*dt;
    guam.step(struct('pos',[p;0;-80],'vel',[v;0;0],'chi',0,'chi_dot',0));
end
z = -guam.state(3);
end

function S = mission(params, traj, brtV, Ki, Kx, mode)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.w_ff_mode = mode;
guam.reset();
rt = guam.refTraj;  n = size(rt.pos,2);
t = rt.time;  u=zeros(1,n); z=zeros(1,n); w=zeros(1,n); ep=zeros(1,n); nv=0;
S.blown = false;
for i = 1:n
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, S.blown = true; return; end
    R = RSLQR.rotm_i2b(s(7),s(8),s(9));
    e = R*(rt.pos(:,i) - s(1:3));
    u(i)=s(4); z(i)=s(3); w(i)=s(6); ep(i)=e(1);
    nv = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
keep = t <= traj.t_node(end) + 3;
S.t=t(keep); S.u=u(keep); S.w=w(keep); S.r=rt.vel(1,keep);
S.alt=-z(keep);  late = S.t >= 15;
S.eu=max(abs(S.u-S.r)); S.ep=max(abs(ep(keep)));
S.ad=max(abs(S.alt-80)); S.swing=max(S.alt(late))-min(S.alt(late));
S.viol=100*nv/n;
end
