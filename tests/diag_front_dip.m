% diag_front_dip - what makes theta dive ~12 deg below the plan right after
% the speed command starts?
%
% Three questions, in order:
%
%   (a) is the dive COMMANDED or is it a tracking failure?
%       theta enters the allocation as a virtual effector, so there is an
%       explicit theta command (act_lon(end), stored as theta_v). Log it
%       against the flown theta. If they agree, the controller asked for it.
%
%   (b) does it scale with the commanded acceleration?
%       The reference slope steps from 0 to 1/T_seg at t = t_hold. Stretch
%       T_seg; if the dive shrinks in proportion, the step is the cause.
%
%   (c) does it shrink when theta is made expensive to use?
%       Wlon0(12) = 0.1 makes theta the cheapest effector by an order of
%       magnitude. Raise it and see.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q  = [0.08 0.01 1000 2.0 0 0]';
base = struct('filter_mode','off');  base.T_seg = 2.0;

g0  = LpC_GUAM(Config('trim_schedule', base));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

%% (a) is the dive commanded?
R = run_one(base, Ki, Kx, 1.0, dt);
fprintf('\n(a) is the dive commanded, or a tracking failure?\n\n');
[~, j] = min(R.dth);
fprintf('  worst point of the front dip, t = %.2f s\n', R.t(j));
fprintf('    plan theta                     %8.2f deg\n', rad2deg(R.thp(j)));
fprintf('    theta the allocator asked for  %8.2f deg\n', rad2deg(R.thp(j)+R.thv(j)));
fprintf('    theta actually flown           %8.2f deg\n', rad2deg(R.th(j)));
fprintf('    -> commanded away from plan    %8.2f deg\n', rad2deg(R.thv(j)));
fprintf('    -> failed to track its own cmd %8.2f deg\n', ...
        rad2deg(R.th(j)-R.thp(j)-R.thv(j)));
fprintf('\n  reference speed slope, ft/s^2:  before t=5  %.3f   after  %.3f\n', ...
        R.dr(round(4.0/dt)), R.dr(round(5.5/dt)));

%% (b) stretch the ramp
fprintf('\n\n(b) stretch the segment time (gentler commanded acceleration)\n\n');
TS = [2.0 3.0 4.0 6.0];
fprintf('%8s %10s | %11s %9s %9s\n', ...
        'T_seg','du/dt','front dip','max V','out [s]');
fprintf('%s\n', repmat('-',1,54));
for a = 1:numel(TS)
    p = base;  p.T_seg = TS(a);
    r = run_one(p, Ki, Kx, 1.0, dt);
    fprintf('%8.1f %10.3f | %8.2f deg %9.3f %9.2f\n', ...
            TS(a), r.dr(round(5.5/dt)), rad2deg(min(r.dth)), r.Vmax, r.tout);
end

%% (c) make theta expensive
fprintf('\n\n(c) raise the allocation weight on theta (Wlon0(12) = 0.1 shipped)\n\n');
WT = [0.1 0.5 1.0 5.0 20.0];
fprintf('%8s | %11s %9s %9s | %9s %9s\n', ...
        'w_theta','front dip','max V','out [s]','max|u-r|','e_pos');
fprintf('%s\n', repmat('-',1,66));
S = cell(1,numel(WT));
for a = 1:numel(WT)
    S{a} = run_one(base, Ki, Kx, WT(a)/0.1, dt);
    r = S{a};
    if r.blown
        fprintf('%8.2f | %s\n', WT(a), 'DIVERGED');
    else
        fprintf('%8.2f | %8.2f deg %9.3f %9.2f | %9.2f %9.1f\n', ...
                WT(a), rad2deg(min(r.dth)), r.Vmax, r.tout, r.eu, r.ep);
    end
end

%% figure
f = figure('Name','front theta dip','Position',[60 60 1150 820]);
ax = gobjects(1,4);
for a = 1:4, ax(a) = subplot(4,1,a); hold on; grid on; end
co = lines(numel(WT));
for a = 1:numel(WT)
    r = S{a};  if r.blown, continue; end
    nm = sprintf('w_\\theta = %.1f', WT(a));
    if abs(WT(a)-0.1) < 1e-9, nm = [nm '  (shipped)']; end %#ok<AGROW>
    plot(ax(1), r.t, rad2deg(r.th), '-','Color',co(a,:),'LineWidth',1.4, ...
         'DisplayName', nm);
    plot(ax(2), r.t, rad2deg(r.thv),'-','Color',co(a,:),'LineWidth',1.4);
    plot(ax(3), r.t, r.pusher,      '-','Color',co(a,:),'LineWidth',1.4);
    plot(ax(4), r.t, r.V,           '-','Color',co(a,:),'LineWidth',1.4);
end
plot(ax(1), R.t, rad2deg(R.thp), 'k--','LineWidth',1.4,'DisplayName','plan');
yline(ax(2),0,'k-');  yline(ax(4),0,'k-');
ylabel(ax(1),'\theta [deg]');            legend(ax(1),'Location','southeast');
ylabel(ax(2),'\theta command - plan [deg]');
ylabel(ax(3),'pusher [rad/s]');
ylabel(ax(4),'V  (>0 = outside)');  xlabel(ax(4),'time [s]');
title(ax(1),'front dip vs the cost of using \theta');
ylim(ax(4), [-1.2 0.4]);
linkaxes(ax,'x');  xlim(ax(1),[0 20]);
exportgraphics(f, fullfile(root,'logger','front_dip.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','front_dip.fig'));
fprintf('\nsaved logger/front_dip(.fig/.png)\n');

% =========================================================================
function R = run_one(params, Ki, Kx, wscale, dt)
root = fileparts(fileparts(mfilename('fullpath')));
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;
bc.LON.W(12,12,:,:) = bc.LON.W(12,12,:,:) * wscale;
guam.reset();

rt = guam.refTraj;  n = size(rt.pos,2);
th=zeros(1,n); thv=zeros(1,n); thp=zeros(1,n); V=zeros(1,n);
pu=zeros(1,n); u=zeros(1,n); ep=zeros(1,n);
R.blown = false;
for i = 1:n
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, R.blown = true; return; end
    e  = brtV.trim_at(rt.vel(1,i));
    Rm = RSLQR.rotm_i2b(s(7),s(8),s(9));
    eb = Rm*(rt.pos(:,i) - s(1:3));
    th(i)=s(8); thp(i)=e(4); u(i)=s(4); ep(i)=eb(1);
    thv(i) = bc.theta_v;                    % allocator's theta command
    pu(i)  = guam.engineDynamics.pos(9);    % pusher speed
    V(i)   = brtV.value(s([4 6 11 8]), rt.vel(1,i));
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
k = rt.time <= traj.t_node(end) + 3;
R.t=rt.time(k); R.th=th(k); R.thv=thv(k); R.thp=thp(k);
R.V=V(k); R.pusher=pu(k); R.u=u(k);
R.dr = [0 diff(rt.vel(1,k))]/dt;
R.eu = max(abs(R.u - rt.vel(1,k)));  R.ep = max(abs(ep(k)));

% the front dip: the window from the hold to a quarter of the way through
t0 = traj.t_node(1);  t1 = t0 + 0.25*(traj.t_node(end)-t0);
m  = R.t >= t0 & R.t <= t1;
R.dth = nan(size(R.t));  R.dth(m) = R.th(m) - R.thp(m);
R.Vmax = max(R.V(m));    R.tout = dt*nnz(R.V(m) >= 0);
end
