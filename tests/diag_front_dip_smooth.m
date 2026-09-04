% diag_front_dip_smooth - the front dip scales with the DEMANDED forward
% acceleration. Stretching T_seg lowers that demand but lengthens the whole
% mission. This asks whether the same relief is available for free by only
% rounding off the corners of the speed schedule.
%
% The schedule steps its slope at every trim node (the hold-to-first-segment
% step at t = t_hold is the biggest). A centred moving average of width tau
% on the speed turns each of those steps into a tau-long ramp. It preserves
% the endpoints and the total distance, so the mission takes exactly as long
% as before.
%
% tau = 0 is the shipped schedule.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q  = [0.08 0.01 1000 2.0 0 0]';
TAU = [0 1.0 2.0 3.0 4.0];
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
g0   = LpC_GUAM(Config('trim_schedule', params));
bc0  = g0.controller.baseline_controller;
cfg  = RSLQRConfig;
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

fprintf('\nT_seg = 2.0 throughout; only the corners of the speed schedule are rounded\n\n');
fprintf('%7s %10s | %11s %9s %9s | %9s %9s %9s\n', 'tau [s]','peak du/dt', ...
        'front dip','max V','out [s]','max|u-r|','e_pos','whole run');
fprintf('%s\n', repmat('-',1,86));

S = cell(1,numel(TAU));
for a = 1:numel(TAU)
    S{a} = run_one(params, traj, brtV, Ki, Kx, TAU(a), dt);
    r = S{a};
    fprintf('%7.1f %10.3f | %8.2f deg %9.3f %9.2f | %9.2f %9.1f %8.2f %%\n', ...
            TAU(a), r.peak_a, rad2deg(r.dip), r.Vmax, r.tout, r.eu, r.ep, r.viol);
end
fprintf('\n"whole run" is the percentage of the whole mission spent outside the tube.\n');

%% figure
f = figure('Name','rounding the speed schedule','Position',[60 60 1150 820]);
ax = gobjects(1,4);
for a = 1:4, ax(a) = subplot(4,1,a); hold on; grid on; end
co = lines(numel(TAU));
for a = 1:numel(TAU)
    r = S{a};
    nm = sprintf('\\tau = %.0f s', TAU(a));
    if TAU(a)==0, nm = [nm '  (shipped)']; end %#ok<AGROW>
    plot(ax(1), r.t, r.r,  '-','Color',co(a,:),'LineWidth',1.4,'DisplayName',nm);
    plot(ax(2), r.t, r.dr, '-','Color',co(a,:),'LineWidth',1.4);
    plot(ax(3), r.t, rad2deg(r.th), '-','Color',co(a,:),'LineWidth',1.4);
    plot(ax(4), r.t, r.V,  '-','Color',co(a,:),'LineWidth',1.4);
end
plot(ax(3), S{1}.t, rad2deg(S{1}.thp), 'k--','LineWidth',1.4,'DisplayName','plan');
yline(ax(4),0,'k-');
ylabel(ax(1),'speed command [ft/s]');  legend(ax(1),'Location','southeast');
ylabel(ax(2),'du/dt [ft/s^2]');
ylabel(ax(3),'\theta [deg]');
ylabel(ax(4),'V  (>0 = outside)');  xlabel(ax(4),'time [s]');
title(ax(1),'rounding the corners of the speed schedule (mission length unchanged)');
ylim(ax(4),[-1.2 0.4]);  ylim(ax(2),[-1 6]);
linkaxes(ax,'x');  xlim(ax(1),[0 25]);
exportgraphics(f, fullfile(root,'logger','front_dip_smooth.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','front_dip_smooth.fig'));
fprintf('saved logger/front_dip_smooth(.fig/.png)\n');

% =========================================================================
function R = run_one(params, traj, brtV, Ki, Kx, tau, dt)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;

rt = guam.refTraj;
if tau > 0
    n = max(3, 2*floor(tau/dt/2)+1);          % odd window, width tau
    v = rt.vel(1,:);
    pad = [repmat(v(1),1,(n-1)/2), v, repmat(v(end),1,(n-1)/2)];
    rt.vel(1,:) = conv(pad, ones(1,n)/n, 'valid');
    rt.pos(1,:) = cumtrapz(rt.time, rt.vel(1,:));
    guam.refTraj = rt;
end
guam.reset();

m = size(rt.pos,2);
th=zeros(1,m); thp=zeros(1,m); V=zeros(1,m); u=zeros(1,m); ep=zeros(1,m);
for i = 1:m
    s = guam.state;
    e  = brtV.trim_at(rt.vel(1,i));
    Rm = RSLQR.rotm_i2b(s(7),s(8),s(9));
    eb = Rm*(rt.pos(:,i) - s(1:3));
    th(i)=s(8); thp(i)=e(4); u(i)=s(4); ep(i)=eb(1);
    V(i) = brtV.value(s([4 6 11 8]), rt.vel(1,i));
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
k = rt.time <= traj.t_node(end) + 3;
R.t=rt.time(k); R.th=th(k); R.thp=thp(k); R.V=V(k); R.u=u(k); R.r=rt.vel(1,k);
R.dr = [0 diff(R.r)]/dt;
R.peak_a = max(R.dr);
R.eu = max(abs(R.u-R.r));  R.ep = max(abs(ep(k)));
R.viol = 100*nnz(R.V >= 0)/numel(R.V);

t0 = traj.t_node(1);  t1 = t0 + 0.25*(traj.t_node(end)-t0);
w = R.t >= t0 & R.t <= t1;
R.dip  = min(R.th(w) - R.thp(w));
R.Vmax = max(R.V(w));
R.tout = dt*nnz(R.V(w) >= 0);
end
