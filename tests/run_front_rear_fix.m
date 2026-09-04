% run_front_rear_fix - the two dips have different causes and different
% remedies. This runs all four combinations so the two are not confused.
%
%   front dip : commanded nose-down to accelerate. Scales with the demanded
%               forward acceleration. Remedy = round the corners of the
%               speed schedule (tau), which does not lengthen the mission.
%   rear dip  : the altitude loop over-correcting a 30 s sag. Remedy = give
%               the controller the level-flight vertical speed directly
%               (w_ff) instead of making it out of position error.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q  = [0.08 0.01 1000 2.0 0 0]';
params = struct('filter_mode','off');  params.T_seg = 2.0;
CASES = { 0, 0, 'neither (shipped)'
          3, 0, 'rounded schedule only'
          0, 2, 'w feedforward only'
          3, 2, 'both' };

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

nc = size(CASES,1);
S = cell(1,nc);
fprintf('\nq1 = 0.08, q4 = 2.0, k_pos = 0.1, T_seg = 2.0 (mission length identical)\n\n');
fprintf('%-24s | %9s %9s | %9s %9s | %9s %9s %8s\n', 'case', ...
        'front V','front s','rear V','rear s', 'max|u-r|','e_pos','out %%');
fprintf('%s\n', repmat('-',1,98));
for a = 1:nc
    S{a} = run_one(params, traj, brtV, Ki, Kx, CASES{a,1}, CASES{a,2}, dt);
    r = S{a};
    fprintf('%-24s | %9.3f %9.2f | %9.3f %9.2f | %9.2f %9.1f %8.2f\n', ...
            CASES{a,3}, r.fV, r.fs, r.rV, r.rs, r.eu, r.ep, r.viol);
end
fprintf(['\n"front"/"rear" are the worst V and the seconds outside the tube in\n' ...
         't = 5..17 s and t = 30..44 s. V < 0 means the tube was never left.\n']);

%% figure
f = figure('Name','front and rear remedies','Position',[60 60 1150 820]);
ax = gobjects(1,4);
for a = 1:4, ax(a) = subplot(4,1,a); hold on; grid on; end
co = lines(nc);
for a = 1:nc
    r = S{a};
    plot(ax(1), r.t, r.V, '-','Color',co(a,:),'LineWidth',1.5,'DisplayName',CASES{a,3});
    plot(ax(2), r.t, rad2deg(r.th),'-','Color',co(a,:),'LineWidth',1.5);
    plot(ax(3), r.t, r.alt,'-','Color',co(a,:),'LineWidth',1.5);
    plot(ax(4), r.t, r.u - r.r,'-','Color',co(a,:),'LineWidth',1.5);
end
plot(ax(2), S{1}.t, rad2deg(S{1}.thp),'k--','LineWidth',1.3,'DisplayName','plan');
yline(ax(1),0,'k-','HandleVisibility','off');
yline(ax(3),80,'k--');  yline(ax(4),0,'k-');
ylabel(ax(1),'V  (>0 = outside)');  legend(ax(1),'Location','southwest');
ylabel(ax(2),'\theta [deg]');  ylabel(ax(3),'altitude [ft]');
ylabel(ax(4),'speed error [ft/s]');  xlabel(ax(4),'time [s]');
title(ax(1),'the front dip and the rear dip need different remedies');
ylim(ax(1),[-0.8 0.3]);
linkaxes(ax,'x');  xlim(ax(1),[0 S{1}.t(end)]);
exportgraphics(f, fullfile(root,'logger','front_rear_fix.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','front_rear_fix.fig'));
fprintf('saved logger/front_rear_fix(.fig/.png)\n');

% =========================================================================
function R = run_one(params, traj, brtV, Ki, Kx, tau, wff, dt)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.w_ff_mode = wff;
rt = guam.refTraj;
if tau > 0
    n = max(3, 2*floor(tau/dt/2)+1);
    v = rt.vel(1,:);
    pad = [repmat(v(1),1,(n-1)/2), v, repmat(v(end),1,(n-1)/2)];
    rt.vel(1,:) = conv(pad, ones(1,n)/n, 'valid');
    rt.pos(1,:) = cumtrapz(rt.time, rt.vel(1,:));
    guam.refTraj = rt;
end
guam.reset();

m = size(rt.pos,2);
th=zeros(1,m); thp=zeros(1,m); V=zeros(1,m); u=zeros(1,m); z=zeros(1,m); ep=zeros(1,m);
for i = 1:m
    s = guam.state;
    e  = brtV.trim_at(rt.vel(1,i));
    Rm = RSLQR.rotm_i2b(s(7),s(8),s(9));
    eb = Rm*(rt.pos(:,i) - s(1:3));
    th(i)=s(8); thp(i)=e(4); u(i)=s(4); z(i)=s(3); ep(i)=eb(1);
    V(i) = brtV.value(s([4 6 11 8]), rt.vel(1,i));
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
k = rt.time <= traj.t_node(end) + 3;
R.t=rt.time(k); R.th=th(k); R.thp=thp(k); R.V=V(k); R.u=u(k);
R.alt=-z(k); R.r=rt.vel(1,k);
R.eu = max(abs(R.u-R.r));  R.ep = max(abs(ep(k)));
R.viol = 100*nnz(R.V>=0)/numel(R.V);
fw = R.t>=5 & R.t<=17;   rw = R.t>=30 & R.t<=44;
R.fV = max(R.V(fw));  R.fs = dt*nnz(R.V(fw)>=0);
R.rV = max(R.V(rw));  R.rs = dt*nnz(R.V(rw)>=0);
end
