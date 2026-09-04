% diag_theta_dip - why does pitch dive the wrong way and come back?
%
% The 3-D comparison showed the flown theta digging well below the scheduled
% trim attitude and recovering. Tracking is what matters here, so this plots
% the time histories that could explain it: the commanded vs flown states,
% and every effector against its limit. If something saturates during the
% dip, that is the cause; if nothing does, it is the loop itself.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  Q1 = 0.08;
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));

g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

Q = [Q1 0.01 1000 0 0 0]';
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

guam = LpC_GUAM(Config('trim_schedule', params));
guam.controller.baseline_controller.LON.Ki = Ki;
guam.controller.baseline_controller.LON.Kx = Kx;
guam.reset();

rt = guam.refTraj;  n = size(rt.pos,2);
L = struct('t',rt.time,'u',zeros(1,n),'w',zeros(1,n),'th',zeros(1,n), ...
           'q',zeros(1,n),'r',rt.vel(1,:),'th_c',zeros(1,n),'w_c',zeros(1,n), ...
           'eng',zeros(9,n),'srf',zeros(5,n),'epos',zeros(1,n));
for i = 1:n
    s = guam.state;
    [X0, ~] = guam.controller.interp_xu0(rt.vel(1,i), rt.vel(3,i));
    Rm = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e  = Rm * (rt.pos(:,i) - s(1:3));
    L.u(i)=s(4); L.w(i)=s(6); L.th(i)=s(8); L.q(i)=s(11);
    L.th_c(i) = X0(11);              % trim pitch at the commanded speed
    L.w_c(i)  = X0(3);               % trim vertical speed at the commanded speed
    L.epos(i) = e(1);
    L.eng(:,i) = guam.engineDynamics.pos;
    L.srf(:,i) = guam.surfaceDynamics.pos;
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
keep = L.t <= traj.t_node(end) + 3;
fn = fieldnames(L);
for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:,keep); end
t = L.t;

%% where is the worst pitch error
eth = rad2deg(L.th - L.th_c);
[~, id] = max(abs(eth));
fprintf('\nworst pitch error %.2f deg at t = %.2f s (u = %.1f ft/s)\n', ...
        eth(id), t(id), L.u(id));
fprintf('   flown theta %.2f deg,  trim theta %.2f deg\n', ...
        rad2deg(L.th(id)), rad2deg(L.th_c(id)));
fprintf('   flown w %.2f ft/s,  trim w %.2f ft/s,  commanded w 0.00\n', ...
        L.w(id), L.w_c(id));
fprintf('   elevator %.2f deg (limit +/-30),  flap %.2f deg\n', ...
        rad2deg(L.srf(3,id)), rad2deg(mean(L.srf(1:2,id))));
fprintf('   lift rotors %.0f..%.0f rpm (limit 1600),  pusher %.0f rpm (limit 2000)\n', ...
        min(L.eng(1:8,id))*60/2/pi, max(L.eng(1:8,id))*60/2/pi, ...
        L.eng(9,id)*60/2/pi);

%% figure
f = figure('Name','pitch dip diagnosis','Position',[60 30 1200 950]);
ax = gobjects(1,5);

ax(1) = subplot(5,1,1); hold on; grid on;
plot(t, L.r, 'k--','LineWidth',1.3,'DisplayName','commanded u');
plot(t, L.u, '-','Color',[.85 .2 .2],'LineWidth',1.5,'DisplayName','flown u');
ylabel('u [ft/s]'); legend('Location','southeast');
title(sprintf('nominal loop, q_1 = %.2f, no governor', Q1));

ax(2) = subplot(5,1,2); hold on; grid on;
plot(t, L.u - L.r, '-','Color',[.85 .2 .2],'LineWidth',1.4);
yline(0,'k-'); ylabel('u error [ft/s]');

ax(3) = subplot(5,1,3); hold on; grid on;
plot(t, rad2deg(L.th_c), 'k--','LineWidth',1.3,'DisplayName','trim \theta at command');
plot(t, rad2deg(L.th),  '-','Color',[0 .45 .75],'LineWidth',1.5,'DisplayName','flown \theta');
plot(t(id), rad2deg(L.th(id)), 'rv','MarkerFaceColor','r','HandleVisibility','off');
ylabel('\theta [deg]'); legend('Location','southwest');

ax(4) = subplot(5,1,4); hold on; grid on;
plot(t, L.w_c, 'k--','LineWidth',1.3,'DisplayName','trim w at command');
plot(t, L.w,  '-','Color',[.4 .2 .7],'LineWidth',1.5,'DisplayName','flown w');
yline(0,'-','Color',[.6 .6 .6],'DisplayName','commanded w = 0');
ylabel('w [ft/s]'); legend('Location','southwest');

ax(5) = subplot(5,1,5); hold on; grid on;
yyaxis left
plot(t, rad2deg(L.srf(3,:)), '-','LineWidth',1.4,'DisplayName','elevator');
yline( 30,'r:','HandleVisibility','off');  yline(-30,'r:','HandleVisibility','off');
ylabel('elevator [deg]'); ylim([-35 35]);
yyaxis right
plot(t, mean(L.eng(1:8,:),1)*60/2/pi, '-','LineWidth',1.4,'DisplayName','lift rotors');
plot(t, L.eng(9,:)*60/2/pi, '-','LineWidth',1.4,'DisplayName','pusher');
yline(1600,'r:','HandleVisibility','off');
ylabel('rotor [rpm]');
xlabel('time [s]'); legend('Location','northwest');

linkaxes(ax,'x'); xlim(ax(1),[0 t(end)]);
exportgraphics(f, fullfile(root,'logger','theta_dip.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','theta_dip.fig'));
fprintf('\nsaved logger/theta_dip(.fig/.png)\n');
