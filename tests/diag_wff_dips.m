% diag_wff_dips - split the BRT departures into the front and the rear
% theta dip and measure each separately, with and without the level-flight
% vertical-speed term.
%
% Reports, per contiguous departure run: when it starts, how long it lasts,
% the worst V, and the speed and attitude there. Also reports the size of
% the feedforward itself at those instants, since the term is proportional
% to u and so cannot be equally strong at both ends of the transition.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
MODE = [0 2];
NM   = {'off (shipped)', 'u*tan(theta) - w_trim'};
Q    = [0.08 0.01 1000 2.0 0 0]';
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

S = cell(1,2);
for p = 1:2
    S{p} = run_mission(params, traj, brtV, Ki, Kx, MODE(p));
end

for p = 1:2
    fprintf('\n=== %s ===\n', NM{p});
    r = S{p};
    [s, e] = runs_of(r.V >= 0);
    if isempty(s), fprintf('  no departures\n'); end
    fprintf('%8s %8s %8s | %8s %8s %8s | %9s\n', ...
            't start','t end','dur [s]','max V','u','theta','w_ff');
    fprintf('%s\n', repmat('-',1,66));
    for k = 1:numel(s)
        idx = s(k):e(k);
        [mv, j] = max(r.V(idx));  j = idx(j);
        fprintf('%8.2f %8.2f %8.2f | %8.3f %8.1f %8.1f | %9.2f\n', ...
                r.t(s(k)), r.t(e(k)), r.t(e(k))-r.t(s(k))+dt, ...
                mv, r.u(j), rad2deg(r.th(j)), r.wff(j));
    end
end

% the two dips, by the time windows the departures fall in
W = {[8 20], [30 44]};
LB = {'front dip', 'rear dip'};
fprintf('\n\nper dip: worst V and how long outside\n\n');
fprintf('%-24s', 'w feedforward');
for d = 1:2, fprintf(' | %11s %8s', [LB{d} ' V'], 'out [s]'); end
fprintf('\n%s\n', repmat('-',1,70));
for p = 1:2
    fprintf('%-24s', NM{p});
    for d = 1:2
        m = S{p}.t >= W{d}(1) & S{p}.t <= W{d}(2);
        fprintf(' | %11.3f %8.2f', max(S{p}.V(m)), dt*nnz(S{p}.V(m) >= 0));
    end
    fprintf('\n');
end

% how big is the term itself, over the mission
fprintf('\n\nsize of the feedforward term over the mission\n\n');
r = S{2};
for d = 1:2
    m = r.t >= W{d}(1) & r.t <= W{d}(2);
    fprintf('%-12s  u = %5.1f..%5.1f ft/s   |w_ff| up to %5.2f ft/s\n', ...
            LB{d}, min(r.u(m)), max(r.u(m)), max(abs(r.wff(m))));
end

%% figure
f = figure('Name','where the departures are','Position',[60 60 1150 820]);
ax = gobjects(1,4);
for a = 1:4, ax(a) = subplot(4,1,a); hold on; grid on; end
co = {[.85 .2 .2],[0 .35 .75]};
for p = 1:2
    r = S{p};
    plot(ax(1), r.t, r.V, '-','Color',co{p},'LineWidth',1.4,'DisplayName',NM{p});
    plot(ax(2), r.t, rad2deg(r.th), '-','Color',co{p},'LineWidth',1.4);
    plot(ax(3), r.t, r.alt, '-','Color',co{p},'LineWidth',1.4);
    plot(ax(4), r.t, r.wff, '-','Color',co{p},'LineWidth',1.4);
end
plot(ax(2), S{1}.t, rad2deg(S{1}.th_plan), 'k--','LineWidth',1.2);
yline(ax(1),0,'k-','HandleVisibility','off');
yline(ax(3),80,'k--');  yline(ax(4),0,'k-');
for d = 1:2
    for a = 1:4
        xline(ax(a), W{d}(1), 'Color',[.6 .6 .6],'HandleVisibility','off');
        xline(ax(a), W{d}(2), 'Color',[.6 .6 .6],'HandleVisibility','off');
    end
end
ylabel(ax(1),'V  (>0 = outside)');  legend(ax(1),'Location','northwest');
ylabel(ax(2),'\theta [deg]');  ylabel(ax(3),'altitude [ft]');
ylabel(ax(4),'w_{ff} [ft/s]');  xlabel(ax(4),'time [s]');
title(ax(1),'BRT value, attitude, altitude, and the feedforward term');
ylim(ax(1), [-1.5 0.5]);
linkaxes(ax,'x');  xlim(ax(1),[0 S{1}.t(end)]);
exportgraphics(f, fullfile(root,'logger','wff_dips.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','wff_dips.fig'));
fprintf('\nsaved logger/wff_dips(.fig/.png)\n');

% =========================================================================
function R = run_mission(params, traj, brtV, Ki, Kx, mode)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.w_ff_mode = mode;
guam.reset();
rt = guam.refTraj;  n = size(rt.pos,2);
V=zeros(1,n); u=zeros(1,n); th=zeros(1,n); z=zeros(1,n);
wff=zeros(1,n); thp=zeros(1,n);
for i = 1:n
    s = guam.state;
    u(i)=s(4); th(i)=s(8); z(i)=s(3);
    V(i) = brtV.value(s([4 6 11 8]), rt.vel(1,i));
    e = brtV.trim_at(rt.vel(1,i));  thp(i) = e(4);
    % the term as the controller computes it: u*tan(theta) minus the trim w
    wff(i) = s(4)*tan(s(8)) - e(2);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
k = rt.time <= traj.t_node(end) + 3;
R.t=rt.time(k); R.V=V(k); R.u=u(k); R.th=th(k); R.alt=-z(k);
R.wff=wff(k); R.th_plan=thp(k);
end

function [s, e] = runs_of(mask)
d = diff([false, mask(:)', false]);
s = find(d==1);  e = find(d==-1)-1;
end
