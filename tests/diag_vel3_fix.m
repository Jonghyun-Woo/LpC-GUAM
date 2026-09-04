% diag_vel3_fix - does commanding the trim's vertical speed fix the pitch error?
%
% ReferenceTrajectory builds vel = zeros(3,N) and fills only row 1, so the
% vertical-speed command is 0 for the whole mission. Its comment says row 3
% is the inertial vertical velocity, where 0 means level flight. But RSLQR
% feeds that number straight into the trim table as the body-frame w
% breakpoint, and at cruise the body w for level flight is u*tan(theta) ~ 13
% to 22 ft/s, not 0. So commanding 0 selects the WH = 0 trim, which at speed
% is a steady CLIMB, while the position reference holds 80 ft. The two fight,
% and the pitch error is the result.
%
% Test: leave everything else alone and command the trim's own w instead.
% Nothing is edited in ReferenceTrajectory - the table is patched here so the
% change is easy to undo if it does not help.
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

% retuned gains, same for both runs
Q = [Q1 0.01 1000 0 0 0]';
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

% Level flight needs body w = u*tan(theta), which ramps from 0 at hover to
% about 13 ft/s at cruise - not the constant 11.67 of the WH3 row, and not 0.
w_level = traj.lon(1,:) .* tan(traj.lon(4,:));

fprintf('\nq1 = %.2f, no governor\n\n', Q1);
fprintf('%-28s %9s %9s %9s %9s\n', 'vertical-speed command', ...
        'max|u-r|', 'max|w|', 'alt dev', 'max|dth|');
fprintf('%s\n', repmat('-', 1, 68));

W3 = {0, traj.lon(2,:), w_level};
NM = {'0  (as shipped)', 'trim w  (11.67 const)', 'u*tan(theta)  (level)'};
S  = cell(1,3);
for c = 1:3
    S{c} = run_case(params, traj, Ki, Kx, W3{c}, dt);
    fprintf('%-28s %9.2f %9.2f %9.1f %9.2f\n', NM{c}, ...
            S{c}.eu, S{c}.wmax, S{c}.ad, S{c}.eth);
end

%% figure
f = figure('Name','vertical-speed command','Position',[60 60 1150 850]);
ax = gobjects(1,4);
lab = {'w_{cmd} = 0 (shipped)','w_{cmd} = 11.67 const','w_{cmd} = u\cdottan\theta'};
col = {[.85 .2 .2], [.55 .35 .1], [0 .5 .2]};
NC  = 3;

ax(1) = subplot(4,1,1); hold on; grid on;
plot(S{1}.t, S{1}.r, 'k--','LineWidth',1.2,'DisplayName','commanded u');
for c = 1:NC, plot(S{c}.t, S{c}.u, '-','Color',col{c},'LineWidth',1.4,'DisplayName',lab{c}); end
ylabel('u [ft/s]'); legend('Location','southeast');
title(sprintf('vertical-speed command (q_1 = %.2f, no governor)', Q1));

ax(2) = subplot(4,1,2); hold on; grid on;
for c = 1:NC, plot(S{c}.t, S{c}.u - S{c}.r, '-','Color',col{c},'LineWidth',1.4); end
yline(0,'k-'); ylabel('u error [ft/s]');

ax(3) = subplot(4,1,3); hold on; grid on;
for c = 1:NC
    plot(S{c}.t, S{c}.w - S{c}.wc, '-','Color',col{c},'LineWidth',1.4);
end
yline(0,'k-'); ylabel('w error [ft/s]');

ax(4) = subplot(4,1,4); hold on; grid on;
for c = 1:NC
    plot(S{c}.t, -S{c}.z, '-','Color',col{c},'LineWidth',1.5,'DisplayName',lab{c});
end
yline(80,'k--','HandleVisibility','off');
ylabel('altitude [ft]'); xlabel('time [s]'); legend('Location','best');

linkaxes(ax,'x'); xlim(ax(1), [0 S{1}.t(end)]);
exportgraphics(f, fullfile(root,'logger','vel3_fix.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','vel3_fix.fig'));
fprintf('\nsaved logger/vel3_fix(.fig/.png)\n');

% =========================================================================
function S = run_case(params, traj, Ki, Kx, w3, dt)
guam = LpC_GUAM(Config('trim_schedule', params));
guam.controller.baseline_controller.LON.Ki = Ki;
guam.controller.baseline_controller.LON.Kx = Kx;

rt = guam.refTraj;
if isscalar(w3), rt.vel(3,:) = w3;  else, rt.vel(3,:) = w3;  end
guam.refTraj = rt;                       % reset() reads the patched table
guam.reset();

n = size(rt.pos,2);
S.t = rt.time;  S.r = rt.vel(1,:);
S.u = zeros(1,n); S.th = zeros(1,n); S.thc = zeros(1,n);
S.w = zeros(1,n); S.z = zeros(1,n); S.wc = rt.vel(3,:);
for i = 1:n
    s = guam.state;
    [X0, ~] = guam.controller.interp_xu0(rt.vel(1,i), rt.vel(3,i));
    S.u(i)=s(4); S.th(i)=s(8); S.w(i)=s(6); S.z(i)=s(3); S.thc(i)=X0(11);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
keep = S.t <= traj.t_node(end) + 3;
fn = {'t','r','u','th','thc','w','z','wc'};
for k = 1:numel(fn), S.(fn{k}) = S.(fn{k})(keep); end
S.eu   = max(abs(S.u - S.r));
S.eth  = max(abs(rad2deg(S.th - S.thc)));
S.ad   = max(abs(-S.z - 80));
S.wmax = max(abs(S.w));
end
