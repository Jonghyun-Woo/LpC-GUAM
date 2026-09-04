% compare_governors - baseline vs option-1 vs the G&K reference governor,
% all scored with the SAME metric.
%
% The earlier scripts each used their own notion of "violation" (V against
% the current anchor, V against the next anchor, ...), which makes the
% numbers incomparable. Here every run is replayed through one metric:
%
%     V(x(t), v(t))  with BRTValue     - the value at the state actually
%                                        flown, for the command actually
%                                        applied at that instant
%
% That is exactly the constraint set C of the reference-governor formulation,
% so "fraction of time with V >= 0" is a fair scoreboard for all three.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
cfg  = Config('trim_schedule', params);

R = struct('name', {}, 'log', {});

%% 1) ungoverned baseline
fprintf('running baseline ...\n');
R(1).name = 'ungoverned';
R(1).log  = run_baseline(cfg, dt_sim);

%% 2) option-1 governor (present-tense gate + re-exit prediction)
fprintf('running option-1 ...\n');
brt_cell = brtV.brt;  gv = brtV.gv;
R(2).name = 'option-1 (anchor gate)';
R(2).log  = run_opt1(traj, brt_cell, gv, cfg, dt_sim);

%% 3) G&K reference governor
fprintf('running RG (slow: rollout prediction) ...\n');
R(3).name = 'RG (G&K 2002)';
R(3).log  = run_rg(traj, brtV, cfg, dt_sim);

%% Score everything with the same metric
fprintf('\n%-24s | mission | V>=0  | worst V | theta range      | alt dev\n', 'governor');
fprintf('%-24s |   [s]   |  [%%]  |         |      [deg]       |  [ft]\n', '');
fprintf('%s\n', repmat('-', 1, 92));
for c = 1:numel(R)
    L = R(c).log;
    n = numel(L.t);
    Vt = zeros(1, n);
    for i = 1:n
        Vt(i) = brtV.value(L.st([4 6 11 8], i), L.v(i));
    end
    R(c).V = Vt;
    fprintf('%-24s | %6.2f  | %5.1f | %+7.3f | %+6.2f .. %+6.2f | %6.1f\n', ...
        R(c).name, L.t_done, 100*nnz(Vt >= 0)/n, max(Vt), ...
        rad2deg(min(L.st(8,:))), rad2deg(max(L.st(8,:))), ...
        max(abs(-L.st(3,:) - 80)));
end
fprintf(['\nV >= 0 means the state has left the BRT of the command being\n' ...
         'applied at that instant - i.e. the liveness certificate is lost.\n']);

%% Figure
f = figure('Name','governor comparison','Position',[50 40 1250 800]);
col = [0.55 0.55 0.55; 0.20 0.40 0.90; 0.85 0.10 0.10];
ax(1) = subplot(2,1,1); hold on; grid on;
plot(traj.time, traj.lon(1,:), 'k--','LineWidth',1.0,'DisplayName','schedule r(t)');
for c = 1:numel(R)
    plot(R(c).log.t, R(c).log.v, '-', 'Color', col(c,:), 'LineWidth',1.5, ...
         'DisplayName', sprintf('%s (%.1f s)', R(c).name, R(c).log.t_done));
end
ylabel('applied command v [ft/s]'); legend('Location','southeast');
title('Applied command');

ax(2) = subplot(2,1,2); hold on; grid on;
for c = 1:numel(R)
    plot(R(c).log.t, R(c).V, '-', 'Color', col(c,:), 'LineWidth',1.4, ...
         'DisplayName', R(c).name);
end
yline(0,'k-','LineWidth',1.2,'HandleVisibility','off');
ylim([-0.5 0.8]); xlabel('time [s]'); ylabel('V(x, v)');
title('Liveness value at the applied command (V \geq 0 = certificate lost)');
legend('Location','northeast');
linkaxes(ax,'x');

% =========================================================================
function L = run_baseline(cfg, dt)
guam = LpC_GUAM(cfg);  rt = guam.refTraj;  M = size(rt.pos,2);
guam.reset();
L = blank(M);
for i = 1:M
    L.t(i) = rt.time(i);  L.st(:,i) = guam.state;  L.v(i) = rt.vel(1,i);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
L.t_done = rt.time(end);
end

function L = run_opt1(traj, brt, gv, cfg, dt)
guam = LpC_GUAM(cfg);
gov  = TrimProgressGovernor(traj, brt, gv, dt, 0.05);
guam.reset();  gov.reset();
Ng = round(90/dt);  L = blank(Ng);  t_done = NaN;
for i = 1:Ng
    tt = (i-1)*dt;
    [ref, gi] = gov.step(guam.state, tt);
    L.t(i) = tt;  L.st(:,i) = guam.state;  L.v(i) = gi.u_ref;
    if isnan(t_done) && gi.progress >= traj.n_trim, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done + 3, break; end
end
L = trim_log(L, i);  L.t_done = t_done;
end

function L = run_rg(traj, brtV, cfg, dt)
guam = LpC_GUAM(cfg);
gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
         struct('T',6.0,'delta',0.3,'eps',0.02,'M',10));
guam.reset();  gov.reset();
N_dec = 10;                       % 10 Hz governor
Ng = round(90/dt);  L = blank(Ng);  t_done = NaN;
v_final = traj.trim_lon(1,end);
for i = 1:Ng
    tt = (i-1)*dt;
    if mod(i-1,N_dec) == 0, ref = gov.step(tt); else, ref = gov.hold(); end
    L.t(i) = tt;  L.st(:,i) = guam.state;  L.v(i) = gov.v;
    if isnan(t_done) && gov.v >= v_final - 1e-6, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done + 3, break; end
end
L = trim_log(L, i);  L.t_done = t_done;
end

function L = blank(N)
L = struct('t',zeros(1,N),'st',zeros(12,N),'v',zeros(1,N),'t_done',NaN);
end
function L = trim_log(L, n)
L.t = L.t(1:n);  L.st = L.st(:,1:n);  L.v = L.v(1:n);
end
