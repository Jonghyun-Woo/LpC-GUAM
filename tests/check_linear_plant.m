% check_linear_plant - is the stored linearisation good enough to predict with?
%
% Step 1 of replacing the nonlinear rollout. This isolates the PLANT
% linearisation from the controller: the nonlinear rollout is run once, its
% effector history is recorded, and that same history is then fed to
%
%   (a) a frozen LTI model      Ap, Bp taken once at the starting trim
%   (b) a scheduled LTV model   Ap, Bp re-interpolated at the current speed
%
% If (b) tracks the nonlinear states over the 6 s horizon, the rollout can be
% replaced. If it does not, no amount of closed-loop assembly will fix it, so
% this has to be answered first.
%
% Xlon = [u; w; q; theta] perturbation from trim;  11 lon effectors are
% U0_idx = [5..13 (8 lift + pusher), 3 (elevator), 1 (flap)].
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  T_HOR = 6.0;  N = round(T_HOR/dt);
T_START = 20;                       % where in the mission to test
% set T_START from the sweep list below; the script loops over these
T_LIST  = [10 20 30 40];
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));
gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
        struct('T',6.0,'delta',0.3,'eps',0.02,'M',40));
guam.reset();  gov.reset();

S  = load(fullfile(root,'tables','trim','trim_table_Poly_ConcatVer4p0.mat'), ...
          'Ap_lon_interp','Bp_lon_interp','UH','WH','XU0_interp');
U0_idx = [5 6 7 8 9 10 11 12 13 3 1];

%% Fly to the test point
N_dec = 10;
for i = 1:round(T_START/dt)
    tt = (i-1)*dt;
    if mod(i-1,N_dec)==0, ref = gov.step(tt); else, ref = gov.hold(); end
    guam.step(ref);
end
v_hold = gov.v;  p0 = gov.pos_n;
fprintf('test point: t = %.1f s, command held at %.1f ft/s\n', T_START, v_hold);

%% (0) Nonlinear rollout, recording states AND effectors
snap = guam.saveState();
XN = zeros(4,N+1);  UE = zeros(11,N);
XN(:,1) = guam.state([4 6 11 8]);
p = p0;
for k = 1:N
    p = p + v_hold*dt;
    guam.step(struct('pos',[p;0;-80],'vel',[v_hold;0;0],'chi',0,'chi_dot',0));
    e = guam.engineDynamics.pos;  s = guam.surfaceDynamics.pos;
    UE(:,k)   = [e(:); s(3); mean(s(1:2))];       % 9 props, elevator, flap
    XN(:,k+1) = guam.state([4 6 11 8]);
end
guam.restoreState(snap);

%% (a) frozen LTI and (b) scheduled LTV, driven by the SAME effectors
XA = zeros(4,N+1);  XB = zeros(4,N+1);
XA(:,1) = XN(:,1);  XB(:,1) = XN(:,1);
[Ap0,Bp0,xe0,ue0] = sched(S, U0_idx, v_hold, 0);
for k = 1:N
    % (a) everything frozen at the starting trim
    dxa = XA(:,k) - xe0;
    XA(:,k+1) = XA(:,k) + dt*(Ap0*dxa + Bp0*(UE(:,k) - ue0));

    % (b) re-scheduled on the current predicted speed
    [Ap,Bp,xe,ue] = sched(S, U0_idx, XB(1,k), 0);
    dxb = XB(:,k) - xe;
    XB(:,k+1) = XB(:,k) + dt*(Ap*dxb + Bp*(UE(:,k) - ue));
end

%% Report
nm = {'u  [ft/s]','w  [ft/s]','q  [deg/s]','th [deg]'};
sc = [1 1 180/pi 180/pi];
fprintf('\nerror against the nonlinear rollout (same effector history)\n');
fprintf('%-12s %10s %10s %10s %10s\n','state','LTI @1s','LTI @6s','LTV @1s','LTV @6s');
fprintf('%s\n', repmat('-',1,56));
i1 = round(1/dt)+1;
for d = 1:4
    fprintf('%-12s %10.3f %10.3f %10.3f %10.3f\n', nm{d}, ...
        sc(d)*(XA(d,i1)-XN(d,i1)), sc(d)*(XA(d,end)-XN(d,end)), ...
        sc(d)*(XB(d,i1)-XN(d,i1)), sc(d)*(XB(d,end)-XN(d,end)));
end

fprintf('\nBRT value at the horizon end (this is what the governor judges on)\n');
fprintf('  nonlinear %+.4f | frozen LTI %+.4f | scheduled LTV %+.4f\n', ...
    brtV.value(XN(:,end), v_hold), brtV.value(XA(:,end), v_hold), ...
    brtV.value(XB(:,end), v_hold));

%% Figure
t = (0:N)*dt;
f = figure('Name','linear vs nonlinear rollout','Position',[80 60 1100 780]);
for d = 1:4
    ax = subplot(4,1,d); hold on; grid on;
    plot(ax, t, sc(d)*XN(d,:), 'k-',  'LineWidth',2.0, 'DisplayName','nonlinear (truth)');
    plot(ax, t, sc(d)*XA(d,:), '--',  'Color',[.85 .2 .2], 'LineWidth',1.5, 'DisplayName','frozen LTI');
    plot(ax, t, sc(d)*XB(d,:), '-',   'Color',[0 .45 .75], 'LineWidth',1.5, 'DisplayName','scheduled LTV');
    ylabel(ax, nm{d});
    if d==1, legend(ax,'Location','best');
        title(ax, sprintf('6 s rollout from t = %.0f s, command held at %.1f ft/s', T_START, v_hold)); end
    if d==4, xlabel(ax,'prediction time [s]'); end
end
exportgraphics(f, fullfile(root,'logger','linear_vs_nonlinear.png'), 'Resolution',150);
savefig(f, fullfile(root,'logger','linear_vs_nonlinear.fig'));

% -------------------------------------------------------------------------
function [Ap,Bp,xe,ue] = sched(S, U0_idx, uh, wh)
% Bilinear interpolation of the stored linearisation and trim at (uh, wh).
uh = min(max(uh, S.UH(1)), S.UH(end));
wh = min(max(wh, S.WH(1)), S.WH(end));
iu = find(S.UH <= uh, 1, 'last');  iu = min(max(iu,1), numel(S.UH)-1);
iw = find(S.WH <= wh, 1, 'last');  iw = min(max(iw,1), numel(S.WH)-1);
ku = (uh - S.UH(iu))/(S.UH(iu+1) - S.UH(iu));
kw = (wh - S.WH(iw))/(S.WH(iw+1) - S.WH(iw));
b = @(M) (1-ku)*(1-kw)*M(:,:,iu,iw)   + ku*(1-kw)*M(:,:,iu+1,iw) + ...
         (1-ku)*   kw *M(:,:,iu,iw+1) + ku*   kw *M(:,:,iu+1,iw+1);
Ap = b(S.Ap_lon_interp);
Bp = b(S.Bp_lon_interp);
X  = (1-ku)*(1-kw)*S.XU0_interp(:,iu,iw)   + ku*(1-kw)*S.XU0_interp(:,iu+1,iw) + ...
     (1-ku)*   kw *S.XU0_interp(:,iu,iw+1) + ku*   kw *S.XU0_interp(:,iu+1,iw+1);
xe = X([1 3 5 11]);              % [u; w; q; theta] trim
ue = X(12 + U0_idx(:));          % 11 lon effector trims
end
