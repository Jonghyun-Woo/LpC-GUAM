% sweep_linear_horizon - how long can a linear model predict before the
% governor would make a different decision?
%
% The governor never looks at states directly; it looks at V. So the quantity
% that matters is |V_linear - V_nonlinear| as a function of the horizon. If
% that stays below the governor's own margin (eps = 0.02) out to some T, the
% rollout can be replaced by a linear propagation with that T.
%
% Run at several points in the mission because the low-speed segments are
% where the governor actually spends its time.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  T_HOR = 6.0;  N = round(T_HOR/dt);
T_LIST = [10 15 20 25 30 40];
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));
gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
        struct('T',6.0,'delta',0.3,'eps',0.02,'M',40));
guam.reset();  gov.reset();

S = load(fullfile(root,'tables','trim','trim_table_Poly_ConcatVer4p0.mat'), ...
         'Ap_lon_interp','Bp_lon_interp','UH','WH','XU0_interp');
U0_idx = [5 6 7 8 9 10 11 12 13 3 1];

probe = [0.5 1 2 3 4 6];                      % horizons to report [s]
EA = nan(numel(T_LIST), numel(probe));        % |dV| frozen LTI
EB = nan(numel(T_LIST), numel(probe));        % |dV| scheduled LTV
VH = nan(1, numel(T_LIST));  UU = nan(1, numel(T_LIST));

N_dec = 10;  i_glob = 0;  t_next = 1;
for i = 1:round(max(T_LIST)/dt)
    tt = (i-1)*dt;
    if t_next <= numel(T_LIST) && tt >= T_LIST(t_next) - 1e-9
        [EA(t_next,:), EB(t_next,:), VH(t_next), UU(t_next)] = ...
            probe_here(guam, gov, brtV, S, U0_idx, dt, N, probe);
        t_next = t_next + 1;
    end
    if mod(i-1,N_dec)==0, ref = gov.step(tt); else, ref = gov.hold(); end
    guam.step(ref);
end

%% Report
fprintf('\n|V_linear - V_nonlinear| against horizon   (governor eps = 0.020)\n');
fprintf('%6s %7s %7s |', 't[s]', 'u', 'V_true');
fprintf('%8.1fs', probe); fprintf('\n');
fprintf('%s\n', repmat('-', 1, 34 + 9*numel(probe)));
for a = 1:numel(T_LIST)
    fprintf('%6.0f %7.1f %7.3f | frozen LTI ', T_LIST(a), UU(a), VH(a));
    fprintf('%8.4f', EA(a,:)); fprintf('\n');
    fprintf('%6s %7s %7s |    sched LTV ', '', '', '');
    fprintf('%8.4f', EB(a,:)); fprintf('\n');
end

fprintf('\nlongest horizon with |dV| <= 0.020\n');
for a = 1:numel(T_LIST)
    ja = find(EA(a,:) <= 0.02, 1, 'last');  jb = find(EB(a,:) <= 0.02, 1, 'last');
    sa = 'none'; if ~isempty(ja), sa = sprintf('%.1f s', probe(ja)); end
    sb = 'none'; if ~isempty(jb), sb = sprintf('%.1f s', probe(jb)); end
    fprintf('  t = %2.0f s (u = %5.1f):  frozen LTI %-7s   scheduled LTV %s\n', ...
            T_LIST(a), UU(a), sa, sb);
end

% -------------------------------------------------------------------------
function [eA, eB, V0, u0] = probe_here(guam, gov, brtV, S, U0_idx, dt, N, probe)
v = gov.v;  u0 = guam.state(4);
snap = guam.saveState();
XN = zeros(4,N+1);  UE = zeros(11,N);
XN(:,1) = guam.state([4 6 11 8]);
p = gov.pos_n;
for k = 1:N
    p = p + v*dt;
    guam.step(struct('pos',[p;0;-80],'vel',[v;0;0],'chi',0,'chi_dot',0));
    e = guam.engineDynamics.pos;  s = guam.surfaceDynamics.pos;
    UE(:,k)   = [e(:); s(3); mean(s(1:2))];
    XN(:,k+1) = guam.state([4 6 11 8]);
end
guam.restoreState(snap);

XA = XN(:,1);  XB = XN(:,1);
[Ap0,Bp0,xe0,ue0] = sched(S, U0_idx, v, 0);
eA = zeros(1,numel(probe));  eB = zeros(1,numel(probe));
jp = 1;
for k = 1:N
    XA = XA + dt*(Ap0*(XA - xe0) + Bp0*(UE(:,k) - ue0));
    [Ap,Bp,xe,ue] = sched(S, U0_idx, XB(1), 0);
    XB = XB + dt*(Ap*(XB - xe) + Bp*(UE(:,k) - ue));
    if jp <= numel(probe) && abs(k*dt - probe(jp)) < dt/2
        Vt = brtV.value(XN(:,k+1), v);
        eA(jp) = abs(brtV.value(XA, v) - Vt);
        eB(jp) = abs(brtV.value(XB, v) - Vt);
        jp = jp + 1;
    end
end
V0 = brtV.value(XN(:,end), v);
end

function [Ap,Bp,xe,ue] = sched(S, U0_idx, uh, wh)
uh = min(max(uh, S.UH(1)), S.UH(end));
wh = min(max(wh, S.WH(1)), S.WH(end));
iu = find(S.UH <= uh, 1, 'last');  iu = min(max(iu,1), numel(S.UH)-1);
iw = find(S.WH <= wh, 1, 'last');  iw = min(max(iw,1), numel(S.WH)-1);
ku = (uh - S.UH(iu))/(S.UH(iu+1) - S.UH(iu));
kw = (wh - S.WH(iw))/(S.WH(iw+1) - S.WH(iw));
b = @(M) (1-ku)*(1-kw)*M(:,:,iu,iw)   + ku*(1-kw)*M(:,:,iu+1,iw) + ...
         (1-ku)*   kw *M(:,:,iu,iw+1) + ku*   kw *M(:,:,iu+1,iw+1);
Ap = b(S.Ap_lon_interp);  Bp = b(S.Bp_lon_interp);
X  = (1-ku)*(1-kw)*S.XU0_interp(:,iu,iw)   + ku*(1-kw)*S.XU0_interp(:,iu+1,iw) + ...
     (1-ku)*   kw *S.XU0_interp(:,iu,iw+1) + ku*   kw *S.XU0_interp(:,iu+1,iw+1);
xe = X([1 3 5 11]);
ue = X(12 + U0_idx(:));
end
