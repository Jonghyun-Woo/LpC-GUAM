% verify_fixed_anchors - the same design, restricted to the 20 trim points
% that actually have a BRT. No interpolated trim points anywhere.
%
% ---------------------------------------------------------------------------
% WHAT CHANGES
% ---------------------------------------------------------------------------
% The placement result (23 trim points, steps 4.4 -> 11.0) puts trim points at
% speeds like 4.23 ft/s, where neither the trim nor the reachable set was ever
% computed - both were interpolated. If interpolated trim points are not
% acceptable, the step length stops being a design variable: the only steps
% available are 8.44 and its multiples, because those are the only places a
% tube exists.
%
% Losing the step leaves exactly one knob. The condition
%
%      a_fwd(u_k) + a_bwd(u_k+1) - D - 2*( e0 + c(u,T)*D/T )  >=  0
%
% with D pinned at 8.44 still has T in it. Longer segments mean a gentler ramp,
% a smaller lag, and a smaller transient - which is the term that actually
% consumes the margin (check_tube_dimensions.m: pitch attitude uses 70 % of its
% allowance, forward speed only 41 %). So the knob that remains is the one
% acting on the binding constraint.
%
% ---------------------------------------------------------------------------
% WHAT THIS SCRIPT DOES
% ---------------------------------------------------------------------------
%  1. extends c(u,T) to long segment times (the stored grid stops at 6 s)
%  2. solves the condition for the shortest admissible T, segment by segment,
%     with the step pinned at 8.44
%  3. flies the resulting schedules and checks them against the tubes
%
% Output: logger/fixed_anchor_time.mat / .png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt     = 0.01;
WH_IDX = 3;
N_TRIM = 20;
U_ROW  = 4;
T_HOLD = 5.0;

O  = load(fullfile(root,'logger','segment_time_opt.mat'));
Lm = load(fullfile(root,'logger','closed_loop_model_ramped.mat'));  M = Lm.M;
Lc = load(fullfile(root,'logger','lag_coefficient.mat'));
u_anchor = O.u_anchor;  a_fwd = O.a_fwd;  a_bwd = O.a_bwd;  e0_u = O.e0_u;
dXE = Lc.dXE;

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj0 = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV  = BRTValue(fullfile(root,'data'), traj0.trim_lon, WH_IDX);

%% ---------------------------------------------------------------------
%  1. c(u,T) out to long T
%% ---------------------------------------------------------------------
T_LONG = [1 2 3 4 6 8 10 12 15 20 25 30];
nT = numel(T_LONG);
C = zeros(nT, N_TRIM);
for a = 1:nT
    for k = 1:N_TRIM
        C(a,k) = ramp_peak(M.Phi{k}, dXE(:,k), T_LONG(a), 8.0, dt, U_ROW, ...
                           ClosedLoopModel.N_XI);
    end
end
fprintf('--- c(u,T) at long segment times [s] ---\n');
fprintf('   T  |');  fprintf(' u=%4.0f', u_anchor([1 3 6 10 15 20]));  fprintf('\n');
for a = 1:nT
    fprintf('%5.0f |', T_LONG(a));  fprintf(' %6.3f', C(a,[1 3 6 10 15 20]));  fprintf('\n');
end
fprintf(['\nc stops growing past ~3 s - the lag has finished building - so\n' ...
         'c/T keeps falling as T grows. That is what buys the margin back.\n']);

%% ---------------------------------------------------------------------
%  2. shortest admissible T per segment, step pinned at 8.44
%% ---------------------------------------------------------------------
fprintf('\n--- shortest admissible T with the step pinned at 8.44 ---\n');
fprintf('seg |    u    | 겹침 overlap | need 2e <= | T_min [s] | g at T_min\n');
T_need = zeros(1, N_TRIM-1);
for k = 1:N_TRIM-1
    D  = u_anchor(k+1) - u_anchor(k);
    ov = a_fwd(k) + a_bwd(k+1) - D;              % the room available
    ok = NaN;  gg = NaN;
    for a = 1:nT
        e = e0_u(k) + C(a,k)*D/T_LONG(a);
        if ov - 2*e >= 0, ok = T_LONG(a);  gg = ov - 2*e;  break; end
    end
    T_need(k) = ok;
    if isnan(ok)
        fprintf('%3d | %7.2f | %12.2f | %10.2f |    ----   | never (even T=%d)\n', ...
                k, u_anchor(k), ov, ov, T_LONG(end));
    else
        fprintf('%3d | %7.2f | %12.2f | %10.2f | %9.0f | %+8.3f\n', ...
                k, u_anchor(k), ov, ov, ok, gg);
    end
end
if all(~isnan(T_need))
    fprintf('\ntotal if each segment takes its own minimum: %.0f s\n', sum(T_need));
else
    fprintf('\nsegments %s cannot be made admissible at ANY segment time.\n', ...
            mat2str(find(isnan(T_need))));
    fprintf(['This is the 1-D condition talking. It demands the whole spread of\n' ...
             'the vehicle fit inside a %0.2f ft/s window, which no amount of time\n' ...
             'fixes because e0 alone is already comparable to it. Section 3 asks\n' ...
             'the real question instead: does the vehicle actually leave a tube?\n'], ...
            a_fwd(1)+a_bwd(2)-8.44);
end

%% ---------------------------------------------------------------------
%  3. the real question: fly the 20 fixed anchors and look at the tubes
%% ---------------------------------------------------------------------
fprintf('\n--- flying the ORIGINAL 20 trim points at a range of segment times ---\n');
fprintf('%6s | %6s | %7s | %8s | %6s | %s\n', ...
        'T [s]', 'time', 'V>=0 %', 'worst V', 'at u', 'verdict');
fprintf('%s\n', repmat('-', 1, 62));
TL = [2.0 2.1 2.2 2.3 2.4 2.5 3.0 4.0];
res = struct([]);
for i = 1:numel(TL)
    [tt, uu] = profile_from(u_anchor, repmat(TL(i),1,N_TRIM-1), dt, T_HOLD);
    R = fly(params, tt, uu, dt, brtV, T_HOLD);
    if R.viol == 0, vd = 'SAFE'; else, vd = 'unsafe'; end
    fprintf('%6.2f | %6.1f | %7.1f | %+8.3f | %6.1f | %s\n', ...
            TL(i), tt(end)-2*T_HOLD, R.viol, R.Vw, R.u_worst, vd);
    res(i).T = TL(i);  res(i).t = tt(end)-2*T_HOLD;  res(i).viol = R.viol;
    res(i).Vw = R.Vw;  res(i).V = R.V;  res(i).tt = tt;  res(i).uu = uu;
end

%% ---------------------------------------------------------------------
%  4. per-segment T: slow only where it is needed
%     The g table showed segments 1-8 failing and 9-19 passing at T = 2, so
%     only the low-speed end needs the extra time.
%% ---------------------------------------------------------------------
fprintf('\n--- per-segment: slow at low speed, quick at high speed ---\n');
fprintf('%-30s | %6s | %7s | %+8s | %s\n', 'schedule', 'time', 'V>=0 %', 'worst V', 'verdict');
fprintf('%s\n', repmat('-', 1, 72));
mixes = { ...
  'T=3.0 s1-4, 2.5 s5-8, 2.0 rest', [repmat(3.0,1,4) repmat(2.5,1,4) repmat(2.0,1,11)]; ...
  'T=3.5 s1-4, 2.5 s5-8, 2.0 rest', [repmat(3.5,1,4) repmat(2.5,1,4) repmat(2.0,1,11)]; ...
  'T=4.0 s1-4, 2.5 s5-8, 2.0 rest', [repmat(4.0,1,4) repmat(2.5,1,4) repmat(2.0,1,11)]; ...
  'T=4.0 s1-4, 3.0 s5-8, 2.0 rest', [repmat(4.0,1,4) repmat(3.0,1,4) repmat(2.0,1,11)]; ...
  'T=3.0 s1-6, 2.0 rest',           [repmat(3.0,1,6) repmat(2.0,1,13)]; ...
  'T=4.0 s1-6, 2.0 rest',           [repmat(4.0,1,6) repmat(2.0,1,13)]; ...
  'T=4.0 s1-8, 2.0 rest',           [repmat(4.0,1,8) repmat(2.0,1,11)]; ...
};
best = struct('t', inf);
for i = 1:size(mixes,1)
    [tt, uu] = profile_from(u_anchor, mixes{i,2}, dt, T_HOLD);
    R = fly(params, tt, uu, dt, brtV, T_HOLD);
    if R.viol == 0, vd = 'SAFE'; else, vd = 'unsafe'; end
    fprintf('%-30s | %6.1f | %7.1f | %+8.3f | %s\n', ...
            mixes{i,1}, tt(end)-2*T_HOLD, R.viol, R.Vw, vd);
    if R.viol == 0 && tt(end)-2*T_HOLD < best.t
        best = struct('t', tt(end)-2*T_HOLD, 'name', mixes{i,1}, 'Vw', R.Vw);
    end
end
for i = 1:numel(res)
    if res(i).viol == 0 && res(i).t < best.t
        best = struct('t', res(i).t, 'name', sprintf('uniform T = %.2f', res(i).T), ...
                      'Vw', res(i).Vw);
    end
end
fprintf('%s\n', repmat('-', 1, 72));
if isfinite(best.t)
    fprintf('FASTEST SAFE WITH THE ORIGINAL 20 POINTS: %s -> %.1f s (worst V %+.3f)\n', ...
            best.name, best.t, best.Vw);
    fprintf('   for comparison: 23 interpolated points, T = 1.90 -> 41.8 s\n');
    fprintf('                   today (20 points, T = 2.0)     -> 38.0 s but UNSAFE\n');
else
    fprintf('nothing with the original 20 points stayed inside the tubes\n');
end

f = figure('Position',[100 100 950 560]);
subplot(2,1,1); hold on;
for i = 1:numel(res), plot(res(i).tt, res(i).uu, 'LineWidth', 1.1); end
ylabel('commanded speed [ft/s]'); grid on;
legend(compose("T = %.1f s", [res.T]), 'Location','southeast', 'FontSize', 8);
title('original 20 trim points, segment time swept');
subplot(2,1,2); hold on;
for i = 1:numel(res), plot(res(i).tt, res(i).V, 'LineWidth', 1.1); end
yline(0,'k-','LineWidth',1.4); ylim([-0.5 0.4]);
xlabel('time [s]'); ylabel('BRT value (V<0 = inside)'); grid on;
saveas(f, fullfile(root,'logger','fixed_anchor_time.png'));
save(fullfile(root,'logger','fixed_anchor_time.mat'), 'res','T_need','C','T_LONG');
fprintf('\nsaved logger/fixed_anchor_time.{mat,png}\n');

% =========================================================================
function pk = ramp_peak(P, ve, T, T_hold, dt, U_ROW, NX)
d = zeros(NX,1);  pk = 0;
for i = 1:round(T/dt)
    d = P*d - ve*dt;
    if abs(d(U_ROW)) > pk, pk = abs(d(U_ROW)); end
end
for j = 1:round(T_hold/dt)
    d = P*d;
    if abs(d(U_ROW)) > pk, pk = abs(d(U_ROW)); end
end
end

function [t, u] = profile_from(us, Ts, dt, T_hold)
u = repmat(us(1), 1, round(T_hold/dt));
for k = 1:numel(Ts)
    n = max(round(Ts(k)/dt), 1);
    u = [u, us(k) + (us(k+1)-us(k))*(0:n-1)/n]; %#ok<AGROW>
end
u = [u, repmat(us(end), 1, round(T_hold/dt)+1)];
t = (0:numel(u)-1)*dt;
end

function R = fly(params, tt, uu, dt, brtV, T_hold)
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
N = numel(tt);  pos = cumtrapz(tt, uu);
V = zeros(1,N);  ua = zeros(1,N);
for k = 1:N
    ua(k) = guam.state(4);
    V(k)  = brtV.value(guam.state([4 6 11 8]), uu(k));
    guam.step(struct('pos',[pos(k);0;-80], 'vel',[uu(k);0;0], 'chi',0, 'chi_dot',0));
end
w = round(T_hold/dt)+1 : N-round(T_hold/dt);
R.V = V;  [R.Vw, iw] = max(V(w));  R.u_worst = uu(w(iw));
R.viol = 100*sum(V(w) >= 0)/numel(w);
R.err  = max(abs(ua(w)-uu(w)));
end
