% check_schedule - judge a per-segment reference-time schedule, and fix it.
%
% ===========================================================================
% WHY THE OLD RULE CANNOT BE REUSED HERE
% ===========================================================================
% With interpolated trim points the rule was
%
%     thickness = a_fwd(u_k) + a_bwd(u_k+1) - step        >=  2*e
%
% and the step was the thing being solved for. With the trim points pinned to
% the 20 that actually have a tube, the step is stuck at 8.44 and the
% thickness is stuck with it - 1.51 ft/s at hover. Solving the same
% inequality for T instead demands T = 15 s on segment 1. Flying it says
% T = 2.3 s is fine. The rule is off by a factor of six.
%
% It is off because of what it asks for. It asks that the WHOLE spread of
% where the vehicle might be fit inside the 1.51 ft/s window where both tubes
% overlap. That is sufficient for safety but nowhere near necessary:
%
%   - the vehicle's speed rises continuously from 0 to 160, so it passes
%     through every handover window on the way whether we plan it or not;
%   - the tubes are 4-D. Being off in speed is survivable if attitude, pitch
%     rate and vertical speed are comfortable, and the 1-D rule cannot see
%     that trade because it has thrown away three of the four dimensions.
%
% So the 1-D thickness rule is a screening test, not a verdict.
%
% ===========================================================================
% THE RULE USED HERE INSTEAD
% ===========================================================================
% Ask the question the design is actually about, and ask it in all four
% dimensions: IS THE VEHICLE INSIDE A REAL TUBE, AT EVERY INSTANT?
%
%   1. build the commanded speed profile from the segment times
%   2. propagate the gap equation to get the deviation d(t)   <- analytic,
%          d(i+1) = Phi(u) d(i) - xi_e'(u) * s(i) * dt           no flying
%   3. the predicted state is  x(t) = trim(u_ref(t)) + d(t)
%   4. look x(t) up in the tubes of the two neighbouring ANCHORS - the stored
%      arrays, not a blend of them - and take the better of the two
%   5. admissible if that value stays below zero for the whole transition
%
% Step 4 is the honest fixed-anchor test: the vehicle has to be inside a tube
% that was actually computed, not inside a guess made between two of them.
%
% A NOTE ON WHAT IS AND IS NOT INTERPOLATED
%   The reference passes continuously between trim points - it has to, the
%   speed cannot teleport - and reading the trim state at an in-between speed
%   is what RSLQR already does on every single timestep (interp_xu0). That is
%   not the thing the team objected to. The objection was to ANCHORING a
%   design point, and its reachable set, at a speed where no set was ever
%   computed. Step 4 keeps that promise: every tube used here is a stored one.
%
% ===========================================================================
% THE CORRECTION RULE - AND THE OBVIOUS VERSION OF IT THAT DOES NOT WORK
% ===========================================================================
% The tempting rule is "whichever segment fails, lengthen that segment". It
% gives nonsense here, and the reason is worth stating because it is a
% property of the problem, not of the code.
%
% The deviation d does not reset at a segment boundary. It carries over. So a
% SHORT segment can look innocent on its own scorecard - the ramp is steep but
% it is over quickly, and d has not finished growing by the time the segment
% ends - while handing the next segment a vehicle that is already far from
% trim. Bisecting each segment against its own worst value therefore drives
% several of them to the floor (0.5 s) and pushes the damage downstream; run
% that loop and it returns a 52 s schedule that is WORSE than the 44 s one it
% started from.
%
% Segments cannot be judged, or fixed, independently.
%
% What is done instead:
%   stage 1  scale the whole schedule up until the WHOLE-schedule value passes
%            - one monotone knob, one bisection, guaranteed to converge
%   stage 2  walk the segments and try to shorten each one, keeping a change
%            only if the WHOLE schedule still passes. Repeat until nothing
%            more can be shaved.
% Every acceptance test is global, so nothing can be bought by pushing cost
% into a neighbour.
%
% Output: logger/schedule_check.mat / .png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt      = 0.01;
WH_IDX  = 3;
N_TRIM  = 20;
N_SEG   = N_TRIM - 1;
T_HOLD  = 5.0;
V_LIMIT = -0.02;      % require this much margin, not merely V < 0

% Floor on the segment time. This is NOT a physical limit - it is the edge of
% where the model may be believed. c(u,T) was checked against a flown mission
% whose segments were 2 s long (verify_lag_coefficient.m) and the verdicts
% here were validated against flights at 2.0-4.0 s. Below about 1.5 s the ramp
% is steep enough that the small-perturbation linearisation behind Phi stops
% describing it, and the predictor starts under-estimating the deviation.
% Measured: without this floor, stage 2 shaves segments to 0.60 s, the
% predictor says the schedule passes at 40.6 s, and flying it puts the vehicle
% outside the tubes 6.6 % of the time (V = +0.316). Extrapolation, not physics.
T_FLOOR = 1.50;

%% ---------------------------------------------------------------------
%  inputs
%% ---------------------------------------------------------------------
Lm = load(fullfile(root,'logger','closed_loop_model_ramped.mat'));  M = Lm.M;
Lc = load(fullfile(root,'logger','lag_coefficient.mat'));           dXE = Lc.dXE;

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj0 = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj0.trim_lon(1,:);
trim_lon = traj0.trim_lon;                       % 4 x 20, the real trims

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4, gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d)); end
TUBE = cell(1,N_TRIM);
for k = 1:N_TRIM
    S = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    TUBE{k} = S.data;
end

%% ---------------------------------------------------------------------
%  1. JUDGE - validate the predictor against schedules already flown
%% ---------------------------------------------------------------------
F = load(fullfile(root,'logger','fixed_anchor_time.mat'));
fprintf('=== JUDGE: predicted vs flown, on schedules already flown ===\n');
fprintf('%7s | %9s | %9s | %8s | %s\n', ...
        'T [s]', 'pred Vw', 'flown Vw', 'diff', 'pred / flown verdict');
fprintf('%s\n', repmat('-',1,62));
for i = 1:numel(F.res)
    Tv = repmat(F.res(i).T, 1, N_SEG);
    P  = predict(Tv, u_anchor, trim_lon, M, dXE, TUBE, gv, dt, T_HOLD);
    vp = ternary(P.Vw < 0, 'SAFE  ', 'unsafe');
    vf = ternary(F.res(i).Vw < 0, 'SAFE  ', 'unsafe');
    fprintf('%7.2f | %+9.3f | %+9.3f | %+8.3f | %s / %s%s\n', ...
            F.res(i).T, P.Vw, F.res(i).Vw, P.Vw-F.res(i).Vw, vp, vf, ...
            ternary(sign(P.Vw)==sign(F.res(i).Vw), '', '   <-- DISAGREE'));
end

%% ---------------------------------------------------------------------
%  2. JUDGE a specific schedule, segment by segment
%% ---------------------------------------------------------------------
T_TRY = [repmat(3.0,1,4) repmat(2.5,1,4) repmat(2.0,1,11)];   % the current proposal
P = predict(T_TRY, u_anchor, trim_lon, M, dXE, TUBE, gv, dt, T_HOLD);
fprintf('\n=== JUDGE: the proposed schedule, segment by segment ===\n');
fprintf('seg |    u range     |  T   | worst V | margin | verdict\n');
fprintf('%s\n', repmat('-',1,62));
for k = 1:N_SEG
    fprintf('%3d | %5.1f -> %5.1f | %4.2f | %+7.3f | %+6.3f | %s\n', ...
            k, u_anchor(k), u_anchor(k+1), T_TRY(k), P.Vseg(k), ...
            -P.Vseg(k), ternary(P.Vseg(k) < V_LIMIT, 'ok', 'TOO SHORT'));
end
fprintf('total %.1f s | worst V %+.3f | %s\n', sum(T_TRY), P.Vw, ...
        ternary(P.Vw < V_LIMIT, 'ADMISSIBLE', 'NOT ADMISSIBLE'));

%% ---------------------------------------------------------------------
%  3. CORRECT - shortest admissible T per segment
%% ---------------------------------------------------------------------
ok_fn = @(Tv) getfield(predict(Tv, u_anchor, trim_lon, M, dXE, TUBE, gv, dt, T_HOLD), 'Vw') < V_LIMIT; %#ok<GFLD>

% ---- stage 1: scale the whole schedule until it passes -----------------
fprintf('\n=== CORRECT stage 1: scale the whole schedule ===\n');
lo = 1.0;  hi = 1.0;
while ~ok_fn(T_TRY*hi) && hi < 8
    hi = hi*1.25;
end
if hi >= 8
    fprintf('even 8x the proposal does not pass - time alone cannot fix this\n');
    T_fix = T_TRY*8;
else
    if hi > 1.0
        lo = hi/1.25;
        for it = 1:14
            mid = (lo+hi)/2;
            if ok_fn(T_TRY*mid), hi = mid; else, lo = mid; end
        end
    end
    T_fix = T_TRY*hi;
    fprintf('scale %.3f -> total %.1f s\n', hi, sum(T_fix));
end

% ---- stage 2: shave whatever the WHOLE schedule still allows ------------
fprintf('\n=== CORRECT stage 2: shave, checking the whole schedule ===\n');
for pass = 1:4
    improved = false;
    for k = N_SEG:-1:1                       % high speed first, most margin
        step = 0.25;
        while T_fix(k) - step >= T_FLOOR
            Tt = T_fix;  Tt(k) = Tt(k) - step;
            if ok_fn(Tt), T_fix = Tt;  improved = true;  else, break; end
        end
    end
    P = predict(T_fix, u_anchor, trim_lon, M, dXE, TUBE, gv, dt, T_HOLD);
    fprintf('pass %d -> total %.1f s | worst V %+.3f\n', pass, sum(T_fix), P.Vw);
    if ~improved, break; end
end

P = predict(T_fix, u_anchor, trim_lon, M, dXE, TUBE, gv, dt, T_HOLD);
fprintf('\nseg |    u range     | T_fix | vs proposal\n');
fprintf('%s\n', repmat('-',1,46));
for k = 1:N_SEG
    fprintf('%3d | %5.1f -> %5.1f | %5.2f | %+6.2f s\n', ...
            k, u_anchor(k), u_anchor(k+1), T_fix(k), T_fix(k)-T_TRY(k));
end
fprintf('\ncorrected schedule: total %.1f s | worst V %+.3f\n', sum(T_fix), P.Vw);
fprintf('proposal was      : total %.1f s\n', sum(T_TRY));

fprintf('\nT_fix = [');  fprintf(' %.2f', T_fix);  fprintf(' ]\n');

%% ---------------------------------------------------------------------
%  4. VERIFY - and the honest state of this tool
%
%  WHERE THE JUDGE IS TRUSTWORTHY
%     Uniform schedules. Checked against eight flown ones from 2.0 to 4.0 s:
%     the verdict was right every time, including the exact boundary (2.2 s
%     unsafe, 2.3 s safe). Use it freely to pick a uniform segment time.
%
%  WHERE IT IS NOT
%     Non-uniform schedules. It is OPTIMISTIC there, and two attempts to patch
%     that failed:
%       - floor the segment time at 1.5 s (the edge of the verified band):
%         the tool returns 40.1 s, flying it gives 3.4 % outside, V = +0.093
%       - demand a wider margin instead (V < -0.15): the predicted value
%         saturates around -0.08 for long segments, so nothing satisfies it
%         until the schedule reaches 352 s. Useless in the other direction.
%     The judge was only ever validated where the segment time is the same
%     everywhere, and mixing segment times is a regime it has not been shown
%     to describe. Anything non-uniform has to be flown.
%
%  SO, TODAY
%     Use the judge to choose a UNIFORM T. Anything non-uniform: fly it.
%     Flight-verified schedules right now: uniform T = 2.30 -> 43.7 s, and
%     3.0/2.5/2.0 by speed band -> 44.0 s with four times the margin.
%
%  TO CLOSE THE GAP
%     Fly three or four non-uniform schedules and compare predicted against
%     flown, exactly as was done for the uniform ones. That is what turns the
%     correction stage into something usable without flying.
%% ---------------------------------------------------------------------
brtV = BRTValue(fullfile(root,'data'), trim_lon, WH_IDX);
fprintf('\n=== VERIFY: fly the corrected schedule ===\n');
fprintf('%-22s | %6s | %7s | %8s | %s\n','schedule','time','V>=0 %','flown Vw','verdict');
fprintf('%s\n', repmat('-',1,62));
for c2 = {{'proposal (44 s)', T_TRY}, {'corrected (40.6 s)', T_fix}}
    nm2 = c2{1}{1};  Tv = c2{1}{2};
    R = fly_sched(Tv, u_anchor, params, dt, T_HOLD, brtV);
    fprintf('%-22s | %6.1f | %7.1f | %+8.3f | %s\n', nm2, sum(Tv), R.viol, R.Vw, ...
            ternary(R.viol==0,'SAFE','unsafe'));
end

f = figure('Position',[100 100 950 560]);
subplot(2,1,1);
stairs([u_anchor(1:end-1) u_anchor(end)], [T_fix T_fix(end)], 'LineWidth', 1.8); hold on;
stairs([u_anchor(1:end-1) u_anchor(end)], [T_TRY T_TRY(end)], '--', 'LineWidth', 1.4);
ylabel('segment time T [s]'); grid on;
legend('corrected', 'proposal', 'Location','northeast');
title(sprintf('corrected %.1f s vs proposal %.1f s', sum(T_fix), sum(T_TRY)));
subplot(2,1,2);
plot(P.t, P.V, 'LineWidth', 1.4); hold on;
yline(0,'k-','LineWidth',1.4); yline(V_LIMIT,'r--');
xlabel('time [s]'); ylabel('predicted tube value'); grid on;
saveas(f, fullfile(root,'logger','schedule_check.png'));
save(fullfile(root,'logger','schedule_check.mat'), 'T_TRY','T_fix','P','V_LIMIT');
fprintf('\nsaved logger/schedule_check.{mat,png}\n');

% =========================================================================
function P = predict(Tv, ua, trim_lon, M, dXE, TUBE, gv, dt, T_hold)
% Predict the tube value along a schedule, without flying.
N_SEG = numel(Tv);
% --- commanded speed profile -------------------------------------------
u = repmat(ua(1), 1, round(T_hold/dt));
seg = zeros(size(u));
for k = 1:N_SEG
    n = max(round(Tv(k)/dt),1);
    u = [u, ua(k) + (ua(k+1)-ua(k))*(0:n-1)/n]; %#ok<AGROW>
    seg = [seg, repmat(k,1,n)]; %#ok<AGROW>
end
u   = [u, repmat(ua(end),1,round(T_hold/dt)+1)];
seg = [seg, repmat(N_SEG,1,round(T_hold/dt)+1)];
t   = (0:numel(u)-1)*dt;
s   = [diff(u)/dt, 0];

% --- gap equation -------------------------------------------------------
NX = ClosedLoopModel.N_XI;
d  = zeros(NX,1);
sel = [4 6 11 8];                              % u, w, q, theta inside xi
X  = zeros(4, numel(u));
for i = 1:numel(u)
    xe = interp_cols(trim_lon, ua, u(i));      % reference trim state
    X(:,i) = xe + d(sel);
    Pm = M.at(u(i));
    ve = interp_cols(dXE, M.u_anchor, u(i));
    d  = Pm*d - ve*s(i)*dt;
end

% --- value in the two neighbouring ANCHOR tubes (no blending) -----------
V = zeros(1, numel(u));
for i = 1:numel(u)
    k = max(min(find(ua <= u(i), 1, 'last'), numel(ua)-1), 1);
    V(i) = min(tube_val(TUBE{k},   X(:,i) - trim_lon(:,k),   gv), ...
               tube_val(TUBE{k+1}, X(:,i) - trim_lon(:,k+1), gv));
end

w = round(T_hold/dt)+1 : numel(u)-round(T_hold/dt);
P.t = t;  P.u = u;  P.V = V;  P.X = X;
P.Vw = max(V(w));
P.Vseg = zeros(1, N_SEG);
for k = 1:N_SEG
    m = (seg == k);  m(1:round(T_hold/dt)) = false;
    if any(m), P.Vseg(k) = max(V(m)); else, P.Vseg(k) = -inf; end
end
end

function v = tube_val(D, xc, gv)
% 4-D linear interpolation of one stored tube at relative state xc
i0 = zeros(1,4);  f = zeros(1,4);
for d = 1:4
    g = gv{d};  n = numel(g);  h = g(2)-g(1);
    tq = (xc(d) - g(1))/h + 1;
    tq = min(max(tq,1), n);
    b  = min(floor(tq), n-1);
    i0(d) = b;  f(d) = tq - b;
end
c = D(i0(1):i0(1)+1, i0(2):i0(2)+1, i0(3):i0(3)+1, i0(4):i0(4)+1);
c = c(:,:,:,1)*(1-f(4)) + c(:,:,:,2)*f(4);
c = c(:,:,1)  *(1-f(3)) + c(:,:,2)  *f(3);
c = c(:,1)    *(1-f(2)) + c(:,2)    *f(2);
v = c(1)      *(1-f(1)) + c(2)      *f(1);
end

function col = interp_cols(A, x, xq)
if xq <= x(1),   col = A(:,1);   return; end
if xq >= x(end), col = A(:,end); return; end
k = find(x <= xq, 1, 'last');  k = min(k, numel(x)-1);
a = (xq - x(k))/(x(k+1)-x(k));
col = (1-a)*A(:,k) + a*A(:,k+1);
end

function R = fly_sched(Tv, ua, params, dt, T_hold, brtV)
u = repmat(ua(1),1,round(T_hold/dt));
for k = 1:numel(Tv)
    n = max(round(Tv(k)/dt),1);
    u = [u, ua(k) + (ua(k+1)-ua(k))*(0:n-1)/n]; %#ok<AGROW>
end
u = [u, repmat(ua(end),1,round(T_hold/dt)+1)];
t = (0:numel(u)-1)*dt;  pos = cumtrapz(t,u);
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
V = zeros(1,numel(u));
for k = 1:numel(u)
    V(k) = brtV.value(guam.state([4 6 11 8]), u(k));
    guam.step(struct('pos',[pos(k);0;-80],'vel',[u(k);0;0],'chi',0,'chi_dot',0));
end
w = round(T_hold/dt)+1 : numel(u)-round(T_hold/dt);
R.Vw = max(V(w));  R.viol = 100*sum(V(w)>=0)/numel(w);
end

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
