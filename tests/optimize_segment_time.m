% optimize_segment_time - the fastest transition that still satisfies the
% safety condition, with the segment time T free instead of pinned at 2 s.
%
% ---------------------------------------------------------------------------
% WHY T HAS TO BE FREED
% ---------------------------------------------------------------------------
% place_trim_points.m held T = 2 s everywhere and moved only the step length.
% That spends the spare margin at high speed on LONGER STEPS and none of it on
% GOING FASTER, so the high-speed half of the transition still takes 2 s per
% segment even though the tubes there overlap two or three times more than
% they need to. Freeing T is what turns spare margin into saved time.
%
% ---------------------------------------------------------------------------
% THE PROBLEM
% ---------------------------------------------------------------------------
% Per segment there are now two design variables, the step D and the time T,
% tied by the same admissibility condition:
%
%   g = a_fwd(u) + a_bwd(u+D) - D - 2*( e0(u) + c(u,T)*D/T )  >=  0      (G)
%
% Total transition time is the sum of the segment times. Locally, the time
% spent per unit of speed gained is
%
%   cost(u,T) = ( T + T_dwell(u) ) / D_max(u,T)                          (C)
%
% where D_max solves g = 0. The objective is separable in u - each speed's
% cost depends only on the T chosen there - so minimising (C) pointwise
% minimises the total. No search over whole schedules is needed.
%
% ---------------------------------------------------------------------------
% WHY THE ANSWER IS NOT SIMPLY "T -> 0"
% ---------------------------------------------------------------------------
% Shrinking T makes the ramp steeper, which makes the lag worse, which forces
% a shorter step - those two cancel to first order and (C) keeps falling.
% Taken literally the model would say an infinitely fast transition is safe,
% which is nonsense. Two things stop it, and only one of them is in the model:
%
%   IN the model:  the dwell. If the vehicle has to settle at each trim point
%       before the next step is commanded, every segment pays T_dwell on top
%       of T, and chopping T into ever smaller pieces multiplies that toll.
%       T_dwell is computed here from the same linear model - it is the time
%       for the speed deviation to fall to DWELL_FRAC of its peak after the
%       ramp ends - so it is not a free parameter.
%
%   NOT in the model:  actuator saturation. Phi was measured with small
%       perturbations, so it describes a loop with unlimited control power. A
%       real ramp steep enough will saturate the rotors and the elevator, and
%       the true error will exceed c*s. That is why the fastest schedule found
%       here MUST be flown before it is believed - see verify_optimized_time.m.
%
% Output: logger/segment_time_opt.mat, logger/segment_time_opt.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt         = 0.01;
WH_IDX     = 3;
N_TRIM     = 20;
E0_WIND    = 0.0;
DWELL_FRAC = 0.10;    % "settled" = speed deviation down to 10 % of its peak
T_GRID     = [0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0 4.0 5.0 6.0];
T_HOLD     = 8.0;     % how long to watch after the ramp when finding the peak
U_ROW      = 4;

%% ---------------------------------------------------------------------
%  1. inputs: a(u) from the tubes, Phi/xi_e from the closed-loop model
%% ---------------------------------------------------------------------
spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4, gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d)); end
ic = [find(gv{1}==0), find(gv{2}==0), find(gv{3}==0), find(gv{4}==0)];

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj.trim_lon(1, :);

a_fwd = zeros(1,N_TRIM);  a_bwd = zeros(1,N_TRIM);
for k = 1:N_TRIM
    S = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    line = squeeze(S.data(:, ic(2), ic(3), ic(4)));
    a_fwd(k) = cross_zero(gv{1}, line, ic(1), +1);
    a_bwd(k) = cross_zero(gv{1}, line, ic(1), -1);
end

Lm = load(fullfile(root,'logger','closed_loop_model_ramped.mat'));  M = Lm.M;
Lc = load(fullfile(root,'logger','lag_coefficient.mat'));
dXE = Lc.dXE;                                   % position rows already zeroed
e0_u = abs(M.xi_e(U_ROW,:) - M.u_anchor) + E0_WIND;

%% ---------------------------------------------------------------------
%  2. c(u,T) and T_dwell(u,T) on a grid
%     Same propagation as compute_lag_coefficient, now swept over T.
%% ---------------------------------------------------------------------
nT = numel(T_GRID);
C  = zeros(nT, N_TRIM);      % lag coefficient [s]
D_dwell = zeros(nT, N_TRIM); % settle time after the ramp ends [s]
for a = 1:nT
    for k = 1:N_TRIM
        [C(a,k), D_dwell(a,k)] = ramp_peak(M.Phi{k}, dXE(:,k), T_GRID(a), ...
                                           T_HOLD, dt, U_ROW, ...
                                           ClosedLoopModel.N_XI, DWELL_FRAC);
    end
end

fprintf('--- c(u,T) [s] ---\n');
fprintf('  T   |');  fprintf(' u=%4.0f', u_anchor([1 4 8 12 16 20]));  fprintf('\n');
for a = 1:nT
    fprintf('%5.2f |', T_GRID(a));  fprintf(' %6.3f', C(a,[1 4 8 12 16 20]));  fprintf('\n');
end
fprintf('\n--- settle time after the ramp ends [s] ---\n');
fprintf('  T   |');  fprintf(' u=%4.0f', u_anchor([1 4 8 12 16 20]));  fprintf('\n');
for a = 1:nT
    fprintf('%5.2f |', T_GRID(a));  fprintf(' %6.2f', D_dwell(a,[1 4 8 12 16 20]));  fprintf('\n');
end

%% ---------------------------------------------------------------------
%  3. D_max(u,T) and the cost, with and without the dwell
%% ---------------------------------------------------------------------
A_f = @(u) interp1(u_anchor, a_fwd, clampu(u,u_anchor));
A_b = @(u) interp1(u_anchor, a_bwd, clampu(u,u_anchor));
E0  = @(u) interp1(u_anchor, e0_u,  clampu(u,u_anchor));
Cf  = @(u,a) interp1(u_anchor, C(a,:), clampu(u,u_anchor));
Df  = @(u,a) interp1(u_anchor, D_dwell(a,:), clampu(u,u_anchor));

% D_max(u, T): the longest step that still satisfies (G) at this speed and
% segment time. g falls monotonically in D, so bisection is safe.
dmax = @(u, a) solve_dmax(u, A_f, A_b, E0, Cf(u,a), T_GRID(a));

for USE_DWELL = [false true]
    if USE_DWELL, tag = 'WITH dwell (settle at every trim point)';
    else,         tag = 'NO dwell (ramp straight into the next segment)'; end
    fprintf('\n\n============ %s ============\n', tag);

    % ---- cost table: seconds spent per ft/s of speed gained --------------
    fprintf('\n--- cost (s per ft/s gained); lower is faster ---\n');
    fprintf('  T   |');  fprintf(' u=%4.0f', u_anchor([1 4 8 12 16 20]));  fprintf('\n');
    cost = zeros(nT, N_TRIM);
    for a = 1:nT
        for k = 1:N_TRIM
            D = dmax(u_anchor(k), a);
            if D <= 1e-4
                cost(a,k) = inf;
            else
                td = USE_DWELL * Df(u_anchor(k), a);
                cost(a,k) = (T_GRID(a) + td) / D;
            end
        end
        fprintf('%5.2f |', T_GRID(a));
        fprintf(' %6.3f', cost(a,[1 4 8 12 16 20]));  fprintf('\n');
    end
    [~, ia_best] = min(cost, [], 1);
    fprintf('\nbest T per anchor:\n');
    fprintf('  u   ');  fprintf('%6.0f', u_anchor);  fprintf('\n');
    fprintf('  T*  ');  fprintf('%6.2f', T_GRID(ia_best));  fprintf('\n');

    % ---- walk with the locally best T ------------------------------------
    u = 0;  us = 0;  ds = [];  Ts = [];  tds = [];
    while u < u_anchor(end) - 1e-6 && numel(us) < 400
        best = inf;  ab = 0;  Db = 0;
        for a = 1:nT
            D = dmax(u, a);
            if D <= 1e-4, continue; end
            td = USE_DWELL * Df(u, a);
            cst = (T_GRID(a) + td) / D;
            if cst < best, best = cst;  ab = a;  Db = D;  end
        end
        if ab == 0, error('no admissible (D,T) at u = %.2f', u); end
        Db = min(Db, u_anchor(end) - u);
        u = u + Db;  us(end+1) = u;  ds(end+1) = Db; %#ok<SAGROW>
        Ts(end+1) = T_GRID(ab);  tds(end+1) = USE_DWELL * Df(us(end-1), ab); %#ok<SAGROW>
    end

    fprintf('\n--- schedule ---\n');
    fprintf('  k |    u    |  step | T [s] | dwell | seg time\n');
    for k = 1:numel(ds)
        fprintf('%3d | %7.2f | %5.2f | %5.2f | %5.2f | %6.2f\n', ...
                k, us(k), ds(k), Ts(k), tds(k), Ts(k)+tds(k));
    end
    tot = sum(Ts) + sum(tds);
    fprintf('\ntrim points %d | segments %d | TOTAL %.2f s\n', ...
            numel(us), numel(ds), tot);

    % ---- reference points ------------------------------------------------
    for Tf = [2.0 1.0 0.5]
        a = find(abs(T_GRID - Tf) < 1e-9, 1);
        uu = 0;  ns = 0;  ttl = 0;  bad = false;
        while uu < u_anchor(end) - 1e-6 && ns < 400
            D = dmax(uu, a);
            if D <= 1e-4, bad = true; break; end
            td = USE_DWELL * Df(uu, a);
            uu = uu + min(D, u_anchor(end)-uu);  ns = ns+1;  ttl = ttl + Tf + td;
        end
        if bad
            fprintf('uniform T = %.2f s : infeasible at u = %.1f\n', Tf, uu);
        else
            fprintf('uniform T = %.2f s : %3d segments, %.2f s\n', Tf, ns, ttl);
        end
    end

    if ~USE_DWELL
        R.nodwell = struct('us',us,'ds',ds,'Ts',Ts,'tot',tot,'cost',cost);
    else
        R.dwell   = struct('us',us,'ds',ds,'Ts',Ts,'tds',tds,'tot',tot,'cost',cost);
    end
end

%% ---------------------------------------------------------------------
%  4. figure
%% ---------------------------------------------------------------------
f = figure('Position',[100 100 1000 720]);
subplot(3,1,1);
plot(R.nodwell.us(1:end-1), R.nodwell.ds, 'o-','LineWidth',1.3); hold on;
plot(R.dwell.us(1:end-1),   R.dwell.ds,   's-','LineWidth',1.3);
yline(8.44,'r--'); ylabel('step [ft/s]'); grid on;
legend('no dwell','with dwell','uniform 8.44','Location','northwest');
title('optimised placement');

subplot(3,1,2);
stairs(R.nodwell.us(1:end-1), R.nodwell.Ts, 'LineWidth',1.4); hold on;
stairs(R.dwell.us(1:end-1),   R.dwell.Ts,   'LineWidth',1.4);
yline(2.0,'r--'); ylabel('segment time T [s]'); grid on;
legend('no dwell','with dwell','fixed 2 s','Location','northeast');

subplot(3,1,3);
plot(R.nodwell.us, [0 cumsum(R.nodwell.Ts)], 'LineWidth',1.4); hold on;
plot(R.dwell.us,   [0 cumsum(R.dwell.Ts + R.dwell.tds)], 'LineWidth',1.4);
plot(linspace(0,160,23), linspace(0,44,23), 'r--','LineWidth',1.2);
xlabel('forward speed [ft/s]'); ylabel('elapsed time [s]'); grid on;
legend('no dwell','with dwell','T = 2 s fixed (44 s)','Location','northwest');
saveas(f, fullfile(root,'logger','segment_time_opt.png'));

save(fullfile(root,'logger','segment_time_opt.mat'), ...
     'R','T_GRID','C','D_dwell','u_anchor','a_fwd','a_bwd','e0_u');
fprintf('\nsaved logger/segment_time_opt.mat and .png\n');

% -------------------------------------------------------------------------
function [pk, t_settle] = ramp_peak(P, ve, T, T_hold, dt, U_ROW, NX, frac)
d = zeros(NX,1);  pk = 0;
for i = 1:round(T/dt)
    d = P*d - ve*dt;
    if abs(d(U_ROW)) > pk, pk = abs(d(U_ROW)); end
end
nh = round(T_hold/dt);  hist = zeros(1,nh);
for j = 1:nh
    d = P*d;  hist(j) = abs(d(U_ROW));
    if hist(j) > pk, pk = hist(j); end
end
idx = find(hist < frac*pk, 1);
if isempty(idx), t_settle = T_hold; else, t_settle = idx*dt; end
end

function D = solve_dmax(u, A_f, A_b, E0, c_u, T)
gg = @(D) A_f(u) + A_b(u+D) - D - 2*(E0(u) + c_u*D/T);
if gg(1e-6) <= 0, D = 0; return; end
lo = 1e-6;  hi = 40;
if gg(hi) > 0, D = hi; return; end
for it = 1:60
    mid = (lo+hi)/2;
    if gg(mid) > 0, lo = mid; else, hi = mid; end
end
D = lo;
end

function d = cross_zero(g, line, i0, sgn)
n = numel(g);  i = i0;
while true
    j = i + sgn;
    if j < 1 || j > n, d = abs(g(i)-g(i0)); return; end
    if line(j) >= 0
        f = -line(i)/(line(j)-line(i));
        d = abs(g(i) + f*(g(j)-g(i)) - g(i0));  return;
    end
    i = j;
end
end

function y = clampu(u, ua), y = min(max(u, ua(1)), ua(end)); end
