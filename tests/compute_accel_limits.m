% compute_accel_limits - the fastest safe acceleration for each transition
% segment, and the transition time that follows from it.
%
% ---------------------------------------------------------------------------
% THE QUESTION
% ---------------------------------------------------------------------------
% The transition is a chain of 20 trim points. Handing over from trim k to
% trim k+1 is only safe if, at some moment, the aircraft is inside BOTH
% reachable sets at once. The aircraft is never on the trim point while it is
% accelerating - it lags, and it must pitch over to accelerate at all - so the
% question is how hard it may accelerate before that offset pushes it out of
% the overlap.
%
% Each segment gets its own constant acceleration. The answer is 19 numbers.
%
% ---------------------------------------------------------------------------
% THE ERROR IS CARRIED, NOT RESET
% ---------------------------------------------------------------------------
% An earlier version measured the offset segment by segment, each time
% starting from a perfectly trimmed aircraft. That is wrong, and it is wrong
% in the unsafe direction at high speed. Phi carries a mode with |lambda| =
% 0.9991, a time constant of 12 s and a period of 59 s, living in the position
% error and the longitudinal integrator - RSLQR's deliberately slow outer
% loop. Over a 33 s transition that mode accumulates, and resetting it at
% every anchor throws the accumulation away: the flown w error at the last
% three segments came out 26-32 % ABOVE what the per-segment measurement
% predicted, enough to leave the tube.
%
% So the error is marched continuously across the whole transition. No
% iteration is needed for this, because the error entering segment k depends
% only on segments 1..k-1, which are already decided. One forward pass:
%
%     segment 1 : start from zero error, solve for a_1, march, keep the error
%     segment 2 : start from THAT error, solve for a_2, march, keep it
%     ...
%
% The acceleration is chosen greedily - the largest that keeps segment k safe,
% given what the earlier segments left behind. That is not guaranteed globally
% optimal (going slower early leaves less error later), and the gap is
% reported at the end so the cost of greediness is visible.
%
% ---------------------------------------------------------------------------
% WHERE THE ERROR COMES FROM
% ---------------------------------------------------------------------------
% ClosedLoopModel supplies, for a command held at v,
%
%       xi(i+1) - xi_e(v) = Phi(v) ( xi(i) - xi_e(v) )
%
% over the whole 34-state closed loop (position error, body velocities,
% attitude, rates, 9 propellers, 5 surfaces, both integrators, both virtual
% attitude effectors). Let the command RAMP at a ft/s^2: the equilibrium the
% gap is measured from also moves, by xi_e'(v)*a*dt per step, so with
% d = xi - xi_e
%
%       d(i+1) = Phi d(i) - xi_e'(v) a dt
%
% Nothing is simulated - this is 34x34 matrix products. The aircraft's actual
% longitudinal state at commanded speed v is then
%
%       x_act(v) = x_trim(v) + d(v)|_[u w q theta]
%
% and the two sets are read at x_act minus each anchor's own trim state. The
% forward-speed lag needs no special handling: it is inside x_act.
%
% ---------------------------------------------------------------------------
% TWO DIFFERENT w's, AND WHICH IS WHICH
% ---------------------------------------------------------------------------
%   BRT w axis - HEADING FRAME. Confirmed by the author. The trim table holds
%       that quantity at 11.667 ft/s for every anchor, so an aircraft sitting
%       on a trim point has zero deviation on this axis.
%
%   Aero-model validity - BODY FRAME. Aerodynamic forces depend on the flow
%       over the airframe, so check_fac_limits in LpC_interp_p_v2 bounds the
%       BODY-frame w at +-10 kt. That is u*sin(theta) + w_h*cos(theta), which
%       reaches 14.75 kt at the cruise trim points - outside the range the
%       model was identified over, at zero acceleration. Pitching the nose
%       down shrinks the u*sin(theta) term, so acceleration brings it back in
%       and there is a MINIMUM acceleration as well as a maximum.
%
%       Separate from both: the reachable sets at anchors 12-20 were themselves
%       computed with the aero model outside its identified range. That is a
%       defect in the data, not something an acceleration choice can repair.
%       It is reported, not corrected - these are the only sets that exist.
%
% ---------------------------------------------------------------------------
% WHY BISECTION AND NOT A FORMULA
% ---------------------------------------------------------------------------
% V is a numerical Hamilton-Jacobi solution stored on a grid, so there is no
% closed form and any method must query the table. Safety is monotone in a
% (more acceleration, more offset, less overlap), so bisection converges with
% a guaranteed bracket. A 2x2 Newton solve on the tangency condition agrees to
% 0.003 but needs more iterations, because the stored value function is
% piecewise-linear and its gradient jumps at every cell boundary.
%
% a_lo and a_hi, by contrast, ARE closed form. With theta(a) = theta_trim -
% a*g_theta, collapse w_body into one sine, w_body = R sin(theta+phi) with
% R = hypot(u, w_h) and phi = atan2(w_h, u); then |w_body| <= W_lim inverts to
%
%       a_lo = ( theta_trim - asin(W_lim/R) + phi ) / g_theta
%       a_hi = ( theta_trim + asin(W_lim/R) + phi ) / g_theta
%
% Output: logger/accel_limits.mat, logger/accel_limits.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 0. constants --------------------------------------------------------
WH_IDX   = 3;          % wbar = +11.667 ft/s, the BRT family that exists
K        = 20;         % trim points
NS       = K - 1;      % segments
KT       = 1.6878;     % ft/s per knot
W_LIM_KT = 10;         % aero model's |v|,|w| range, every branch
W_LIM    = W_LIM_KT * KT;

A_MAX    = 12;         % bisection bracket, ft/s^2
A_MIN    = 0.15;       % below this a segment is treated as infeasible
N_BISECT = 22;

SEL = [4 6 11 8];      % u, w, q, theta inside the 34-state xi
LON = [1 3 5 11];      % u, w, q, theta inside the 25-state XU0

line = @(c) fprintf('%s\n', repmat(c, 1, 94));

%% 1. trim points ------------------------------------------------------
T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:K);   UH = UH(:)';
WH3 = T.WH(WH_IDX);
XU0 = T.XU0_interp(:, 1:K, WH_IDX);
trim_lon = XU0(LON, :);               % 4 x K, [u; w; q; theta]
th_trim  = XU0(11, :);
du_seg   = diff(UH);

line('='); fprintf('TRIM POINTS\n'); line('=');
fprintf('  WH index %d -> wbar = %+.4f ft/s (%.2f kt), constant over all anchors\n', ...
        WH_IDX, WH3, WH3/KT);
fprintf('  %d anchors, %.2f to %.2f ft/s, step %.2f ft/s\n\n', ...
        K, UH(1), UH(end), du_seg(1));

%% 2. closed-loop model, with health checks ----------------------------
Sm = load(fullfile(root, 'logger', 'closed_loop_model_ramped.mat'));  M = Sm.M;
Sl = load(fullfile(root, 'logger', 'lag_coefficient.mat'));
dt  = Sl.dt;   dXE = Sl.dXE;

assert(numel(Sl.u_anchor) == K && max(abs(Sl.u_anchor(:)' - UH)) < 1e-6, ...
    'compute_accel_limits:anchorMismatch', ...
    'The stored closed-loop model was built at different anchors.');

line('='); fprintf('CLOSED-LOOP MODEL HEALTH\n'); line('=');

% (a) the recursion must be exactly linear in the acceleration
p1 = ramp_peak(M.Phi{3}, dXE(:,3), dt, SEL, 2, 6, 1);
p3 = ramp_peak(M.Phi{3}, dXE(:,3), dt, SEL, 2, 6, 3);
fprintf('  linear in a          : |3*d(a=1) - d(a=3)| = %.2e  (want ~1e-14)\n', ...
        max(abs(3*p1 - p3)));
assert(max(abs(3*p1 - p3)) < 1e-9, 'compute_accel_limits:notLinear', ...
    'The ramp recursion is not linear in the acceleration.');

% (b) commanding 1 ft/s more must settle 1 ft/s faster
fprintf('  du_e/dv over anchors : %.4f .. %.4f  (want ~1)\n', ...
        min(dXE(4,:)), max(dXE(4,:)));

% (c) the settled equilibrium must agree with the trim table it is meant to be
dev = M.xi_e(SEL, :) - trim_lon;
fprintf('  |xi_e - trim table|  : u %.2f ft/s, w %.2f ft/s, theta %.2f deg (max)\n', ...
        max(abs(dev(1,:))), max(abs(dev(2,:))), rad2deg(max(abs(dev(4,:)))));

% (d) the slow mode that forces the error to be carried
lamx = zeros(1,K); taux = zeros(1,K);
for k = 1:K
    lam = eig(M.Phi{k});
    tau = -dt ./ log(abs(lam));
    ok  = abs(lam) < 0.99995 & abs(imag(lam)) > 1e-6;
    [taux(k), j] = max(tau .* ok);   lamx(k) = abs(lam(j));
end
fprintf('  slowest osc. mode    : |lambda| %.4f, tau %.1f s (this is why the\n', ...
        max(lamx), max(taux));
fprintf('                         error is marched, not reset per segment)\n\n');

%% 3. BRT value functions ----------------------------------------------
spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
V = cell(1, K);
for k = 1:K
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    V{k} = griddedInterpolant(gv, S.data, 'linear', 'nearest');
end

%% 4. one forward pass: solve each segment carrying the error in --------
line('='); fprintf('FORWARD PASS  -  error carried across segments\n'); line('=');
fprintf('   k   u_trim    a [ft/s^2]   T_k [s]   handover u   error IN (w, th)   error OUT\n');

a_use = zeros(1, NS);   T_k = zeros(1, NS);
u_hand = nan(1, NS);    feasible = true(1, NS);
d_in  = zeros(34, NS);  d_out = zeros(34, NS);
pk_seg = zeros(4, NS);  wb_seg = zeros(2, NS);

d = zeros(34, 1);                    % trimmed at hover before the transition
for k = 1:NS
    d_in(:, k) = d;

    % Largest a that keeps the aircraft covered the whole way across the
    % segment. The bracket starts at A_MIN, never at zero: a vanishing a makes
    % the segment take forever and the march step count explode.
    if march(V, M, dXE, trim_lon, UH, k, A_MAX, d, dt, SEL, W_LIM)
        lo = A_MAX;                              % safe at the top of the bracket
    elseif ~march(V, M, dXE, trim_lon, UH, k, A_MIN, d, dt, SEL, W_LIM)
        lo = 0;                                  % unsafe even at the bottom
    else
        lo = A_MIN;  hi = A_MAX;
        for i = 1:N_BISECT
            m = (lo + hi)/2;
            if march(V, M, dXE, trim_lon, UH, k, m, d, dt, SEL, W_LIM)
                lo = m;
            else
                hi = m;
            end
        end
    end

    if lo < A_MIN
        feasible(k) = false;  a_use(k) = NaN;  T_k(k) = NaN;
        fprintf('%4d %8.1f %12s %9s %12s   *** no safe acceleration ***\n', ...
                k, UH(k), 'none', '-', '-');
        break;
    end

    a_use(k) = lo;
    T_k(k)   = du_seg(k) / lo;
    [~, d, u_hand(k), pk_seg(:,k), wb_seg(:,k)] = ...
        march(V, M, dXE, trim_lon, UH, k, lo, d, dt, SEL, W_LIM);
    d_out(:, k) = d;

    fprintf('%4d %8.1f %12.2f %9.2f %12.1f %8.2f %6.2f %8.2f %6.2f\n', ...
            k, UH(k), a_use(k), T_k(k), u_hand(k), ...
            d_in(6,k), rad2deg(d_in(8,k)), pk_seg(2,k), rad2deg(pk_seg(4,k)));
end

T_min = sum(T_k(feasible));

%% 5. aero-model validity, as a separate diagnostic --------------------
line('='); fprintf('AERO-MODEL VALIDITY  -  body-frame w against +-%.0f kt\n', W_LIM_KT); line('=');
fprintf('  the reachable sets at the OUTSIDE anchors were themselves built with\n');
fprintf('  the model out of range. That is a data defect, reported not repaired.\n\n');

% g_theta only for the closed-form window; the flown value comes from the march
g_th = zeros(1, K);
for k = 1:K
    pk = ramp_peak(M.Phi{k}, dXE(:,k), dt, SEL, 3, 8, 1);
    g_th(k) = pk(4);
end
[alo_pt, ahi_pt, wb0] = deal(zeros(1, K));
for k = 1:K
    [alo_pt(k), ahi_pt(k)] = validity_window(UH(k), WH3, th_trim(k), g_th(k), W_LIM);
    wb0(k) = UH(k)*sin(th_trim(k)) + WH3*cos(th_trim(k));
end
a_lo = max(alo_pt(1:NS), alo_pt(2:K));
a_hi = min(ahi_pt(1:NS), ahi_pt(2:K));

fprintf('   k   u_trim   w_body at trim [kt]   sets built in range?     a_lo     a_hi\n');
for k = 1:K
    if abs(wb0(k)) > W_LIM, st = 'NO  <-- defect'; else, st = 'yes           '; end
    fprintf('%4d %8.1f %17.2f   %s %10.2f %8.2f\n', ...
            k, UH(k), wb0(k)/KT, st, alo_pt(k), ahi_pt(k));
end

fprintf('\n   k   u_trim   a used    a_lo    a_hi   |w_body| flown [kt]   verdict\n');
for k = 1:NS
    if ~feasible(k), continue; end
    wmax = max(abs(wb_seg(:,k)))/KT;
    if wmax > W_LIM_KT, vd = 'OUT OF MODEL RANGE'; else, vd = 'in range'; end
    fprintf('%4d %8.1f %8.2f %7.2f %7.2f %18.2f   %s\n', ...
            k, UH(k), a_use(k), a_lo(k), a_hi(k), wmax, vd);
end

%% 6. what greediness cost --------------------------------------------
% Re-run the whole pass with every acceleration scaled by a common factor.
% If the greedy schedule were leaving time on the table, some factor below 1
% would buy enough late-segment speed to beat it.
line('='); fprintf('IS THE GREEDY CHOICE COSTING ANYTHING?\n'); line('=');
fprintf('  scale    total T [s]   safe?\n');
best_scale = 1;  best_T = T_min;
for sc = [0.70 0.80 0.90 0.95 1.00]
    [Ts, okall] = replay(V, M, dXE, trim_lon, UH, du_seg, a_use*sc, dt, SEL, W_LIM);
    if okall, tag = 'yes'; else, tag = 'NO'; end
    fprintf('%7.2f %13.2f   %s\n', sc, Ts, tag);
    if okall && Ts < best_T, best_T = Ts;  best_scale = sc; end
end
if best_scale < 1
    fprintf('\n  a uniformly slower schedule is FASTER overall (scale %.2f, %.2f s).\n', ...
            best_scale, best_T);
    fprintf('  greedy is leaving %.2f s on the table.\n', T_min - best_T);
else
    fprintf('\n  no uniform slow-down beats greedy.\n');
end

%% 7. result -----------------------------------------------------------
line('='); fprintf('RESULT\n'); line('=');
fprintf('  segments with no safe acceleration : %d\n', sum(~feasible));
fprintf('  transition time                    : %.2f s\n', T_min);
[amin, kmin] = min(a_use(feasible));
fprintf('  slowest segment                    : k = %d, a = %.2f ft/s^2, %.2f s\n', ...
        kmin, amin, T_k(kmin));
fprintf('  peak error over the transition     : w %.2f ft/s, theta %.2f deg\n', ...
        max(pk_seg(2,feasible)), rad2deg(max(pk_seg(4,feasible))));

%% 8. figure -----------------------------------------------------------
f = figure('Position', [80 80 1180 720], 'Color', 'w');

subplot(2,2,1);
stairs(UH(1:NS), a_use, 'LineWidth', 2.2); hold on;
plot(UH(1:NS), a_lo, '--', 'LineWidth', 1.4);
plot(UH(1:NS), a_hi, '--', 'LineWidth', 1.4);
grid on; xlabel('u [ft/s]'); ylabel('a [ft/s^2]');
title('acceleration per segment');
legend({'a used','a_{lo} validity','a_{hi} validity'}, 'Location','northwest');

subplot(2,2,2);
bar(UH(1:NS), T_k, 'FaceColor', [.35 .55 .62]);
grid on; xlabel('u [ft/s]'); ylabel('T_k [s]');
title(sprintf('segment time,  total %.1f s', T_min));

subplot(2,2,3);
plot(UH(1:NS), pk_seg(2,:), '-o', 'LineWidth', 1.8); hold on;
plot(UH(1:NS), d_in(6,1:NS), '-s', 'LineWidth', 1.4);
grid on; xlabel('u [ft/s]'); ylabel('w error [ft/s]');
title('w error: carried in vs peak within segment');
legend({'peak in segment','carried in'}, 'Location','best');

subplot(2,2,4);
plot(UH, wb0/KT, '-o', 'LineWidth', 1.8); hold on;
plot(UH(1:NS), max(abs(wb_seg))/KT, '-^', 'LineWidth', 1.6);
yline( W_LIM_KT, 'r--', 'LineWidth', 1.4);
grid on; xlabel('u [ft/s]'); ylabel('|w_{body}| [kt]');
title('body-frame w vs the aero model range');
legend({'at trim','flown'}, 'Location','best');

saveas(f, fullfile(root, 'logger', 'accel_limits.png'));

save(fullfile(root, 'logger', 'accel_limits.mat'), ...
     'UH','WH3','trim_lon','th_trim','du_seg','a_use','T_k','T_min','u_hand', ...
     'd_in','d_out','pk_seg','wb_seg','a_lo','a_hi','wb0','feasible','g_th', ...
     'W_LIM','WH_IDX','best_scale','best_T');
fprintf('\nsaved logger/accel_limits.mat and logger/accel_limits.png\n');

%% ---------------------------------------------------------------------
%  helpers
%% ---------------------------------------------------------------------
function pk = ramp_peak(P, dx, dt, SEL, Tramp, Thold, a)
% Peak |d| over a ramp of slope a held for Tramp, then Thold frozen. Used only
% for the health check and for g_theta in the closed-form validity window.
pk = zeros(4,1);  d = zeros(34,1);
for i = 1:round(Tramp/dt), d = P*d - dx*a*dt;  pk = max(pk, abs(d(SEL))); end
for i = 1:round(Thold/dt), d = P*d;            pk = max(pk, abs(d(SEL))); end
end

function [ok, d_end, u_h, pk, wb] = march(V, M, dXE, trim_lon, UH, k, a, d0, dt, SEL, W_LIM) %#ok<INUSD>
% March segment k at acceleration a, starting from the error d0 that the
% earlier segments left behind.
%
% The segment is safe when the aircraft is COVERED THE WHOLE WAY: inside set k
% from the start until some handover step, and inside set k+1 from that step
% through to the end of the segment. Merely touching the overlap somewhere is
% not enough - the aircraft enters the segment already inside both sets, so
% that weaker test is satisfied at the first step no matter how hard it then
% accelerates, and it lets a run away to the bracket.
%
% The aircraft's absolute state is the interpolated trim plus the error, so
% the forward-speed lag needs no separate treatment: it is inside x_act.
P  = M.Phi{k};   dx = dXE(:, k);
du = du_of(UH, k);
N  = min(max(round(du / a / dt), 1), 40000);   % cap: a tiny a must not hang
d  = d0;
pk = abs(d0(SEL));

% --- march the 34-state recursion, recording the absolute LON state -------
% Sequential by nature. The value-function lookups are done afterwards in one
% vectorised call each; querying them inside this loop costs ~50x more.
XA = zeros(4, N);   vv = zeros(1, N);
dtl = trim_lon(:,k+1) - trim_lon(:,k);
for i = 1:N
    v    = UH(k) + a * (i*dt);                       % commanded speed
    xe   = trim_lon(:,k) + ((v - UH(k)) / du) * dtl;
    XA(:,i) = xe + d(SEL);                           % actual [u; w; q; theta]
    vv(i)   = v;
    pk = max(pk, abs(d(SEL)));
    d  = P*d - dx*a*dt;
end

wbv   = XA(1,:) .* sin(XA(4,:)) + XA(2,:) .* cos(XA(4,:));
wb    = [min(wbv); max(wbv)];

E1  = XA - trim_lon(:,k);
E2  = XA - trim_lon(:,k+1);
in1 = V{k}(  E1(1,:)', E1(2,:)', E1(3,:)', E1(4,:)') <= 0;
in2 = V{k+1}(E2(1,:)', E2(2,:)', E2(3,:)', E2(4,:)') <= 0;

% h1 = last step of the leading run inside set k
j = find(~in1, 1, 'first');
if isempty(j), h1 = N; else, h1 = j - 1; end
% h2 = first step of the trailing run inside set k+1
j = find(~in2, 1, 'last');
if isempty(j), h2 = 1; else, h2 = j + 1; end

ok = (h2 <= N) && (h2 <= h1 + 1);
if ok, u_h = vv(max(h2,1)); else, u_h = NaN; end
d_end = d;
end

function du = du_of(UH, k)
du = UH(k+1) - UH(k);
end

function [Ttot, okall] = replay(V, M, dXE, trim_lon, UH, du_seg, a_sched, dt, SEL, W_LIM)
% Fly a given schedule end to end, carrying the error, and report whether every
% segment still finds a handover speed.
d = zeros(34,1);  Ttot = 0;  okall = true;
for k = 1:numel(a_sched)
    if ~isfinite(a_sched(k)) || a_sched(k) <= 0, okall = false; return; end
    [ok, d] = march(V, M, dXE, trim_lon, UH, k, a_sched(k), d, dt, SEL, W_LIM);
    okall = okall && ok;
    Ttot  = Ttot + du_seg(k)/a_sched(k);
end
end

function [a_lo, a_hi] = validity_window(u, w_h, th, g_th, W_LIM)
% Closed form. w_body(theta) = u sin(theta) + w_h cos(theta) = R sin(theta+phi),
% so |w_body| <= W_LIM inverts to an interval in theta, and
% theta(a) = th - a*g_th turns that into an interval in a.
R   = hypot(u, w_h);
phi = atan2(w_h, u);
if W_LIM >= R || g_th <= 0
    a_lo = 0;  a_hi = inf;  return;
end
s = asin(W_LIM / R);
a_lo = max((th - s + phi) / g_th, 0);
a_hi = (th + s + phi) / g_th;
end
