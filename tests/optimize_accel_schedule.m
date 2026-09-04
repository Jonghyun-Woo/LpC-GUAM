% optimize_accel_schedule - push every segment's acceleration to the largest
% value that keeps the WHOLE transition inside the reachable sets.
%
% ---------------------------------------------------------------------------
% WHY A SEARCH IS NEEDED AT ALL
% ---------------------------------------------------------------------------
% Each segment would like its own largest acceleration, and if the segments
% were independent that would be 19 separate one-line answers. They are not.
% The closed loop carries a mode with |lambda| = 0.9991, a time constant of
% 12 s and a period of 59 s, living in the position error and the longitudinal
% integrator. Over a 33-40 s transition that mode accumulates across segment
% boundaries, so the acceleration chosen at segment 11 decides whether segment
% 17 survives.
%
% Measured, not argued: with the per-segment answers in place, segment 17 fails
% no matter how gently it is flown - lowering a(17) from 3.41 all the way to
% 0.80 does not help. Slowing segments 11-16 by 30 % does. The failure is
% upstream of where it shows up.
%
% Taking the largest acceleration segment by segment as the march goes forward
% fails for the same reason, and fails hard: segment 1 grabs 8.93 ft/s^2,
% finishes in 0.95 s, and the error peak - which arrives 2.75 s after a rate
% change - lands inside segment 2, which then has no feasible acceleration at
% all. Greedy is blind to damage it does downstream.
%
% So the schedule is searched, not derived. This is NOT the circular
% definition that was rejected earlier (offset needs schedule, schedule needs
% offset): given a schedule, the error is computed exactly in one forward
% pass. What is left is an ordinary search over 19 numbers.
%
% ---------------------------------------------------------------------------
% THE SEARCH
% ---------------------------------------------------------------------------
% Coordinate ascent from a known-feasible schedule. Segments are visited
% slowest-first, because the time saved by raising a segment is
% d(du/a)/da = -du/a^2, largest where a is smallest. Raising one segment is
% accepted only if the whole remainder of the transition still passes.
%
% The error entering segment k is cached, so testing a change at segment k
% only re-marches segments k..19 rather than the whole transition.
%
% ---------------------------------------------------------------------------
% THE SAFETY TEST
% ---------------------------------------------------------------------------
% A segment passes when the aircraft is COVERED THE WHOLE WAY: inside set k
% from the first step until some handover step, and inside set k+1 from that
% step through the last. Merely touching the overlap somewhere is not enough -
% the aircraft enters a segment already inside both sets, so the weaker test
% is satisfied at step one however hard it then accelerates.
%
% The state is the 34-state closed loop marched by d(i+1) = Phi d(i) -
% xi_e'(v) a dt, and the aircraft's absolute longitudinal state is the
% interpolated SETTLED state xi_e plus that error. Nothing is simulated.
%
% ---------------------------------------------------------------------------
% THE FRAME, AND WHY xi_e AND NOT THE TRIM TABLE
% ---------------------------------------------------------------------------
% xi_e and the closed-loop state are BODY frame; the trim table and the BRT
% axes are HEADING frame. They are rotated before any lookup. Adding the two
% directly - which an earlier version did - silently mixes an 11.7 ft/s
% standing offset into w and 2-4 deg into theta.
%
% The offset is real, not an artefact. The mission holds 80 ft (level), while
% the trim points and their reachable sets were built at WH3, a 11.667 ft/s
% descent. Flown level, the aircraft's heading-frame w measures 0.00 to -0.23
% ft/s across the whole envelope - confirmed against the real simulator - so
% it sits a constant -11.7 ft/s off every anchor before it accelerates at all.
% The sets absorb it (V stays near -0.20) but it costs 10-21 %% of the usable
% speed half-width. Starting from xi_e rather than the trim table brings that
% offset in automatically.
%
% Output: logger/accel_schedule.mat, logger/accel_schedule.png
% EPS_SCALE multiplies the measured model uncertainty. It exists to ANSWER
% "how much is the uncertainty costing us", not to make the answer look good:
% shrinking it does not make the aircraft faster, it removes the accounting
% that says the answer means anything. The honest ways to reduce the
% uncertainty are to measure Phi at realistic perturbation sizes, or to
% re-measure eps on schedules shaped like the ones actually flown.
if ~exist('EPS_SCALE','var') || isempty(EPS_SCALE), EPS_SCALE = 1.0; end
clearvars -except EPS_SCALE;  close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 0. constants --------------------------------------------------------
WH_IDX = 3;   K = 20;   NS = K - 1;   KT = 1.6878;
% The search may not go past the largest acceleration the model uncertainty
% was MEASURED at. Beyond a = 8 there is no measurement, and holding the
% uncertainty at its last measured value would understate it - the search then
% buys speed with an error nobody has bounded. An earlier run put segment 8 at
% exactly this bracket, which is the symptom to watch for.
A_MAX  = 8;
A_FLOOR = 0.20;
N_BIS  = 11;              % bisection steps per segment visit
N_SWEEP = 3;
SEL = [4 6 11 8];
LON = [1 3 5 11];
line = @(c) fprintf('%s\n', repmat(c, 1, 88));

%% 1. inputs -----------------------------------------------------------
T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:K);  UH = UH(:)';
WH3 = T.WH(WH_IDX);
XU0 = T.XU0_interp(:, 1:K, WH_IDX);
trim_lon = XU0(LON, :);
du_seg   = diff(UH);

% The DESCENDING model. Its equilibrium sits on the trim points to within
% 0.010 ft/s and 0.017 deg, so xi_e and the anchors agree and the only offset
% left is the tracking error. Built by tests/build_descending_model.m.
Sm = load(fullfile(root,'logger','closed_loop_model_descending.mat'));
M = Sm.M;  dXE = Sm.dXE;  dt = M.dt;
fprintf('model built at w_ref = %.3f ft/s\n', M.w_ref);

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
V = cell(1,K);
for k = 1:K
    S = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat',k,WH_IDX)), 'data');
    V{k} = griddedInterpolant(gv, S.data, 'linear', 'nearest');
end

% MODEL UNCERTAINTY, MEASURED NOT GUESSED
% ----------------------------------------
% tests/measure_model_validity flew the simulator and the linear model side by
% side on the same descending command and recorded how far apart they ended
% up, per channel, at four uniform accelerations. Those numbers are the
% columns below: rows are [u; w; q; theta] in ft/s, ft/s, rad/s, rad.
%
% Read as a fraction of the reachable set's half-width the disagreement is
% modest where it matters - at most 8 %% at a = 2 and 18 %% at a = 3 - and only
% becomes serious at a = 5 (42 %%), by which point the true error has already
% passed the half-width (108 %% in u) and the schedule is infeasible anyway.
%
% An earlier version capped the error at a hand-picked 4 ft/s instead. Every
% segment then sat exactly on that cap and the answer, 133.8 s, was a property
% of the guess rather than of the vehicle.
EPS_A   = [2      3      5      8    ];
EPS_TAB = [0.41   0.71   2.05   5.87 ;      % u     [ft/s]
           1.98   3.72  10.45  19.23 ;      % w     [ft/s]
           0.0147 0.0227 0.0456 0.0892;     % q     [rad/s]
           0.0220 0.0545 0.0918 0.1026];    % theta [rad]
% EPS_SCALE may be a scalar (all channels) or a 4-vector (per channel:
% u, w, q, theta). The per-channel form is how the 19 s the uncertainty costs
% gets attributed - zero one channel and the difference is what fixing the
% model in that channel would buy.
EPS_TAB = EPS_TAB .* EPS_SCALE(:);

ctx = struct('V',{V}, 'Phi',{M.Phi}, 'dXE',dXE, 'TL',trim_lon, 'XE',M.xi_e(SEL,:), ...
             'UH',UH, 'du',du_seg, 'dt',dt, 'SEL',SEL, 'NS',NS, ...
             'GMIN',spec.grid_min(:), 'GMAX',spec.grid_max(:), ...
             'EPS_A',EPS_A, 'EPS_TAB',EPS_TAB);

%% 2. a feasible starting point ---------------------------------------
% Found rather than assumed: the largest UNIFORM acceleration the whole
% transition survives. Uniform is a poor schedule but a reliable seed, and a
% hardcoded one dies silently every time the closed-loop model changes.
line('='); fprintf('STARTING POINT  (largest uniform acceleration)\n'); line('=');
lo = 0.2;  hi = 8;
assert(check_from(ctx, 1, zeros(34,1), lo*ones(1,NS)), ...
    'optimize_accel_schedule:noSeed', ...
    'Even %.2f ft/s^2 uniform fails; nothing is feasible.', lo);
for i = 1:12
    m = (lo+hi)/2;
    if check_from(ctx, 1, zeros(34,1), m*ones(1,NS)), lo = m; else, hi = m; end
end
a = lo * ones(1, NS);
[ok0, dent] = check_from(ctx, 1, zeros(34,1), a);
assert(ok0, 'optimize_accel_schedule:seedInfeasible', 'Seed re-check failed.');
fprintf('  uniform %.2f ft/s^2 passes, T = %.2f s\n\n', lo, sum(du_seg./a));

%% 3. coordinate ascent ------------------------------------------------
line('='); fprintf('COORDINATE ASCENT  (slowest segment first)\n'); line('=');
t0 = tic;
for sw = 1:N_SWEEP
    [~, order] = sort(du_seg ./ a, 'descend');       % biggest time sink first
    moved = 0;
    for k = order
        base = a(k);
        % try the ceiling first; if that passes there is nothing to bisect
        trial = a;  trial(k) = A_MAX;
        if check_from(ctx, k, dent(:,k), trial)
            best = A_MAX;
        else
            lo = base;  hi = A_MAX;
            for i = 1:N_BIS
                m = (lo+hi)/2;
                trial = a;  trial(k) = m;
                if check_from(ctx, k, dent(:,k), trial), lo = m; else, hi = m; end
            end
            best = lo;
        end
        if best > base + 1e-3
            a(k) = best;
            [okk, dent] = check_from(ctx, 1, zeros(34,1), a);
            assert(okk, 'optimize_accel_schedule:lostFeasibility', ...
                'Accepted a raise that does not pass a full re-check.');
            moved = moved + 1;
        end
    end
    fprintf('  sweep %d : %2d segments raised, T = %.2f s   (%.0f s elapsed)\n', ...
            sw, moved, sum(du_seg./a), toc(t0));
    if moved == 0, break; end
end

T_min = sum(du_seg ./ a);

%% 4. is every segment now at its ceiling? -----------------------------
% For each segment, confirm that a further nudge upward breaks the transition.
% If any segment can still be raised, the ascent has not converged.
line('='); fprintf('CEILING CHECK  -  can any single segment go faster?\n'); line('=');
fprintf('   k   u_trim      a      T_k     +1%% still safe?   time it would save\n');
at_ceiling = true(1, NS);
for k = 1:NS
    trial = a;  trial(k) = a(k) * 1.01;
    can = check_from(ctx, k, dent(:,k), trial);
    at_ceiling(k) = ~can;
    if can
        save_s = du_seg(k)/a(k) - du_seg(k)/trial(k);
        st = sprintf('YES  <-- not converged   %.3f s', save_s);
    else
        st = 'no (at the limit)';
    end
    fprintf('%4d %8.1f %7.2f %8.2f   %s\n', k, UH(k), a(k), du_seg(k)/a(k), st);
end
fprintf('\n  segments still at their ceiling : %d / %d\n', sum(at_ceiling), NS);

%% 5. margins ----------------------------------------------------------
% V is nearly flat in the interior (about -0.21 at the trim point, rising to 0
% only over the outer quarter), so a raw V value says little. Each segment's
% worst V is converted to a distance from the boundary along the speed axis.
line('='); fprintf('MARGIN AT THE WORST MOMENT OF EACH SEGMENT\n'); line('=');
fprintf('   k   u_trim      a     worst V    dist to boundary [ft/s]   peak w   peak th\n');
[Vw, dist, pkw, pkt] = deal(zeros(1, NS));
for k = 1:NS
    [~, ~, Vw(k), pk] = seg_march(ctx, k, a(k), dent(:,k));
    dist(k) = v_to_dist(V{k}, Vw(k));
    pkw(k)  = pk(2);   pkt(k) = rad2deg(pk(4));
    fprintf('%4d %8.1f %7.2f %10.3f %19.2f %11.2f %8.2f\n', ...
            k, UH(k), a(k), Vw(k), dist(k), pkw(k), pkt(k));
end
[dmin, kd] = min(dist);

%% 6. result -----------------------------------------------------------
line('='); fprintf('RESULT\n'); line('=');
fprintf('  transition time            : %.2f s\n', T_min);
fprintf('  thinnest margin            : %.2f ft/s at segment %d\n', dmin, kd);
fprintf('  peak error over transition : w %.2f ft/s, theta %.2f deg\n', ...
        max(pkw), max(pkt));
fprintf('\n  schedule (ft/s^2):\n   ');
fprintf('%.2f ', a);  fprintf('\n');

f = figure('Position',[80 80 1120 420], 'Color','w');
subplot(1,3,1);
stairs(UH(1:NS), a, 'LineWidth', 2.2); grid on;
xlabel('u [ft/s]'); ylabel('a [ft/s^2]'); title('acceleration per segment');
subplot(1,3,2);
bar(UH(1:NS), du_seg./a, 'FaceColor', [.35 .55 .62]); grid on;
xlabel('u [ft/s]'); ylabel('T_k [s]'); title(sprintf('segment time, total %.1f s', T_min));
subplot(1,3,3);
plot(UH(1:NS), dist, '-o', 'LineWidth', 1.8); grid on;
xlabel('u [ft/s]'); ylabel('distance to boundary [ft/s]');
title('margin at the worst moment');
saveas(f, fullfile(root,'logger','accel_schedule.png'));

save(fullfile(root,'logger','accel_schedule.mat'), ...
     'a','UH','du_seg','T_min','dent','Vw','dist','pkw','pkt','at_ceiling','trim_lon');
fprintf('\nsaved logger/accel_schedule.mat and logger/accel_schedule.png\n');

%% ---------------------------------------------------------------------
%  helpers
%% ---------------------------------------------------------------------
function [ok, dent] = check_from(ctx, k0, d0, a)
% Does the transition survive from segment k0 to the end, entering k0 with the
% error d0? Stops at the first failure - only the verdict is needed.
dent = zeros(34, ctx.NS);
d = d0;  ok = true;
for k = k0:ctx.NS
    dent(:,k) = d;
    [okk, d] = seg_march(ctx, k, a(k), d);
    if ~okk, ok = false; return; end
end
end

function [ok, d_end, Vworst, pk] = seg_march(ctx, k, a, d0)
% March one segment and test coverage. The value functions are queried once,
% vectorised over the whole segment; querying inside the loop costs ~50x more.
P  = ctx.Phi{k};   dx = ctx.dXE(:,k);   dt = ctx.dt;   du = ctx.du(k);
N  = min(max(round(du/a/dt), 1), 40000);
d  = d0;   pk = abs(d0(ctx.SEL));
XB = zeros(4, N);   d_hist = zeros(4, N);
dxe = ctx.XE(:,k+1) - ctx.XE(:,k);
for i = 1:N
    v = ctx.UH(k) + a*(i*dt);
    d_hist(:,i) = d(ctx.SEL);
    XB(:,i) = ctx.XE(:,k) + ((v - ctx.UH(k))/du)*dxe + d(ctx.SEL);
    pk = max(pk, abs(d(ctx.SEL)));
    d  = P*d - dx*a*dt;
end
% xi_e and the closed-loop state are BODY frame; the trim table and the BRT
% axes are HEADING frame. Rotate before looking anything up. Skipping this
% mixes an 11.7 ft/s standing offset into w and 2-4 deg into theta.
th = XB(4,:);
XA = [ XB(1,:).*cos(th) + XB(2,:).*sin(th) ;
      -XB(1,:).*sin(th) + XB(2,:).*cos(th) ;
       XB(3,:) ;
       th ];
E1 = XA - ctx.TL(:,k);
E2 = XA - ctx.TL(:,k+1);

% The linear model does not land exactly where the aircraft does. Push the
% predicted deviation OUTWARD by the measured disagreement before looking it
% up, so a point only passes if it would still be inside had the model been
% wrong by as much as it was measured to be. Outward is the worst direction
% for a set that contains its own centre, so one inflated lookup replaces the
% sixteen corners of the uncertainty box.
epsm = model_eps(ctx, a);
E1 = E1 + sign(E1) .* epsm;
E2 = E2 + sign(E2) .* epsm;

v1 = ctx.V{k}(  E1(1,:)', E1(2,:)', E1(3,:)', E1(4,:)');
v2 = ctx.V{k+1}(E2(1,:)', E2(2,:)', E2(3,:)', E2(4,:)');

% Anything outside the stored grid is not certifiable. griddedInterpolant is
% set to 'nearest' past the edges, so it silently returns the boundary value;
% that value is strongly positive over most of each face (9.15 at anchor 1) but
% 2 %% of anchor 10's w face is negative, so the clamp cannot be relied on.
bad1 = any(E1 < ctx.GMIN | E1 > ctx.GMAX, 1)';
bad2 = any(E2 < ctx.GMIN | E2 > ctx.GMAX, 1)';
v1(bad1) = 1;   v2(bad2) = 1;

j = find(v1 > 0, 1, 'first');   if isempty(j), h1 = N; else, h1 = j-1; end
j = find(v2 > 0, 1, 'last');    if isempty(j), h2 = 1; else, h2 = j+1; end
ok = (h2 <= N) && (h2 <= h1 + 1);
Vworst = max(min(v1, v2));
d_end  = d;
end

function e = model_eps(ctx, a)
% Measured model-vs-simulator disagreement at acceleration a, per channel.
% Below the measured range it is scaled down in proportion to a - the error is
% driven by how fast the command moves - and above it, held at the last
% column, though the search never gets there because the set runs out first.
A = ctx.EPS_A;  Tb = ctx.EPS_TAB;
if a <= A(1)
    e = Tb(:,1) * (a / A(1));
elseif a >= A(end)
    e = Tb(:,end);
else
    e = interp1(A, Tb', a)';
end
end

function dist = v_to_dist(Vk, Vval)
% How far from the boundary, measured along the speed axis from the trim
% point, does the level set V = Vval sit? Reported as the remaining distance
% to V = 0. V is not a signed distance, so this conversion is what makes the
% margins comparable between segments.
if Vval >= 0, dist = 0; return; end
r0 = ray_hit(Vk, 0);
rv = ray_hit(Vk, Vval);
dist = max(r0 - rv, 0);
end

function r = ray_hit(Vk, level)
lo = 0; hi = 16;
for i = 1:40
    m = (lo+hi)/2;
    if Vk(m,0,0,0) <= level, lo = m; else, hi = m; end
end
r = lo;
end
