% optimize_schedule_sim - search the acceleration schedule against the REAL
% simulator instead of a linear surrogate.
%
% ---------------------------------------------------------------------------
% WHY DROP THE LINEAR MODEL
% ---------------------------------------------------------------------------
% The linear search returns ~43 s, of which about 19 s is margin carried for
% the surrogate's own error rather than for anything the aircraft does. Three
% attempts to shrink that failed, all against the same wall:
%
%   secant linearisation   the marginal mode (|lambda| = 0.999) flips unstable
%                          at twice the shipped perturbation - 0.2 ft/s
%   dropping that mode     prediction gets 10x worse; it is load-bearing
%   restoring the dropped  no help: the discarded -dq*dw and gravity terms are
%   nonlinear terms        OPEN-loop physics, and Phi is a CLOSED-loop map that
%                          already contains the controller's response to them
%
% The literature confirms the structural reason. Tischler & Tobias, RDMR-AF-
% 16-01, Section 2.3.2, builds its A_aero and B_aero to hold "only the
% aerodynamic dimensional stability and control derivatives... not gravity,
% Coriolis, or kinematic terms", restoring those exactly in Sections 2.3.5 and
% 2.3.7. Crucially their A is the BARE AIRCRAFT - 9 rigid-body states, control
% loop outside. Ours collapses controller, allocation, actuators, aero and
% rigid body into one 34x34 map, so there is nowhere to put a nonlinear term.
%
% That paper also solves the opposite problem from ours: it builds a simulator
% out of point models when no nonlinear aero model exists. We have the
% nonlinear aero model. Stitching would be a downgrade.
%
% So the surrogate is dropped and the search is run against the simulator. The
% margin for surrogate error disappears because there is no surrogate.
%
% This does not break the standing constraint on the guidance law: what was
% ruled out was prediction INSIDE the controller, and fitting theory backwards
% to simulation results. The safety criterion here is unchanged - the same
% reachable sets, read at each anchor's own trim point - and the simulator is
% used to evaluate candidates, not to invent the criterion.
%
% ---------------------------------------------------------------------------
% WHAT MAKES IT AFFORDABLE
% ---------------------------------------------------------------------------
% A naive search would re-fly the whole mission per candidate. Two things cut
% that down:
%
%   - the settle at the descending hover trim (45 s, needed because the hover
%     trim point is itself a descent, so a vehicle at rest is 11.7 ft/s off it)
%     is identical for every candidate and is flown ONCE;
%   - the state entering segment k depends only on segments 1..k-1, so it is
%     cached and a change at segment k re-flies only k..19.
%
% ---------------------------------------------------------------------------
% THE SAFETY TEST - unchanged
% ---------------------------------------------------------------------------
% A segment passes when the aircraft is covered the whole way: inside anchor
% k's set from the first step until some handover step, and inside anchor
% k+1's set from there through the last. Each set is read at the deviation
% from its OWN trim point, in the heading frame the sets were built in, and
% never as a blend of two anchors - that is a set nobody computed.
%
% Output: logger/schedule_sim.mat, logger/schedule_sim.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 0. constants --------------------------------------------------------
dt       = 0.01;
WDOT     = 11.667;        % WH3 descent the sets were built at [ft/s]
T_SETTLE = 45.0;
ALT0     = 1600;
A_MAX    = 16;   % raised: segment 19 pinned the old bracket at 12
A_MIN    = 0.30;
N_BIS    = 8;
N_SWEEP  = 2;
% The transition does not end at the last anchor. A schedule that accelerates
% hard into segment 19 and stops overshoots afterwards - measured 10.2 ft/s
% past cruise, against a 9.67 ft/s half-width at anchor 20 - and an earlier
% version never saw it, because the check stopped at the last segment. With
% nothing downstream to constrain it, segment 19 ran straight to the search
% bracket. The hold below is part of the safety test, not a postscript.
T_CRUISE = 25.0;          % hold cruise and keep checking [s]
SEL      = [4 6 11 8];
LON      = [1 3 5 11];
line = @(c) fprintf('%s\n', repmat(c, 1, 90));

T   = load('trim_table_Poly_ConcatVer4p0.mat');
UH  = T.UH(1:20);  UH = UH(:)';
trim_lon = T.XU0_interp(LON, 1:20, 3);
du_seg   = diff(UH);
NS = 19;

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
V = cell(1,20);
for k = 1:20
    S = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat',k,3)), 'data');
    V{k} = griddedInterpolant(gv, S.data, 'linear', 'nearest');
end

params = struct('filter_mode','off');  params.T_seg = 2.0;
guam = LpC_GUAM(Config('trim_schedule', params));

ctx = struct('guam',guam, 'V',{V}, 'TL',trim_lon, 'UH',UH, 'du',du_seg, ...
             'dt',dt, 'WDOT',WDOT, 'NS',NS, 'T_CRUISE',T_CRUISE, ...
             'GMIN',spec.grid_min(:), 'GMAX',spec.grid_max(:));

%% 1. settle once at the descending hover trim -------------------------
line('='); fprintf('SETTLING AT THE DESCENDING HOVER TRIM (%.0f s)\n', T_SETTLE); line('=');
guam.reset();
s = guam.saveState();  s.state(3) = -ALT0;  guam.restoreState(s);
p = [0; 0; -ALT0];
tic;
for i = 1:round(T_SETTLE/dt)
    p = p + [0; 0; WDOT*dt];
    guam.step(struct('pos',p, 'vel',[0;0;WDOT], 'chi',0, 'chi_dot',0));
end
st = guam.state;
uh =  st(4)*cos(st(8)) + st(6)*sin(st(8));
wh = -st(4)*sin(st(8)) + st(6)*cos(st(8));
fprintf('  settled: u_h %.2f, w_h %.2f (trim %.3f), theta %.2f deg  [%.0f s wall]\n\n', ...
        uh, wh, trim_lon(2,1), rad2deg(st(8)), toc);
seed = struct('S', guam.saveState(), 'p', p);

%% 2. timing, to size the search --------------------------------------
a0 = 2.5 * ones(1, NS);
tic;  [ok0, ~] = fly_from(ctx, 1, seed, a0);  t_full = toc;
fprintf('  one full-mission evaluation: %.1f s wall, uniform 2.5 %s\n', ...
        t_full, ternary(ok0, 'PASSES', 'FAILS'));
fprintf('  budget: ~%d evaluations planned -> ~%.0f min\n\n', ...
        NS*(N_BIS+1)*N_SWEEP, NS*(N_BIS+1)*N_SWEEP*t_full*0.6/60);

%% 3. seed: largest uniform acceleration that passes -------------------
line('='); fprintf('SEED  (largest uniform acceleration)\n'); line('=');
lo = A_MIN;  hi = 8;
assert(fly_from(ctx, 1, seed, lo*ones(1,NS)), 'optimize_schedule_sim:noSeed', ...
    'Even %.2f ft/s^2 uniform fails.', lo);
for i = 1:6
    m = (lo+hi)/2;
    if fly_from(ctx, 1, seed, m*ones(1,NS)), lo = m; else, hi = m; end
    fprintf('  tried %.2f -> %s\n', m, ternary(lo==m,'pass','fail'));
end
a = lo * ones(1, NS);
[ok, cache] = fly_from(ctx, 1, seed, a);
assert(ok, 'optimize_schedule_sim:seed', 'Seed re-check failed.');
fprintf('  uniform %.2f passes, T = %.2f s\n\n', lo, sum(du_seg./a));

%% 4. coordinate ascent ------------------------------------------------
line('='); fprintf('COORDINATE ASCENT AGAINST THE SIMULATOR\n'); line('=');
t0 = tic;
for sw = 1:N_SWEEP
    [~, order] = sort(du_seg ./ a, 'descend');
    moved = 0;
    for k = order
        base = a(k);
        trial = a;  trial(k) = A_MAX;
        if fly_from(ctx, k, cache(k), trial)
            best = A_MAX;
        else
            plo = base;  phi = A_MAX;
            for i = 1:N_BIS
                m = (plo+phi)/2;
                trial = a;  trial(k) = m;
                if fly_from(ctx, k, cache(k), trial), plo = m; else, phi = m; end
            end
            best = plo;
        end
        if best > base + 5e-3
            a(k) = best;
            [okk, cache] = fly_from(ctx, 1, seed, a);
            assert(okk, 'optimize_schedule_sim:lost', 'Accepted a raise that fails.');
            moved = moved + 1;
            fprintf('    seg %2d : %.2f -> %.2f   (T = %.2f s, %.0f s wall)\n', ...
                    k, base, best, sum(du_seg./a), toc(t0));
        end
    end
    fprintf('  sweep %d done: %d raised, T = %.2f s\n', sw, moved, sum(du_seg./a));
    if moved == 0, break; end
end
T_min = sum(du_seg ./ a);

%% 5. final flight and margins ----------------------------------------
line('='); fprintf('FINAL SCHEDULE, FLOWN\n'); line('=');
[okf, ~, diag] = fly_from(ctx, 1, seed, a);
fprintf('   k   u_trim      a     T_k [s]   worst V   dist to boundary [ft/s]\n');
for k = 1:NS
    fprintf('%4d %8.1f %7.2f %9.2f %10.4f %19.2f\n', k, UH(k), a(k), ...
            du_seg(k)/a(k), diag.Vw(k), v_to_dist(V{k}, diag.Vw(k)));
end
fprintf('\n  all segments covered : %d\n', okf);
fprintf('  transition time      : %.2f s\n', T_min);
fprintf('  schedule [ft/s^2]    : ');  fprintf('%.2f ', a);  fprintf('\n');

% accel_schedule.mat holds whatever the last surrogate run wrote, and the
% per-channel sweeps overwrite it - reading T_min from there once reported the
% zero-uncertainty run (24.18 s) as if it were the real surrogate answer. The
% number to compare against is pinned here so the comparison cannot drift.
T_SURROGATE = 42.99;   % linear surrogate carrying the MEASURED uncertainty
sL = load(fullfile(root,'logger','accel_schedule.mat'));
fprintf('\n  surrogate, measured uncertainty : %6.2f s\n', T_SURROGATE);
fprintf('  surrogate, perfect model assumed: %6.2f s   (last accel_schedule.mat run)\n', sL.T_min);
fprintf('  simulator search                : %6.2f s\n', T_min);
fprintf('  recovered %.2f s that the surrogate had to carry as margin\n', T_SURROGATE - T_min);

save(fullfile(root,'logger','schedule_sim.mat'), 'a','T_min','UH','du_seg','diag');

f = figure('Position',[80 80 1100 400], 'Color','w');
subplot(1,3,1); stairs(UH(1:NS), a, 'LineWidth',2.2); hold on;
stairs(UH(1:NS), sL.a, '--', 'LineWidth',1.6); grid on;
xlabel('u [ft/s]'); ylabel('a [ft/s^2]'); title('acceleration per segment');
legend({'simulator search','linear surrogate'}, 'Location','northwest');
subplot(1,3,2); bar(UH(1:NS), du_seg./a, 'FaceColor',[.35 .55 .62]); grid on;
xlabel('u [ft/s]'); ylabel('T_k [s]'); title(sprintf('total %.1f s', T_min));
subplot(1,3,3); plot(UH(1:NS), diag.Vw, '-o', 'LineWidth',1.8); hold on;
yline(0,'r--'); grid on; xlabel('u [ft/s]'); ylabel('worst V');
title('margin at the worst moment');
saveas(f, fullfile(root,'logger','schedule_sim.png'));
fprintf('\nsaved logger/schedule_sim.mat and logger/schedule_sim.png\n');

%% ---------------------------------------------------------------------
function [ok, cache, diag] = fly_from(ctx, k0, start, a)
% Fly segments k0..19 from a cached simulator state. Returns whether every
% segment stayed covered, plus the state entering each later segment so a
% subsequent test can restart there instead of from hover.
g = ctx.guam;   dt = ctx.dt;
g.restoreState(start.S);   p = start.p;
cache = repmat(struct('S',[],'p',[]), 1, ctx.NS);
ok = true;   diag.Vw = nan(1, ctx.NS);
for k = k0:ctx.NS
    cache(k).S = g.saveState();   cache(k).p = p;
    N = max(round(ctx.du(k)/a(k)/dt), 1);
    XB = zeros(4, N);
    for i = 1:N
        s = g.state;
        XB(:, i) = [s(4); s(6); s(11); s(8)];
        v = ctx.UH(k) + a(k)*(i*dt);
        p = p + [v*dt; 0; ctx.WDOT*dt];
        g.step(struct('pos',p, 'vel',[v;0;ctx.WDOT], 'chi',0, 'chi_dot',0));
    end
    th = XB(4,:);
    XA = [ XB(1,:).*cos(th) + XB(2,:).*sin(th) ;      % body -> heading
          -XB(1,:).*sin(th) + XB(2,:).*cos(th) ;
           XB(3,:) ; th ];
    E1 = XA - ctx.TL(:,k);
    E2 = XA - ctx.TL(:,k+1);
    v1 = ctx.V{k}(  E1(1,:)', E1(2,:)', E1(3,:)', E1(4,:)');
    v2 = ctx.V{k+1}(E2(1,:)', E2(2,:)', E2(3,:)', E2(4,:)');
    v1(any(E1 < ctx.GMIN | E1 > ctx.GMAX, 1)') = 1;   % outside the grid is
    v2(any(E2 < ctx.GMIN | E2 > ctx.GMAX, 1)') = 1;   % not certifiable
    j = find(v1 > 0, 1, 'first');  if isempty(j), h1 = N; else, h1 = j-1; end
    j = find(v2 > 0, 1, 'last');   if isempty(j), h2 = 1; else, h2 = j+1; end
    diag.Vw(k) = max(min(v1, v2));
    if ~((h2 <= N) && (h2 <= h1 + 1)) || any(~isfinite(XB(:)))
        ok = false;  return;
    end
end

% --- hold cruise and keep checking against the last anchor's own set -----
% Only when the whole transition was flown; a partial re-fly from segment k
% still ends at the last anchor, so this runs either way.
Nc = round(ctx.T_CRUISE/ctx.dt);
XC = zeros(4, Nc);
for i = 1:Nc
    s = g.state;
    XC(:, i) = [s(4); s(6); s(11); s(8)];
    p = p + [ctx.UH(end)*dt; 0; ctx.WDOT*dt];
    g.step(struct('pos',p, 'vel',[ctx.UH(end);0;ctx.WDOT], 'chi',0, 'chi_dot',0));
end
thc = XC(4,:);
XAc = [ XC(1,:).*cos(thc) + XC(2,:).*sin(thc) ;
       -XC(1,:).*sin(thc) + XC(2,:).*cos(thc) ;
        XC(3,:) ; thc ];
Ec = XAc - ctx.TL(:,end);
vc = ctx.V{end}(Ec(1,:)', Ec(2,:)', Ec(3,:)', Ec(4,:)');
vc(any(Ec < ctx.GMIN | Ec > ctx.GMAX, 1)') = 1;
diag.Vcruise = max(vc);
if diag.Vcruise > 0 || any(~isfinite(XC(:)))
    ok = false;  return;
end
end

function dist = v_to_dist(Vk, Vval)
% V is nearly flat inside (about -0.21 at the trim point, rising to 0 only
% over the outer quarter), so a raw value says little. Convert to a distance
% from the boundary along the speed axis.
if ~isfinite(Vval) || Vval >= 0, dist = 0; return; end
dist = max(ray_hit(Vk, 0) - ray_hit(Vk, Vval), 0);
end

function r = ray_hit(Vk, level)
lo = 0; hi = 16;
for i = 1:40
    m = (lo+hi)/2;
    if Vk(m,0,0,0) <= level, lo = m; else, hi = m; end
end
r = lo;
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
