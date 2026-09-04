% place_trim_points - where to put the trim points so that consecutive safety
% tubes overlap by enough for the vehicle to actually pass through the overlap.
%
% ---------------------------------------------------------------------------
% THE CONDITION
% ---------------------------------------------------------------------------
% Two consecutive tubes overlap along the connecting direction by
%
%       overlap = a_fwd(u_k) + a_bwd(u_k+1) - step
%
% where a_fwd is how far the k-th tube reaches forward in forward speed from
% its own trim point and a_bwd is how far the (k+1)-th reaches back. A
% negative overlap is a gap of that size, so one expression covers both the
% overlapping and the non-overlapping case.
%
% Overlap on its own only says a handover REGION exists. The vehicle does not
% fly the reference exactly - it is behind by e - so the region has to be wide
% enough to contain the whole spread of where the vehicle might be:
%
%       g(step) = a_fwd(u) + a_bwd(u+step) - step - 2*e   >=  0            (G)
%
%       e = e0(u) + c(u) * step / T
%
% e0 is the error that survives at a held command; the second term is the
% ramp-proportional lag, both from compute_lag_coefficient.m. g > 0 means
% margin left over, g < 0 means the vehicle can miss the handover.
%
% Placement is then: stand at u_k, take the LONGEST step that still satisfies
% (G), and repeat. That is a walk over a table, done before any flight - the
% output is a list of speeds.
%
% ---------------------------------------------------------------------------
% WHY g HAS EXACTLY ONE ROOT
% ---------------------------------------------------------------------------
% Growing the step costs 1 ft/s of overlap per ft/s of step directly, plus
% 2*c/T more through the error term. It buys back only da/du per ft/s, and the
% measured a(u) rises about 5 ft/s over the whole 160 ft/s envelope - a slope
% near 0.03. So g falls monotonically and crosses zero once. Bisection is safe.
%
% ---------------------------------------------------------------------------
% THE HARD LIMIT
% ---------------------------------------------------------------------------
% Shrinking the step to nothing sends the overlap to 2*a(u) and the ramp term
% to zero, leaving  2*a(u) - 2*e0(u). So no placement whatsoever exists where
%
%       a(u)  <=  e0(u)
%
% That is the band the vehicle cannot cross on the speed axis alone, and it is
% where a second scheduling parameter (tilt, on a vehicle that has one) would
% have to come in. Note this limit involves only e0 - the part of the error
% that does NOT shrink with the step.
%
% Output: logger/trim_placement.mat, logger/trim_placement.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt      = 0.01;
T_SEG   = 2.0;      % segment duration [s]
WH_IDX  = 3;
N_TRIM  = 20;
E0_WIND = 0.0;      % extra irreducible error to allow for wind/model error
                    % [ft/s]. 0 = nominal. Sweep at the end.

%% ---------------------------------------------------------------------
%  1. a(u) - how far each tube reaches in forward speed
%     Read straight off the stored 4-D value arrays: walk the u axis through
%     the tube's own centre (w, q, theta perturbations = 0, which are exact
%     grid points) and find where the value crosses zero.
%% ---------------------------------------------------------------------
spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
ic = [find(gv{1}==0), find(gv{2}==0), find(gv{3}==0), find(gv{4}==0)];

params = struct('filter_mode','off');  params.T_seg = T_SEG;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj.trim_lon(1, :);

a_fwd = zeros(1, N_TRIM);  a_bwd = zeros(1, N_TRIM);
for k = 1:N_TRIM
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    line = squeeze(S.data(:, ic(2), ic(3), ic(4)));   % V along the u axis
    a_fwd(k) = cross_zero(gv{1}, line, ic(1), +1);
    a_bwd(k) = cross_zero(gv{1}, line, ic(1), -1);
end
fprintf('--- tube half-width in forward speed ---\n');
fprintf('  k |    u    | a_fwd | a_bwd\n');
for k = 1:N_TRIM
    fprintf('%3d | %7.2f | %5.2f | %5.2f\n', k, u_anchor(k), a_fwd(k), a_bwd(k));
end
fprintf('a grows %.2f -> %.2f ft/s across the envelope (slope ~%.3f per ft/s)\n\n', ...
        a_fwd(1), a_fwd(end), (a_fwd(end)-a_fwd(1))/(u_anchor(end)-u_anchor(1)));

%% ---------------------------------------------------------------------
%  2. c(u) and e0(u)
%% ---------------------------------------------------------------------
L = load(fullfile(root, 'logger', 'lag_coefficient.mat'));
S = load(fullfile(root, 'logger', 'closed_loop_model_ramped.mat'));
c_u  = L.c_pk;                                   % [s]
e0_u = abs(S.M.xi_e(4, :) - S.M.u_anchor) + E0_WIND;   % [ft/s]

A_f = @(u) interp1(u_anchor, a_fwd, clampu(u, u_anchor), 'linear');
A_b = @(u) interp1(u_anchor, a_bwd, clampu(u, u_anchor), 'linear');
C_u = @(u) interp1(u_anchor, c_u,   clampu(u, u_anchor), 'linear');
E0  = @(u) interp1(u_anchor, e0_u,  clampu(u, u_anchor), 'linear');

G = @(u, D) A_f(u) + A_b(u+D) - D - 2*(E0(u) + C_u(u)*D/T_SEG);

%% ---------------------------------------------------------------------
%  3. is the CURRENT uniform schedule admissible?
%% ---------------------------------------------------------------------
fprintf('--- current uniform schedule, step = 8.44 ft/s, T = %.1f s ---\n', T_SEG);
fprintf('  k |    u    | overlap | 2*e   | g       | verdict\n');
fprintf('%s\n', repmat('-', 1, 56));
for k = 1:N_TRIM-1
    D  = u_anchor(k+1) - u_anchor(k);
    ov = A_f(u_anchor(k)) + A_b(u_anchor(k+1)) - D;
    ee = 2*(E0(u_anchor(k)) + C_u(u_anchor(k))*D/T_SEG);
    g  = ov - ee;
    if g >= 0, vd = 'ok'; else, vd = 'FAILS'; end
    fprintf('%3d | %7.2f | %7.2f | %5.2f | %+7.2f | %s\n', ...
            k, u_anchor(k), ov, ee, g, vd);
end

%% ---------------------------------------------------------------------
%  4. the hard limit:  a(u) > e0(u) ?
%% ---------------------------------------------------------------------
lim = a_fwd - e0_u;
fprintf('\n--- crossable at all?  a(u) - e0(u) must be > 0 ---\n');
fprintf('min %.2f ft/s at u = %.1f  -> %s\n\n', min(lim), ...
        u_anchor(find(lim==min(lim),1)), ...
        ternary(min(lim) > 0, 'every speed is crossable', ...
                              'THERE IS AN UNCROSSABLE BAND'));

%% ---------------------------------------------------------------------
%  5. walk
%% ---------------------------------------------------------------------
u_end = u_anchor(end);
u = 0;  us = 0;  steps = [];
while u < u_end - 1e-6 && numel(us) < 200
    if G(u, 1e-6) <= 0, error('no admissible step at u = %.2f', u); end
    lo = 1e-6;  hi = 40;
    if G(u, hi) > 0
        D = hi;
    else
        for it = 1:60
            mid = (lo+hi)/2;
            if G(u, mid) > 0, lo = mid; else, hi = mid; end
        end
        D = lo;
    end
    D = min(D, u_end - u);
    u = u + D;  us(end+1) = u;  steps(end+1) = D; %#ok<SAGROW>
end

fprintf('--- placement ---\n');
fprintf('  k |    u    |  step  | vs 8.44\n');
fprintf('%s\n', repmat('-', 1, 40));
for k = 1:numel(steps)
    fprintf('%3d | %7.2f | %6.2f | %+6.2f\n', k, us(k), steps(k), steps(k)-8.44);
end
fprintf('%3d | %7.2f |   --   |  (cruise)\n', numel(us), us(end));
fprintf('\ntrim points: %d  (uniform schedule uses %d)\n', numel(us), N_TRIM);
fprintf('segments   : %d  ->  transition time %.1f s at T = %.1f s/segment\n', ...
        numel(steps), numel(steps)*T_SEG, T_SEG);
fprintf('             uniform: %d segments -> %.1f s\n', ...
        N_TRIM-1, (N_TRIM-1)*T_SEG);

%% ---------------------------------------------------------------------
%  6. how much wind can this survive?
%% ---------------------------------------------------------------------
fprintf('\n--- sensitivity to the irreducible error e0 ---\n');
fprintf(' extra e0 | trim pts | segments | time [s] | note\n');
for w = [0 0.5 1.0 1.5 2.0 3.0 4.0]
    E0w = @(u) E0(u) + w;
    Gw  = @(u, D) A_f(u) + A_b(u+D) - D - 2*(E0w(u) + C_u(u)*D/T_SEG);
    uu = 0;  ns = 0;  bad = false;
    while uu < u_end - 1e-6 && ns < 400
        if Gw(uu, 1e-6) <= 0, bad = true; break; end
        lo = 1e-6;  hi = 40;
        if Gw(uu, hi) > 0, D = hi; else
            for it = 1:60
                mid = (lo+hi)/2;
                if Gw(uu, mid) > 0, lo = mid; else, hi = mid; end
            end
            D = lo;
        end
        uu = uu + min(D, u_end - uu);  ns = ns + 1;
    end
    if bad
        fprintf('%9.1f |    --    |    --    |    --    | uncrossable at u = %.1f\n', w, uu);
    else
        fprintf('%9.1f | %8d | %8d | %8.1f |\n', w, ns+1, ns, ns*T_SEG);
    end
end

%% ---------------------------------------------------------------------
%  7. figure
%% ---------------------------------------------------------------------
f = figure('Position', [100 100 1000 700]);
uq = linspace(0, u_end, 400);

subplot(3,1,1);
plot(u_anchor, a_fwd, 'o-', 'LineWidth', 1.4); hold on;
plot(u_anchor, e0_u, 's-', 'LineWidth', 1.4);
ylabel('[ft/s]'); grid on;
legend('tube half-width a(u)', 'irreducible error e_0(u)', 'Location','northwest');
title('what the tube gives, and what is lost before the ramp even starts');

subplot(3,1,2);
plot(uq, arrayfun(@(x) G(x, 8.44), uq), 'LineWidth', 1.6); hold on;
yline(0, 'k-'); grid on;
ylabel('g at step = 8.44'); xlabel('');
title('the current uniform schedule: g < 0 means the handover can be missed');

subplot(3,1,3);
stem(us(1:end-1), steps, 'filled'); hold on;
yline(8.44, 'r--', 'LineWidth', 1.2);
xlabel('forward speed [ft/s]'); ylabel('step [ft/s]'); grid on;
legend('placed step', 'uniform 8.44', 'Location','northwest');
title(sprintf('placement: %d trim points, %d segments', numel(us), numel(steps)));
saveas(f, fullfile(root, 'logger', 'trim_placement.png'));

save(fullfile(root, 'logger', 'trim_placement.mat'), ...
     'u_anchor', 'a_fwd', 'a_bwd', 'c_u', 'e0_u', 'us', 'steps', 'T_SEG');
fprintf('\nsaved logger/trim_placement.mat and logger/trim_placement.png\n');

% -------------------------------------------------------------------------
function d = cross_zero(g, line, i0, sgn)
% distance from grid point i0 to where `line` crosses zero, going in
% direction sgn. Linear interpolation between the bracketing grid points.
n = numel(g);  i = i0;
while true
    j = i + sgn;
    if j < 1 || j > n, d = abs(g(i) - g(i0)); return; end   % never crosses
    if line(j) >= 0
        f = -line(i) / (line(j) - line(i));
        d = abs(g(i) + f*(g(j)-g(i)) - g(i0));  return;
    end
    i = j;
end
end

function y = clampu(u, ua)
y = min(max(u, ua(1)), ua(end));
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
