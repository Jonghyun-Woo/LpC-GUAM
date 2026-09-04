% accel_trim_envelope - steady-acceleration envelope at each trim point.
%
% For each of the 20 transition trim points, compute:
%   (1) the range of forward acceleration a that the actuators can hold in a
%       steady sense (wdot = 0, qdot = 0, q = 0), and
%   (2) for each a in that range, the range of pitch-attitude perturbation
%       dtheta that is compatible with it.
%
% The equilibrium is linearised about the trim point, so with u, w and q
% perturbations pinned at zero the only free state is dtheta:
%
%     Bp(1:3,:) * ddelta  +  Ap(1:3,4) * dtheta  =  [a; 0; 0]
%
%   unknowns : 11 actuator perturbations + dtheta   = 12
%   equations: 3
%   -> 9 degrees of freedom. Many actuator combinations give the same a, and
%      they differ in how much of the job the pitch attitude does. The LP
%      walks that freedom to its two extremes.
%
% No reachable-set data is used here. The output table is meant to be read
% later, when the BRT overlap constraint is laid on top of it.
%
% WH is fixed at index 3 (wbar = +11.667 ft/s, a 700 ft/min descent). No
% other row of the trim table is touched.
%
% Output: data/accel_trim_envelope_WH3.mat, logger/accel_envelope*.png

%% 0. guards and constants ------------------------------------------------
% Bail out early if the LP solver is missing, and fix every constant here.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

if isempty(which('linprog'))
    error('accel_trim_envelope:noToolbox', ...
        ['MATLAB linprog was not found. This script needs the Optimization\n' ...
         'Toolbox. Install/licence it, or replace the LP calls with an\n' ...
         'equivalent solver, before running again.']);
end

WH_IDX  = 3;        % wbar = +11.667 ft/s
K       = 20;       % transition trim points (UH 1..20)
NA      = 41;       % points in the acceleration grid
NU      = 11;       % actuators: lift 1-8, pusher, elevator, flap
G       = 32.174;   % ft/s^2
A_PROBE = [3.38 4.69 5.63];   % accelerations the band plot / summary use
DTH_LIN_FLAG = 0.2; % |dtheta| beyond this (rad) is outside linear validity
SPLICE_WARN  = 17:19;         % segments known to sit near the table splice

% Actuator absolute limits. Matches RSLQRConfig (eng_max = [1600x8; 2000] rpm,
% ele/flp = +-30 deg); repeated here so this script states its own envelope.
LIM_LO = [zeros(8,1);  0;      -deg2rad(30); -deg2rad(30)];
LIM_HI = [167.55*ones(8,1); 209.44; deg2rad(30);  deg2rad(30)];
ACT_NAME = {'lift1','lift2','lift3','lift4','lift5','lift6','lift7','lift8', ...
            'pusher','elevator','flap'};

opts = optimoptions('linprog','Display','none');

%% 1. load the trim table -------------------------------------------------
% Read-only. The file lives in tables/trim/ and is on the path via genpath.
S = load('trim_table_Poly_ConcatVer4p0.mat');
UH_used = S.UH(1:K);
XU0     = S.XU0_interp(:, 1:K, WH_IDX);      % 25 x 20
AP      = S.Ap_lon_interp(:, :, 1:K, WH_IDX);   % 4 x 4  x 20   [u w q theta]
BP      = S.Bp_lon_interp(:, :, 1:K, WH_IDX);   % 4 x 11 x 20
WL      = S.W_lon_interp(:, :, 1:K, WH_IDX);    % 12 x 12 x 20
BL      = S.B_lon_interp(:, :, 1:K, WH_IDX);    % 3 x 12 x 20

fprintf('trim table  : %s\n', which('trim_table_Poly_ConcatVer4p0.mat'));
fprintf('WH index %d  : wbar = %+.4f ft/s\n', WH_IDX, S.WH(WH_IDX));
fprintf('trim points : %d   (UH %.1f .. %.1f ft/s)\n\n', K, UH_used(1), UH_used(end));

%% 2. actuator trim values and box bounds in perturbation coordinates ------
% Bp column order is [lift1..8, pusher, elevator, flap]; the matching trim
% values sit in XU0 rows 17:24, 25, 15 and 13 respectively.
u_trim = zeros(NU, K);
dd_lb  = zeros(NU, K);
dd_ub  = zeros(NU, K);
for k = 1:K
    u_trim(:,k) = [XU0(17:24,k); XU0(25,k); XU0(15,k); XU0(13,k)];
    dd_lb(:,k)  = LIM_LO - u_trim(:,k);
    dd_ub(:,k)  = LIM_HI - u_trim(:,k);
end

%% 3. STEP 1 - largest and smallest steady acceleration --------------------
% Decision vector z = [ddelta(11); dtheta; a]; a enters the equality as a
% free variable so the LP can push it to its extreme.
a_min = nan(K,1);  a_max = nan(K,1);
ddelta_at_amin = nan(K,NU);  ddelta_at_amax = nan(K,NU);
dtheta_at_amin = nan(K,1);   dtheta_at_amax = nan(K,1);
flag_amin = zeros(K,1);      flag_amax = zeros(K,1);

for k = 1:K
    Aeq = [BP(1:3,:,k), AP(1:3,4,k), [-1;0;0]];
    beq = zeros(3,1);
    lb  = [dd_lb(:,k); -Inf; -Inf];
    ub  = [dd_ub(:,k);  Inf;  Inf];

    [z, ~, ef] = linprog([zeros(NU+1,1); -1], [], [], Aeq, beq, lb, ub, opts);
    flag_amax(k) = ef;
    if ef == 1
        a_max(k) = z(end);  dtheta_at_amax(k) = z(NU+1);
        ddelta_at_amax(k,:) = z(1:NU)';
    end

    [z, ~, ef] = linprog([zeros(NU+1,1); +1], [], [], Aeq, beq, lb, ub, opts);
    flag_amin(k) = ef;
    if ef == 1
        a_min(k) = z(end);  dtheta_at_amin(k) = z(NU+1);
        ddelta_at_amin(k,:) = z(1:NU)';
    end
end

%% 4. STEP 2 - acceleration grid per trim point ----------------------------
% Each trim point gets its own grid because the reachable a range differs.
a_grid = nan(K, NA);
for k = 1:K
    if ~isnan(a_min(k)) && ~isnan(a_max(k))
        a_grid(k,:) = linspace(a_min(k), a_max(k), NA);
    end
end

%% 5. STEP 3 + STEP 5 - pitch band and the actuator solution behind it -----
% a is a constant now, so z = [ddelta(11); dtheta] and the LP walks dtheta to
% both extremes. The actuator vector at each extreme is kept so the binding
% actuator can be identified afterwards.
dtheta_min = nan(K,NA);  dtheta_max = nan(K,NA);
ddelta_at_min = nan(K,NA,NU);  ddelta_at_max = nan(K,NA,NU);
active_at_min = zeros(K,NA,NU,'int8');   % -1 at lower bound, +1 at upper, 0 free
active_at_max = zeros(K,NA,NU,'int8');
resid_max = 0;

for k = 1:K
    Aeq = [BP(1:3,:,k), AP(1:3,4,k)];
    lb  = [dd_lb(:,k); -Inf];
    ub  = [dd_ub(:,k);  Inf];
    tolb = 1e-6 * max(1, max(abs([dd_lb(:,k); dd_ub(:,k)])));
    for i = 1:NA
        a = a_grid(k,i);
        if isnan(a), continue; end
        beq = [a; 0; 0];

        [z, ~, ef] = linprog([zeros(NU,1); -1], [], [], Aeq, beq, lb, ub, opts);
        if ef == 1
            dtheta_max(k,i) = z(end);
            ddelta_at_max(k,i,:) = z(1:NU);
            active_at_max(k,i,:) = int8( (z(1:NU) >= dd_ub(:,k)-tolb) ...
                                       - (z(1:NU) <= dd_lb(:,k)+tolb) );
            resid_max = max(resid_max, norm(Aeq*z - beq));
        end

        [z, ~, ef] = linprog([zeros(NU,1); +1], [], [], Aeq, beq, lb, ub, opts);
        if ef == 1
            dtheta_min(k,i) = z(end);
            ddelta_at_min(k,i,:) = z(1:NU);
            active_at_min(k,i,:) = int8( (z(1:NU) >= dd_ub(:,k)-tolb) ...
                                       - (z(1:NU) <= dd_lb(:,k)+tolb) );
            resid_max = max(resid_max, norm(Aeq*z - beq));
        end
    end
end

%% 6. STEP 4 - what the shipped GUAM allocation picks ----------------------
% M = W^-1 B' (B W^-1 B')^-1; its 12th row is the virtual pitch effector, so
% M(12,1)*a is the pitch the current allocation spends on forward accel.
M_theta_row = nan(K,1);
dtheta_alloc = nan(K,NA);
for k = 1:K
    W = WL(:,:,k);  B = BL(:,:,k);
    M = (W \ B') / (B * (W \ B'));       % 12 x 3
    M_theta_row(k) = M(12,1);
    dtheta_alloc(k,:) = M_theta_row(k) * a_grid(k,:);
end

%% 7. probe accelerations - bands at the three schedule accelerations ------
% Solved directly (not interpolated) so the summary and the band plot agree.
NP = numel(A_PROBE);
dth_min_probe = nan(K,NP);  dth_max_probe = nan(K,NP);
dth_alloc_probe = nan(K,NP);  probe_feasible = false(K,NP);
for k = 1:K
    Aeq = [BP(1:3,:,k), AP(1:3,4,k)];
    lb  = [dd_lb(:,k); -Inf];
    ub  = [dd_ub(:,k);  Inf];
    for p = 1:NP
        a = A_PROBE(p);
        if isnan(a_min(k)) || a < a_min(k) || a > a_max(k), continue; end
        beq = [a; 0; 0];
        [z1, ~, e1] = linprog([zeros(NU,1); -1], [], [], Aeq, beq, lb, ub, opts);
        [z2, ~, e2] = linprog([zeros(NU,1); +1], [], [], Aeq, beq, lb, ub, opts);
        if e1 == 1 && e2 == 1
            dth_max_probe(k,p) = z1(end);
            dth_min_probe(k,p) = z2(end);
            dth_alloc_probe(k,p) = M_theta_row(k) * a;
            probe_feasible(k,p) = true;
        end
    end
end

%% 8. CHECK 1 - hover must reduce to dtheta = -a/g -------------------------
% At UH = 0 the pusher makes no thrust, so tilting is the only way to
% accelerate: the band should collapse onto -a/g with no width.
a_chk = 1.0;
Aeq = [BP(1:3,:,1), AP(1:3,4,1)];
lb  = [dd_lb(:,1); -Inf];  ub = [dd_ub(:,1); Inf];
[z1, ~, e1] = linprog([zeros(NU,1); -1], [], [], Aeq, [a_chk;0;0], lb, ub, opts);
[z2, ~, e2] = linprog([zeros(NU,1); +1], [], [], Aeq, [a_chk;0;0], lb, ub, opts);
chk1_expect = -a_chk / G;
chk1_ok = false;
fprintf('--- CHECK 1: hover reduces to dtheta = -a/g ---\n');
if e1 == 1 && e2 == 1
    chk1_hi = z1(end);  chk1_lo = z2(end);
    err = max(abs([chk1_hi chk1_lo] - chk1_expect));
    chk1_ok = err < 2e-3;
    fprintf('  a = %.2f ft/s^2 at UH = %.1f\n', a_chk, UH_used(1));
    fprintf('  expected      : %+.6f rad\n', chk1_expect);
    fprintf('  dtheta_max    : %+.6f rad\n', chk1_hi);
    fprintf('  dtheta_min    : %+.6f rad\n', chk1_lo);
    fprintf('  worst error   : %.2e rad   -> %s\n', err, ternary(chk1_ok,'PASS','FAIL'));
else
    fprintf('  LP infeasible at hover (exitflags %d / %d) -> FAIL\n', e1, e2);
end
if ~chk1_ok
    fprintf(['\n  *** CHECK 1 FAILED. A sign or index error in Ap/Bp is the\n' ...
             '      likely cause. Do NOT trust any other number in this run\n' ...
             '      until this is resolved. ***\n']);
end

%% 9. CHECK 2 - equality residual -----------------------------------------
% Every stored solution must satisfy the equilibrium to solver tolerance.
fprintf('\n--- CHECK 2: equality residual ---\n');
fprintf('  max |Bp*ddelta + Ap(:,4)*dtheta - [a;0;0]| = %.3e   -> %s\n', ...
        resid_max, ternary(resid_max <= 1e-8, 'PASS', 'FAIL (>1e-8)'));

%% 10. CHECK 3 - the shipped allocation must lie inside the band -----------
% If it does not, either the allocation matrix or the box conversion is wrong.
viol = 0;  viol_list = [];
for k = 1:K
    for i = 1:NA
        if isnan(dtheta_min(k,i)) || isnan(dtheta_max(k,i)), continue; end
        tol = 1e-6 + 1e-6*abs(dtheta_alloc(k,i));
        if dtheta_alloc(k,i) < dtheta_min(k,i)-tol || ...
           dtheta_alloc(k,i) > dtheta_max(k,i)+tol
            viol = viol + 1;
            viol_list(end+1,:) = [k, i, dtheta_alloc(k,i), ...
                                  dtheta_min(k,i), dtheta_max(k,i)]; %#ok<SAGROW>
        end
    end
end
fprintf('\n--- CHECK 3: allocation inside [dtheta_min, dtheta_max] ---\n');
fprintf('  violations: %d of %d   -> %s\n', viol, sum(~isnan(dtheta_min(:))), ...
        ternary(viol == 0, 'PASS', 'FAIL'));
if viol > 0
    fprintf('  first few (k, i, alloc, min, max):\n');
    disp(viol_list(1:min(8,end),:));
    fprintf(['  A non-zero count means the allocation matrix or the box\n' ...
             '  conversion is wrong - the band and the allocation disagree.\n']);
end

%% 11. flags - linearity range and the table splice ------------------------
% Nothing is corrected here; the suspect points are only reported.
lin_flag = (abs(dtheta_min) > DTH_LIN_FLAG) | (abs(dtheta_max) > DTH_LIN_FLAG);
fprintf('\n--- FLAGS ---\n');
fprintf('  |dtheta| > %.2f rad (%.1f deg), outside linear validity: %d of %d grid points\n', ...
        DTH_LIN_FLAG, rad2deg(DTH_LIN_FLAG), sum(lin_flag(:)), sum(~isnan(dtheta_min(:))));
kf = find(any(lin_flag,2));
if ~isempty(kf)
    fprintf('    affected trim points: %s\n', mat2str(kf'));
end

% jump detector on a_max: flag any step more than 3x the median step
da = abs(diff(a_max));
jump_thr = 3 * median(da(~isnan(da)));
jump_k = find(da > jump_thr);
fprintf('  a_max jump detector (threshold %.2f ft/s^2 between neighbours):\n', jump_thr);
if isempty(jump_k)
    fprintf('    none\n');
else
    for j = jump_k'
        fprintf('    UH %d -> %d  (%.1f -> %.1f ft/s):  a_max %.2f -> %.2f  (jump %.2f)\n', ...
                j, j+1, UH_used(j), UH_used(j+1), a_max(j), a_max(j+1), da(j));
    end
end
if any(ismember(jump_k, SPLICE_WARN))
    fprintf(['    NOTE: a jump falls in segments %s, near the trim-table\n' ...
             '    splice (UH(20) = %.4f, UH(21) = %.4f). Reported as found,\n' ...
             '    NOT smoothed.\n'], mat2str(SPLICE_WARN), S.UH(20), S.UH(21));
end

%% 12. console summary table ----------------------------------------------
% One row per trim point, with the pitch band at the low-band acceleration.
p_ref = 1;   % A_PROBE(1) = 3.38 ft/s^2
fprintf('\n--- summary (pitch band shown at a = %.2f ft/s^2) ---\n', A_PROBE(p_ref));
fprintf('%3s | %8s | %8s | %8s | %19s | %9s | %s\n', ...
        'k','UH[ft/s]','a_min','a_max','dtheta band [deg]','alloc[deg]','note');
fprintf('%s\n', repmat('-',1,86));
for k = 1:K
    if probe_feasible(k,p_ref)
        band = sprintf('%+7.2f .. %+7.2f', rad2deg(dth_min_probe(k,p_ref)), ...
                                           rad2deg(dth_max_probe(k,p_ref)));
        alc  = sprintf('%+8.2f', rad2deg(dth_alloc_probe(k,p_ref)));
    else
        band = sprintf('%19s', 'a out of range');
        alc  = sprintf('%8s', '-');
    end
    note = '';
    if any(lin_flag(k,:)), note = [note '|dtheta|>0.2 somewhere ']; end %#ok<AGROW>
    if ismember(k, jump_k) || ismember(k-1, jump_k), note = [note 'a_max jump ']; end %#ok<AGROW>
    fprintf('%3d | %8.1f | %8.2f | %8.2f | %19s | %9s | %s\n', ...
            k, UH_used(k), a_min(k), a_max(k), band, alc, note);
end

% which actuator limits the pitch band, at the reference acceleration
fprintf('\n--- binding actuators at the dtheta extremes (a = %.2f) ---\n', A_PROBE(p_ref));
fprintf('%3s | %8s | %-34s | %s\n', 'k','UH','at dtheta_max','at dtheta_min');
fprintf('%s\n', repmat('-',1,86));
for k = 1:K
    [~, ia] = min(abs(a_grid(k,:) - A_PROBE(p_ref)));
    if isnan(dtheta_max(k,ia)), continue; end
    fprintf('%3d | %8.1f | %-34s | %s\n', k, UH_used(k), ...
            active_str(squeeze(active_at_max(k,ia,:)), ACT_NAME), ...
            active_str(squeeze(active_at_min(k,ia,:)), ACT_NAME));
end

%% 13. figure (a) - acceleration envelope ---------------------------------
% What the actuators can hold, before any reachable-set constraint.
f1 = figure('Position',[80 80 900 480]);
fill([UH_used; flipud(UH_used)], [a_max; flipud(a_min)], [0.85 0.90 0.98], ...
     'EdgeColor','none'); hold on;
plot(UH_used, a_max, 'b-o', 'LineWidth', 1.4, 'MarkerSize', 4);
plot(UH_used, a_min, 'r-o', 'LineWidth', 1.4, 'MarkerSize', 4);
yline(0, 'k-');
for p = 1:NP, yline(A_PROBE(p), 'k--', sprintf('a = %.2f', A_PROBE(p))); end
xlabel('trim speed UH [ft/s]'); ylabel('steady forward acceleration [ft/s^2]');
title('steady-acceleration envelope at WH3 (actuator limits only, no BRT)');
legend({'reachable','a_{max}','a_{min}'}, 'Location','best'); grid on;
saveas(f1, fullfile(root,'logger','accel_envelope.png'));

%% 14. figure (b) - pitch band at the three probe accelerations ------------
% Band = what the actuators allow; line = what the shipped allocation picks.
f2 = figure('Position',[80 80 980 760]);
for p = 1:NP
    subplot(NP,1,p); hold on;
    ok = probe_feasible(:,p);
    if any(ok)
        uu = UH_used(ok);
        lo = rad2deg(dth_min_probe(ok,p));  hi = rad2deg(dth_max_probe(ok,p));
        fill([uu; flipud(uu)], [hi; flipud(lo)], [0.88 0.94 0.88], 'EdgeColor','none');
        plot(uu, hi, 'g-',  'LineWidth', 1.1);
        plot(uu, lo, 'g-',  'LineWidth', 1.1);
        plot(uu, rad2deg(dth_alloc_probe(ok,p)), 'm-s', 'LineWidth', 1.5, 'MarkerSize', 4);
    end
    yline(0,'k-');
    yline( rad2deg(DTH_LIN_FLAG), 'r:', 'linearity limit');
    yline(-rad2deg(DTH_LIN_FLAG), 'r:');
    ylabel('d\theta [deg]'); grid on;
    title(sprintf('a = %.2f ft/s^2   (%d of %d trim points feasible)', ...
                  A_PROBE(p), sum(ok), K));
    if p == 1, legend({'allowed band','','','GUAM allocation'}, 'Location','best'); end
    if p == NP, xlabel('trim speed UH [ft/s]'); end
end
saveas(f2, fullfile(root,'logger','accel_envelope_bands.png'));

%% 15. save ---------------------------------------------------------------
% Everything a later BRT-constrained pass needs to look up.
active_constraint = struct( ...
    'at_dtheta_max', active_at_max, ...
    'at_dtheta_min', active_at_min, ...
    'actuator_names', {ACT_NAME}, ...
    'code_meaning', '-1 = at lower limit, +1 = at upper limit, 0 = interior');

meta = struct('WH_IDX', WH_IDX, 'wbar', S.WH(WH_IDX), 'g', G, ...
              'lim_lo', LIM_LO, 'lim_hi', LIM_HI, 'a_probe', A_PROBE, ...
              'dtheta_linear_flag', DTH_LIN_FLAG, ...
              'trim_table', which('trim_table_Poly_ConcatVer4p0.mat'));

out = fullfile(root,'data','accel_trim_envelope_WH3.mat');
save(out, 'UH_used','a_min','a_max','a_grid','dtheta_min','dtheta_max', ...
     'dtheta_alloc','ddelta_at_max','ddelta_at_min','active_constraint', ...
     'ddelta_at_amax','ddelta_at_amin','dtheta_at_amax','dtheta_at_amin', ...
     'M_theta_row','u_trim','dd_lb','dd_ub','lin_flag', ...
     'dth_min_probe','dth_max_probe','dth_alloc_probe','probe_feasible','meta');

fprintf('\nsaved %s\n', out);
fprintf('saved logger/accel_envelope.png, logger/accel_envelope_bands.png\n');

%% ------------------------------------------------------------------------
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function s = active_str(code, names)
hi = find(code > 0);  lo = find(code < 0);
parts = {};
for j = hi(:)', parts{end+1} = [names{j} '(hi)']; end %#ok<AGROW>
for j = lo(:)', parts{end+1} = [names{j} '(lo)']; end %#ok<AGROW>
if isempty(parts), s = '(none binding)'; else, s = strjoin(parts, ', '); end
end
