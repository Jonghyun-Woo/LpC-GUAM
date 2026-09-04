% mu_theta_nonlinear - is the linear pitch-to-acceleration relation optimistic
% or conservative once the pitch gets large?
%
% The transition calculation uses  dtheta = mu_theta(u) * a  with mu_theta
% taken from the allocation matrix, i.e. from a small perturbation about the
% trim point. The accelerations that calculation wants demand 13-19 deg of
% pitch. This script asks only for the SIGN of the error there, not for a
% replacement model.
%
% PART A  the hover relation, exactly, with no aero model at all
% PART B  nonlinear force balance from the polynomial aero/propulsion model
% PART C  what the measured validity range does to the schedule
%
% ---------------------------------------------------------------------------
% ONE DEPARTURE FROM THE BRIEF, STATED UP FRONT
% ---------------------------------------------------------------------------
% B3 asks to restore vertical equilibrium with  Z + m*g*cos(theta) = 0. That
% is the BODY-axis condition (w_body held), and it contradicts two other parts
% of the same brief: B3's own text says ubar and wbar are held, and CHECK 3
% demands agreement with a = -g*tan(theta).
%
%   body-axis   wdot_b = 0        ->  T = m*g*cos(theta),  a = -g*sin(theta)
%   heading     wdot_bar = 0      ->  T = m*g/cos(theta),  a = -g*tan(theta)
%
% Only the second one holds the descent rate fixed, which is what "ubar, wbar
% fixed" means and what the commanded speed (= ubar) is measured against. So
% the heading-frame condition is used as the primary result, and the literal
% body-axis condition is computed alongside as a cross-check. Both are saved.
%
% Reads only. Output: data/mu_theta_nonlinear.mat, logger/mu_theta_nl_*.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

K       = 20;
WH_IDX  = 3;
DTH     = (-0.35:0.025:0)';      % 15 pitch perturbations [rad]
TOL     = 0.10;                  % linear valid while |ratio-1| <= TOL
W_BODY_LIM = 17;                 % ft/s, roughly +-10 kt - model-range warning
S_LO = 0.0;  S_HI = 4.0;         % rotor scale bracket

line = @(c) fprintf('%s\n', repmat(c,1,100));

%% ========================================================================
%  0. is the aero model actually here?
%% ========================================================================
need = {'poly_aero_wrapper_Mod_du','LpC_aero_p_v2','LpC_interp_p_v2', ...
        'LpC_model_parameters','check_fac_limits'};
missing = need(cellfun(@(f) isempty(which(f)), need));
HAVE_MODEL = isempty(missing);
line('='); fprintf('0. AERO MODEL AVAILABILITY\n'); line('=');
for i = 1:numel(need)
    fprintf('  %-28s %s\n', need{i}, ternstr(isempty(which(need{i})),'MISSING','found'));
end
if ~HAVE_MODEL
    fprintf('\n  Missing: %s\n', strjoin(missing,', '));
    fprintf('  PART A will run; PART B and PART C cannot and will be skipped.\n');
end

cfgV = VehicleConfig();  m = cfgV.mass;
cfgR = RSLQRConfig();    rho = cfgR.rho;  G = cfgR.grav;
Wt = m*G;
fprintf('\n  mass %.3f slug, weight %.1f lbf, rho %.6f slug/ft^3\n', m, Wt, rho);

%% ========================================================================
%  PART A - hover, exactly, without any model
%% ========================================================================
line('='); fprintf('PART A - exact vs linear pitch at hover\n'); line('=');
aA = (1:12)';
th_exact  = -atan(aA/G);
th_linear = -aA/G;
ratioA    = th_exact ./ th_linear;

fprintf('%6s | %14s | %14s | %8s\n','a','theta_exact[d]','theta_lin[d]','ratio');
line('-');
for i = 1:numel(aA)
    fprintf('%6.1f | %14.3f | %14.3f | %8.4f\n', ...
            aA(i), rad2deg(th_exact(i)), rad2deg(th_linear(i)), ratioA(i));
end
fprintf(['\nratio < 1 everywhere: the exact relation needs LESS pitch than the\n' ...
         'linear one for the same acceleration, so at hover the linear model is\n' ...
         'CONSERVATIVE - it under-reports what the vehicle can do.\n' ...
         'At a = 12 the linear model over-states the required pitch by %.1f %%.\n'], ...
        100*(1/ratioA(end) - 1));

if ~HAVE_MODEL
    save(fullfile(dataDir,'mu_theta_nonlinear.mat'), 'aA','th_exact','th_linear','ratioA');
    fprintf('\nSTOPPED after PART A - the aero model is not available.\n');
    return;
end

%% ========================================================================
%  PART B setup
%% ========================================================================
T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:K);  WH3 = T.WH(WH_IDX);
XU0 = T.XU0_interp(:, 1:K, WH_IDX);
th_trim = XU0(11,:)';

mu_theta = zeros(K,1);
for k = 1:K
    W = T.W_lon_interp(:,:,k,WH_IDX);  Bl = T.B_lon_interp(:,:,k,WH_IDX);
    M = (W \ Bl') / (Bl * (W \ Bl'));
    mu_theta(k) = M(12,1);
end

CM = LpC_model_parameters();
Units.deg = pi/180;  Units.knot = 1.6878;  a_snd = 1125.33;  blend = 2;
% force(vb, del, omp) -> [X;Y;Z] and the model's own validity flags
forceof = @(vb, del, omp) call_aero(vb, del, omp, rho, a_snd, blend, Units, CM);

%% ========================================================================
%  PART B2 - CHECK 1: does the stored trim actually balance?
%% ========================================================================
line('='); fprintf('PART B2 / CHECK 1 - residual of the stored trim points\n'); line('=');
fprintf('%3s | %8s | %10s %10s | %10s %10s | %8s | %s\n', ...
        'k','UH','u_body','w_body','resX[lbf]','resZ[lbf]','%% of W','validity');
line('-');
resX = nan(K,1); resZ = nan(K,1); ubody = nan(K,1); wbody = nan(K,1);
vld0 = cell(K,1);
for k = 1:K
    th = th_trim(k);
    ubody(k) = UH(k)*cos(th) - WH3*sin(th);
    wbody(k) = UH(k)*sin(th) + WH3*cos(th);
    [Fk, vk] = forceof([ubody(k);0;wbody(k)], XU0(13:16,k), XU0(17:25,k));
    resX(k) = Fk(1) - Wt*sin(th);
    resZ(k) = Fk(3) + Wt*cos(th);
    vld0{k} = vstr(vk);
    fprintf('%3d | %8.1f | %10.3f %10.3f | %10.3f %10.3f | %7.3f%% | %s\n', ...
            k, UH(k), ubody(k), wbody(k), resX(k), resZ(k), ...
            100*max(abs([resX(k) resZ(k)]))/Wt, vld0{k});
end
worst = max(max(abs([resX resZ])))/Wt;
fprintf('\n  worst residual %.3f %% of weight -> %s\n', 100*worst, ...
        ternstr(worst < 0.01, 'PASS (input convention is right)', 'FAIL'));
if worst >= 0.01
    fprintf(['\n  *** CHECK 1 FAILED. The input convention is wrong; every number\n' ...
             '      below would be meaningless. Stopping. ***\n']);
    return;
end
nflag = sum(abs(wbody) > W_BODY_LIM);
fprintf('  trim points with |w_body| > %.0f ft/s (model-range warning): %d\n', ...
        W_BODY_LIM, nflag);

%% ========================================================================
%  PART B3-B4 - tilt, restore vertical equilibrium, read the acceleration
%  Rotor scaling: all EIGHT lift rotors multiplied by one common factor s,
%  pusher and surfaces left at their trim values. s is found by bisection.
%% ========================================================================
ND = numel(DTH);
a_nl   = nan(K,ND);   % heading-frame condition (primary)
a_nl_b = nan(K,ND);   % body-frame condition (the literal B3 wording)
a_li   = nan(K,ND);
s_used = nan(K,ND);   rotor_max = nan(K,ND);
oor    = false(K,ND); % model reported the point out of range

for k = 1:K
    del = XU0(13:16,k);  omp0 = XU0(17:25,k);
    for i = 1:ND
        thn = th_trim(k) + DTH(i);
        ub = UH(k)*cos(thn) - WH3*sin(thn);
        wb = UH(k)*sin(thn) + WH3*cos(thn);

        % heading-frame: wdot_bar = 0  <=>  wdot_b = udot_b * tan(theta)
        gh = @(s) resid_head(s, ub, wb, thn, del, omp0, forceof, m, G);
        s1 = bisect(gh, S_LO, S_HI);
        if ~isnan(s1)
            omp = omp0; omp(1:8) = omp0(1:8)*s1;
            [Fk, vk] = forceof([ub;0;wb], del, omp);
            udot_b = Fk(1)/m - G*sin(thn);
            a_nl(k,i) = udot_b / cos(thn);          % = d(ubar)/dt
            s_used(k,i) = s1;  rotor_max(k,i) = max(omp(1:8));
            oor(k,i) = any(vflags(vk));
        end

        % body-frame cross-check: Z + m*g*cos(theta) = 0
        gb = @(s) resid_body(s, ub, wb, thn, del, omp0, forceof, m, G);
        s2 = bisect(gb, S_LO, S_HI);
        if ~isnan(s2)
            omp = omp0; omp(1:8) = omp0(1:8)*s2;
            Fk = forceof([ub;0;wb], del, omp);
            a_nl_b(k,i) = Fk(1)/m - G*sin(thn);
        end

        a_li(k,i) = DTH(i) / mu_theta(k);
    end
end
ratio_mu = a_nl ./ a_li;          % against mu_theta, as the brief specifies

% ---------------------------------------------------------------------------
% Second comparator, added after ratio_mu turned out not to approach 1.
% mu_theta is the share of the job the ALLOCATION assigns to pitch, with the
% pusher taking the rest; the sweep above tilts the aircraft with the pusher
% held at trim. Away from hover those are different questions, so ratio_mu
% carries a constant offset that has nothing to do with large-angle behaviour.
% The question actually asked - is the curve super- or sub-linear - is answered
% by comparing the curve with ITS OWN tangent at the origin.
% ---------------------------------------------------------------------------
slope0 = nan(K,1);
for k = 1:K
    if ~isnan(a_nl(k,ND-1)), slope0(k) = a_nl(k,ND-1) / DTH(ND-1); end
end
a_tan = slope0 .* DTH';           % K x ND
ratio = a_nl ./ a_tan;            % 1 at the origin by construction
mu_implied = 1 ./ slope0;         % what the model says mu_theta should be

%% ========================================================================
%  CHECK 2 and CHECK 3
%% ========================================================================
line('='); fprintf('CHECK 2 - does the model agree with mu_theta as dtheta -> 0?\n'); line('=');
i_small = ND-1;                                    % dtheta = -0.025
fprintf('  at dtheta = %+.3f rad (pusher held at trim):\n', DTH(i_small));
fprintf('%3s | %8s | %12s | %12s | %10s | %12s | %10s\n', ...
        'k','UH','a_nonlinear','a_linear','ratio_mu','mu_implied','mu_theta');
line('-');
for k = 1:K
    fprintf('%3d | %8.1f | %12.4f | %12.4f | %10.4f | %12.6f | %10.6f\n', ...
            k, UH(k), a_nl(k,i_small), a_li(k,i_small), ratio_mu(k,i_small), ...
            mu_implied(k), mu_theta(k));
end
r0 = ratio_mu(:,i_small);
fprintf('\n  ratio_mu range %.4f .. %.4f  -> %s\n', min(r0), max(r0), ...
        ternstr(all(abs(r0-1) < 0.10), 'PASS', 'FAIL'));
fprintf(['\n  DIAGNOSIS. The gap is not a large-angle effect: it is already there at\n' ...
         '  the smallest step. mu_theta is the pitch share the ALLOCATION assigns\n' ...
         '  while the pusher takes the remainder; this sweep tilts with the pusher\n' ...
         '  frozen, so it measures pitch acting alone. The two agree at hover\n' ...
         '  (ratio %.4f, where the pusher makes no thrust) and separate as the\n' ...
         '  pusher takes over. mu_implied above is what the model alone says.\n' ...
         '  Because of this, ratio_mu cannot answer the large-angle question, and\n' ...
         '  everything below uses the curve against its own tangent instead.\n'], ...
        ratio_mu(1,i_small));

line('='); fprintf('CHECK 3 - hover against the exact tan relation\n'); line('=');
% Referenced to the trim attitude: the trim point itself has a = 0, so the
% exact hover relation to compare against is the CHANGE in -g*tan.
a_tan_exact = -G*(tan(th_trim(1)+DTH) - tan(th_trim(1)));
fprintf('%12s | %14s | %18s | %10s\n','dtheta[rad]','a_nonlinear','-g*dtan(theta)','diff');
line('-');
for i = 1:ND
    fprintf('%12.3f | %14.4f | %18.4f | %10.4f\n', ...
            DTH(i), a_nl(1,i), a_tan_exact(i), a_nl(1,i) - a_tan_exact(i));
end
dmax = max(abs(a_nl(1,:)' - a_tan_exact));
fprintf('\n  worst difference %.4f ft/s^2 -> %s\n', dmax, ...
        ternstr(dmax < 0.5, 'PASS', 'differs - see note'));
fprintf(['  note: WH3 is a %.2f ft/s descent, so hover here is not truly zero\n' ...
         '  airspeed. A small aero contribution is expected and is not an error.\n'], WH3);

%% ========================================================================
%  PART B5 - how far the linear relation stays within TOL
%% ========================================================================
dth_valid = zeros(K,1);
for k = 1:K
    dth_valid(k) = 0;
    for i = ND-1:-1:1                       % walk outward from dtheta = 0
        if isnan(ratio(k,i)) || abs(ratio(k,i)-1) > TOL, break; end
        dth_valid(k) = abs(DTH(i));
    end
end

line('='); fprintf('TABLE 2 - curvature: the curve against its own tangent (|ratio-1| <= %.0f %%)\n', 100*TOL); line('=');
fprintf('%3s | %8s | %12s | %12s %10s | %12s | %s\n', ...
        'k','UH','mu_implied','dth_valid[r]','[deg]','ratio@-0.2','notes');
line('-');
i20 = find(abs(DTH + 0.20) < 1e-9);
for k = 1:K
    nt = '';
    if any(oor(k,:)), nt = [nt sprintf('out-of-range at %d pts ', sum(oor(k,:)))]; end %#ok<AGROW>
    if abs(wbody(k)) > W_BODY_LIM, nt = [nt 'w_body high ']; end %#ok<AGROW>
    fprintf('%3d | %8.1f | %12.6f | %12.3f %10.1f | %12.4f | %s\n', ...
            k, UH(k), mu_implied(k), dth_valid(k), rad2deg(dth_valid(k)), ...
            ratio(k,i20), nt);
end

nopt = sum(ratio(:,i20) > 1);  npes = sum(ratio(:,i20) < 1);
fprintf(['\nSIGN VERDICT at dtheta = -0.20 rad, curve vs its own tangent:\n' ...
         '  ratio > 1 at %d of %d trim points, ratio < 1 at %d.\n' ...
         '  ratio > 1 means the straight line UNDER-states what the tilt delivers,\n' ...
         '  i.e. the linear relation is CONSERVATIVE at that pitch.\n'], nopt, K, npes);

%% ========================================================================
%  PART C - what that does to the schedule
%% ========================================================================
line('='); fprintf('PART C - schedule under the measured validity range\n'); line('=');
prev = fullfile(dataDir,'min_transition_time_WH3.mat');
if ~isfile(prev)
    fprintf('  data/min_transition_time_WH3.mat not found - PART C skipped.\n');
    a_star = [];
else
    P = load(prev);  a_star = P.a_star;  du_seg = P.du_seg;  NS = numel(a_star);
    a_lin_arb = 0.20 ./ abs(mu_theta(1:NS));       % the earlier arbitrary 0.2 rad
    a_lin_new = dth_valid(1:NS) ./ abs(mu_theta(1:NS));
    a_safe_arb = min(a_star, a_lin_arb);
    a_safe_new = min(a_star, a_lin_new);
    T_arb = du_seg ./ a_safe_arb;   T_arb(a_safe_arb <= 0) = Inf;
    T_new = du_seg ./ a_safe_new;   T_new(a_safe_new <= 0) = Inf;
    if any(~isfinite(T_new))
        fprintf(['  ** dth_valid = 0 on segments %s: the curvature already exceeds\n' ...
                 '     %.0f %% at the first step there, so no positive linear range\n' ...
                 '     exists and T is infinite. Reported, not patched.\n\n'], ...
                mat2str(find(~isfinite(T_new))'), 100*TOL);
    end

    fprintf('%3s | %9s | %13s | %11s | %12s | %9s\n', ...
            'k','a_star','a_lin(0.2 arb)','a_lin_new','a_safe_new','T_new [s]');
    line('-');
    for k = 1:NS
        fprintf('%3d | %9.2f | %13.2f | %11.2f | %12.2f | %9.2f\n', ...
                k, a_star(k), a_lin_arb(k), a_lin_new(k), a_safe_new(k), T_new(k));
    end
    fprintf('%3s | %9s | %13s | %11s | %12s | %9.2f\n','sum','','','','', sum(T_new));
    fprintf('\n  T_total with the arbitrary 0.2 rad cap : %.2f s\n', sum(T_arb));
    fprintf('  T_total with the measured range        : %.2f s\n', sum(T_new));
    fprintf(['  Both are eta = 1 and carry no actuator limit. The second is not a\n' ...
             '  target being aimed at - it is whatever the measured range gives.\n']);
end

%% ========================================================================
%  FIGURES
%% ========================================================================
show = [1 6 11 16 20];  show = show(show <= K);
f1 = figure('Position',[60 60 1000 620]); hold on;
cm = lines(numel(show));
for j = 1:numel(show)
    k = show(j);
    plot(DTH, a_li(k,:), '--', 'Color', cm(j,:), 'LineWidth', 1.2);
    plot(DTH, a_nl(k,:), '-o', 'Color', cm(j,:), 'LineWidth', 1.6, 'MarkerSize', 4);
end
xlabel('\Delta\theta [rad]'); ylabel('forward acceleration [ft/s^2]');
title('dashed = linear (\mu_\theta), solid = nonlinear force balance');
legend(reshape([arrayfun(@(k) sprintf('UH%d lin',k), show, 'UniformOutput',false); ...
                arrayfun(@(k) sprintf('UH%d nl', k), show, 'UniformOutput',false)],1,[]), ...
       'Location','northeast','FontSize',8);
grid on;
saveas(f1, fullfile(root,'logger','mu_theta_nl_curves.png'));

f2 = figure('Position',[60 60 900 480]);
plot(UH, rad2deg(dth_valid), 'b-o', 'LineWidth', 1.6); hold on;
yline(rad2deg(0.20), 'r--', 'the 0.20 rad cap used earlier', 'LineWidth', 1.3);
xlabel('trim speed UH [ft/s]'); ylabel('valid |\Delta\theta| [deg]');
title(sprintf('range over which the linear relation stays within %.0f %%', 100*TOL));
grid on;
saveas(f2, fullfile(root,'logger','mu_theta_nl_valid.png'));

%% ========================================================================
%  SAVE
%% ========================================================================
meta = struct('WH_IDX', WH_IDX, 'WH3', WH3, 'UH', UH, 'dtheta', DTH, 'tol', TOL, ...
              'rotor_scaling', 'all 8 lift rotors by one common factor, bisection', ...
              'vertical_condition_primary', 'heading frame: wdot_bar = 0', ...
              'vertical_condition_crosscheck', 'body frame: Z + m*g*cos(theta) = 0', ...
              'note', 'no flight data or previous schedule was used');
save(fullfile(dataDir,'mu_theta_nonlinear.mat'), ...
     'aA','th_exact','th_linear','ratioA','UH','th_trim','mu_theta','DTH', ...
     'a_nl','a_nl_b','a_li','ratio','ratio_mu','slope0','mu_implied','dth_valid','s_used','rotor_max','oor', ...
     'resX','resZ','ubody','wbody','a_star','meta');
fprintf('\nsaved data/mu_theta_nonlinear.mat\n');
fprintf('saved logger/mu_theta_nl_{curves,valid}.png\n');

%% ------------------------------------------------------------------------
function [F, V] = call_aero(vb, del, omp, rho, a_snd, blend, Units, CM)
srf = [del(1)+del(2), del(1)-del(2), del(3), del(3), del(4)];
[X,Y,Z,~,~,~,V] = LpC_aero_p_v2([vb(:); 0; 0; 0]', omp(:)', srf, rho, a_snd, blend, Units, CM);
F = [X; Y; Z];
end

function r = resid_head(s, ub, wb, thn, del, omp0, forceof, m, G)
omp = omp0; omp(1:8) = omp0(1:8)*s;
F = forceof([ub;0;wb], del, omp);
udot_b = F(1)/m - G*sin(thn);
wdot_b = F(3)/m + G*cos(thn);
r = wdot_b - udot_b*tan(thn);        % zero  <=>  d(wbar)/dt = 0
end

function r = resid_body(s, ub, wb, thn, del, omp0, forceof, m, G)
omp = omp0; omp(1:8) = omp0(1:8)*s;
F = forceof([ub;0;wb], del, omp);
r = F(3)/m + G*cos(thn);             % zero  <=>  d(w_body)/dt = 0
end

function s = bisect(g, lo, hi)
a = g(lo); b = g(hi);
if isnan(a) || isnan(b) || a*b > 0, s = NaN; return; end
for i = 1:60
    m = 0.5*(lo+hi);  v = g(m);
    if isnan(v), s = NaN; return; end
    if a*v <= 0, hi = m; b = v; else, lo = m; a = v; end
end
s = 0.5*(lo+hi);
end

function s = vstr(V)
% LpC_interp_p_v2 returns Validity = [invalid_speed invalid_prop_speed].
v = vflags(V);
nm = {'invalid_speed','invalid_prop_speed'};
on = nm(v(1:min(numel(v),2)) ~= 0);
if isempty(on), s = 'ok'; else, s = strjoin(on, ','); end
end

function v = vflags(V)
if isstruct(V)
    f = fieldnames(V);  v = false(1,numel(f));
    for i = 1:numel(f), v(i) = any(V.(f{i})(:)); end
else
    v = logical(V(:))';
end
end

function s = ternstr(c,a,b), if c, s = a; else, s = b; end, end
