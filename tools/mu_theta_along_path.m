% mu_theta_along_path - sweep the acceleration along the manoeuvre the vehicle
% actually flies, and check every evaluation point against the aero model's
% own validity ranges.
%
% WHY THE PREVIOUS SWEEP WAS WRONG
% --------------------------------
% tools/mu_theta_nonlinear.m held ubar and wbar fixed and moved theta. With the
% flight-path angle pinned, angle of attack then moves one-for-one with pitch:
% at cruise a 11.5 deg nose-down took alpha from +8.8 deg to -2.6 deg and the
% wing unloaded completely. The vehicle never visits that state. In a real
% acceleration the speed rises and the trim attitude rises with it, cancelling
% part of that alpha swing. This script follows the segment instead.
%
% WHAT IS NEW HERE
%   theta(u) = theta_trim(u) + mu_theta(u)*a      <- theta_trim(u) varies
%   trim effectors interpolated along the segment as well
%   every evaluation point audited against the model branch's own limits
%
% VERTICAL CONDITION - as written in the brief (A4): Z + m*g*cos(theta) = 0,
% the BODY-axis condition, w_body held. The acceleration that falls out is the
% body-axis one, and the conversion to the horizontal plane is
%     a_horizontal = a_body*cos(theta) + wdot_b*sin(theta)
%                  = a_body*cos(theta)          (wdot_b = 0 under this condition)
% a_horizontal is the one used everywhere below, because the commanded speed is
% ubar. Note this condition does NOT hold the descent rate: it gives
% wdot_bar = -a_body*sin(theta). The heading-frame alternative is computed
% alongside as a_horiz_hd so the two can be compared.
%
% No flight log, schedule or previous total is read anywhere in this file.
%
% Reads only. Output: data/mu_theta_along_path.mat, logger/path_sweep_*.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

K       = 20;  NS = K-1;
WH_IDX  = 3;
NU      = 21;                    % samples along a segment
A_SWEEP = 0.5:0.25:12;           % ft/s^2
ACH_TOL = 0.90;                  % a_achieved must reach this fraction of a
KT      = 1.6878;                % ft/s per knot
RPM     = 60/(2*pi);             % rad/s -> rpm
ROT_ACT_MAX = 1600;              % actuator limit [rpm], = RSLQRConfig eng_max
S_LO = 0.05; S_HI = 3.0;
A_REF = 3.0;

line = @(c) fprintf('%s\n', repmat(c,1,104));

%% ========================================================================
%  0. model, trim table, mu_theta, and the model's own validity ranges
%% ========================================================================
cfgV = VehicleConfig();  m = cfgV.mass;
cfgR = RSLQRConfig();    rho = cfgR.rho;  G = cfgR.grav;
CM = LpC_model_parameters();
Units.deg = pi/180;  Units.knot = KT;  a_snd = 1125.33;  blend = 2;

T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:K);  WH3 = T.WH(WH_IDX);
XU0 = T.XU0_interp(:, 1:K, WH_IDX);
th_trim = XU0(11,:)';
du_seg  = diff(UH);

mu_theta = zeros(K,1);
for k = 1:K
    W = T.W_lon_interp(:,:,k,WH_IDX);  Bl = T.B_lon_interp(:,:,k,WH_IDX);
    M = (W \ Bl') / (Bl * (W \ Bl'));
    mu_theta(k) = M(12,1);
end

xlsxPath = fullfile(root,'GUAM','AeroPolynomial','LpC_Polynomial_Models_v2p1', ...
                    'LpC_p_Model_Ranges_v2p1.xlsx');
line('='); fprintf('0. MODEL VALIDITY RANGES\n'); line('=');
fprintf('  range spreadsheet : %s\n', ternstr(isfile(xlsxPath), 'present', 'NOT FOUND'));
fprintf(['  The ranges used below are read from the check_fac_limits calls in\n' ...
         '  LpC_interp_p_v2.m, which is what the model actually enforces:\n' ...
         '    u        [-5, 45] kt hover branch,  [-5, 95] kt transition/cruise\n' ...
         '    v, w     [-10, 10] kt in EVERY branch\n' ...
         '    surfaces [-30, 30] deg\n' ...
         '    lift rotors [550, 1550] rpm,  pusher [750, 1750] rpm\n' ...
         '  A clamped variable has ZERO sensitivity in that direction, so any\n' ...
         '  force change attributed to it there is spurious.\n\n']);

wb_trim = UH.*sin(th_trim) + WH3*cos(th_trim);
fprintf('  w at the trim points, BEFORE any tilt is applied:\n');
fprintf('%5s | %8s | %14s | %11s | %s\n','k','UH','w_body [ft/s]','w_body [kt]','status');
line('-');
for k = 1:K
    fprintf('%5d | %8.1f | %14.3f | %11.2f | %s\n', k, UH(k), wb_trim(k), ...
            wb_trim(k)/KT, ternstr(abs(wb_trim(k)/KT) > 10, '*** CLAMPED ***', 'inside'));
end
n_cl = sum(abs(wb_trim/KT) > 10);
fprintf('\n  %d of %d trim points are already clamped in w at zero acceleration.\n', n_cl, K);
if n_cl > 0
    fprintf(['  Those trim points sit outside the range the aero model was identified\n' ...
             '  over. Nothing downstream can repair that.\n']);
end
fprintf('\n');

evalseg = @(k,a) sweep_segment(k, a, UH, th_trim, mu_theta, XU0, WH3, NU, ...
                               m, G, rho, a_snd, blend, Units, CM, KT, RPM, S_LO, S_HI);

%% ========================================================================
%  CHECK 1 - at a = 0 the path is the trim path, so the bias must be ~0
%% ========================================================================
line('='); fprintf('CHECK 1 - residual acceleration at a = 0 (must be ~0)\n'); line('=');
a_bias = nan(NS,1);
for k = 1:NS
    R0 = evalseg(k, 0);
    a_bias(k) = mean(R0.a_horiz);
end
fprintf('%3s | %8s | %14s      %3s | %8s | %14s\n','k','UH','a at a=0','k','UH','a at a=0');
line('-');
for k = 1:ceil(NS/2)
    fprintf('%3d | %8.1f | %14.4f', k, UH(k), a_bias(k));
    k2 = k + ceil(NS/2);
    if k2 <= NS, fprintf('      %3d | %8.1f | %14.4f', k2, UH(k2), a_bias(k2)); end
    fprintf('\n');
end
fprintf('\n  worst |bias| = %.4f ft/s^2 -> %s\n', max(abs(a_bias)), ...
        ternstr(max(abs(a_bias)) < 0.25, 'PASS', ...
                'nonzero: interpolated trim between anchors is not an exact trim'));
fprintf(['  The bias is subtracted from a_achieved below, so the ratio measures the\n' ...
         '  response to the commanded acceleration and not the interpolation error.\n']);

%% ========================================================================
%  PART A / TABLE 1 - what a demanded acceleration actually delivers
%% ========================================================================
line('='); fprintf('TABLE 1 - demanded vs delivered at a = %.2f ft/s^2\n', A_REF); line('=');
fprintf('%3s | %8s | %12s | %8s | %10s | %10s | %s\n', ...
        'k','UH','a_achieved','ratio','a_worst','clamped %','variables clamped');
line('-');
ach1 = nan(NS,1); rat1 = nan(NS,1); wst1 = nan(NS,1); vio1 = nan(NS,1);
for k = 1:NS
    R = evalseg(k, A_REF);
    ach1(k) = mean(R.a_horiz) - a_bias(k);
    wst1(k) = min(R.a_horiz)  - a_bias(k);
    rat1(k) = ach1(k)/A_REF;
    vio1(k) = 100*mean(R.anyclamp);
    fprintf('%3d | %8.1f | %12.4f | %8.4f | %10.4f | %9.1f%% | %s\n', ...
            k, UH(k), ach1(k), rat1(k), wst1(k), vio1(k), clampstr(R.clampcount, NU));
end

%% ========================================================================
%  PART B - the largest acceleration that is delivered AND legal
%% ========================================================================
% Feasibility is NOT monotone in a. Raising a pitches the nose down, which
% lowers w_body, which can pull a trim point that starts OUTSIDE the model's
% w range back INSIDE it. So the whole sweep is evaluated and the maximum
% feasible a taken; stopping at the first failure would report 0 for every
% segment whose trim point is already clamped.
NA = numel(A_SWEEP);
a_feas = zeros(NS,1);  fail_why = cell(NS,1);
feas_tab = false(NS,NA);  c1_tab = false(NS,NA);
c2_tab = false(NS,NA);    c3_tab = false(NS,NA);
a_feas_lo = nan(NS,1);    feas_contig = true(NS,1);
for k = 1:NS
    for ia = 1:NA
        a = A_SWEEP(ia);
        R = evalseg(k, a);
        ach = mean(R.a_horiz) - a_bias(k);
        c1_tab(k,ia) = ach >= ACH_TOL*a;
        c2_tab(k,ia) = ~any(R.anyclamp);
        c3_tab(k,ia) = all(R.rot_rpm_max <= ROT_ACT_MAX) && all(isfinite(R.scale));
        feas_tab(k,ia) = c1_tab(k,ia) && c2_tab(k,ia) && c3_tab(k,ia);
    end
    ok = feas_tab(k,:);
    if any(ok)
        i2 = find(ok, 1, 'last');  i1 = find(ok, 1, 'first');
        a_feas(k) = A_SWEEP(i2);  a_feas_lo(k) = A_SWEEP(i1);
        feas_contig(k) = all(ok(i1:i2));
        if i2 == NA
            fail_why{k} = 'sweep ceiling';
        else
            why = {};
            if ~c1_tab(k,i2+1), why{end+1} = 'delivery < 90%'; end %#ok<AGROW>
            if ~c2_tab(k,i2+1), why{end+1} = 'model range'; end %#ok<AGROW>
            if ~c3_tab(k,i2+1), why{end+1} = 'rotor limit'; end %#ok<AGROW>
            fail_why{k} = strjoin(why, ' + ');
        end
    else
        why = {};
        if ~any(c1_tab(k,:)), why{end+1} = 'delivery < 90% everywhere'; end %#ok<AGROW>
        if ~any(c2_tab(k,:)), why{end+1} = 'model range at every a'; end %#ok<AGROW>
        if ~any(c3_tab(k,:)), why{end+1} = 'rotor limit at every a'; end %#ok<AGROW>
        if isempty(why), why = {'no a satisfies all three at once'}; end
        fail_why{k} = strjoin(why, ' + ');
    end
end
if any(a_feas > 0 & a_feas_lo > A_SWEEP(1))
    kk = find(a_feas > 0 & a_feas_lo > A_SWEEP(1));
    fprintf(['\n** NOTE: on segments %s the SMALL accelerations are infeasible and the\n' ...
             '   larger ones are feasible. Pitching down lowers w_body, so acceleration\n' ...
             '   pulls those points back inside the model range. Lowest feasible a:\n'], ...
             mat2str(kk'));
    for k = kk(:)'
        fprintf('     segment %2d : feasible from a = %.2f to %.2f%s\n', ...
                k, a_feas_lo(k), a_feas(k), ternstr(feas_contig(k),'',' (with gaps)'));
    end
end

%% ========================================================================
%  PART C / TABLE 3 - model branch use and what got clamped
%% ========================================================================
line('='); fprintf('TABLE 3 - branch and clamping, evaluated at each a_feasible\n'); line('=');
fprintf('%3s | %8s | %9s | %10s | %-16s | %s\n', ...
        'k','UH','a_feas','clamped %','branch','clamped variables');
line('-');
branch_use = cell(NS,1);  clamp_pct = nan(NS,1);
for k = 1:NS
    R = evalseg(k, max(a_feas(k), A_SWEEP(1)));
    clamp_pct(k) = 100*mean(R.anyclamp);
    branch_use{k} = strjoin(unique(R.branch), ',');
    fprintf('%3d | %8.1f | %9.2f | %9.1f%% | %-16s | %s\n', ...
            k, UH(k), a_feas(k), clamp_pct(k), branch_use{k}, clampstr(R.clampcount, NU));
end

%% ========================================================================
%  CHECK 2 - hover segment, path sweep against the linear prediction
%% ========================================================================
line('='); fprintf('CHECK 2 - segment 1: the two sweep styles must nearly agree\n'); line('=');
Rp = evalseg(1, A_REF);
fprintf('  path sweep, mean a : %.4f ft/s^2\n', mean(Rp.a_horiz) - a_bias(1));
fprintf('  demanded           : %.4f ft/s^2\n', A_REF);
d1 = mean(Rp.a_horiz) - a_bias(1) - A_REF;
fprintf('  difference         : %+.4f -> %s\n', d1, ...
        ternstr(abs(d1) < 0.6, 'PASS', 'differs - check the implementation'));
fprintf(['  theta_trim moves only %.5f rad across segment 1, so the new term is\n' ...
         '  small there and the two sweep styles cannot differ much.\n'], ...
        abs(th_trim(2)-th_trim(1)));

%% ========================================================================
%  CHECK 3 - angle of attack, path sweep against fixed-ubar sweep
%% ========================================================================
line('='); fprintf('CHECK 3 - alpha during acceleration, both sweep styles\n'); line('=');
fprintf('%3s | %8s | %12s | %13s | %13s | %s\n', ...
        'k','UH','alpha_trim','alpha path','alpha fixed','recovered');
line('-');
for k = [1 6 12 16 19]
    R = evalseg(k, A_REF);
    ubt = UH(k)*cos(th_trim(k)) - WH3*sin(th_trim(k));
    al_tr = rad2deg(atan2(wb_trim(k), ubt));
    al_pa = rad2deg(R.alpha(end));
    thf = th_trim(k) + mu_theta(k)*A_REF;
    ubf = UH(k)*cos(thf) - WH3*sin(thf);
    wbf = UH(k)*sin(thf) + WH3*cos(thf);
    al_fx = rad2deg(atan2(wbf, ubf));
    fprintf('%3d | %8.1f | %12.2f | %13.2f | %13.2f | %+8.2f deg\n', ...
            k, UH(k), al_tr, al_pa, al_fx, al_pa - al_fx);
end
fprintf(['\n  "path" is evaluated at the END of the segment, where the speed and the\n' ...
         '  trim attitude have both risen; "fixed" holds ubar at the segment start,\n' ...
         '  as the earlier sweep did. A positive difference is alpha the path\n' ...
         '  sweep recovers.\n']);

%% ========================================================================
%  PART D / TABLE 2 - side by side with the earlier answers
%% ========================================================================
a_star = nan(NS,1);  a_lin_old = nan(NS,1);
p1 = fullfile(dataDir,'min_transition_time_WH3.mat');
p2 = fullfile(dataDir,'mu_theta_nonlinear.mat');
if isfile(p1), P = load(p1); a_star = P.a_star(:); end
if isfile(p2), Q = load(p2); a_lin_old = Q.dth_valid(1:NS) ./ abs(Q.mu_theta(1:NS)); end

a_final = a_feas;
if all(isfinite(a_star)), a_final = min(a_star, a_feas); end
T_k = du_seg ./ a_final;  T_k(a_final <= 0) = Inf;

line('='); fprintf('TABLE 2 - three independent limits and the combination\n'); line('=');
fprintf('%3s | %8s | %9s | %10s | %11s | %9s | %9s | %s\n', ...
        'k','UH','a_star','a_lin_old','a_feasible','a_final','T(k)','stopped by');
line('-');
for k = 1:NS
    fprintf('%3d | %8.1f | %9.2f | %10.2f | %11.2f | %9.2f | %9.2f | %s\n', ...
            k, UH(k), a_star(k), a_lin_old(k), a_feas(k), a_final(k), T_k(k), fail_why{k});
end
fprintf('%3s | %8s | %9s | %10s | %11s | %9s | %9.2f |\n', 'sum','','','','','', sum(T_k));
fprintf(['\n  a_star     reachable-set gate (grid only, no vehicle limits)\n' ...
         '  a_lin_old  fixed-ubar sweep, linear-validity range\n' ...
         '  a_feasible this file: delivered along the real path, inside the model\n' ...
         '             ranges, and inside the rotor limits\n' ...
         '  a_final    min(a_star, a_feasible)\n']);

%% ========================================================================
%  FIGURES
%% ========================================================================
shw = [1 6 12 19];
f1 = figure('Position',[50 50 1000 700]);
for j = 1:numel(shw)
    k = shw(j);  subplot(2,2,j); hold on;
    R0 = evalseg(k, 0);  R3 = evalseg(k, A_REF);
    plot(R0.u, rad2deg(R0.alpha), 'k-',  'LineWidth', 1.8);
    plot(R3.u, rad2deg(R3.alpha), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 3);
    thf = th_trim(k) + mu_theta(k)*A_REF;
    af = rad2deg(atan2(UH(k)*sin(thf)+WH3*cos(thf), UH(k)*cos(thf)-WH3*sin(thf)));
    yline(af, 'b--', 'fixed-ubar sweep', 'LineWidth', 1.2);
    yline(10*KT/max(UH(k),1)*0, 'w-');   % keep axes stable
    xlabel('u [ft/s]'); ylabel('\alpha [deg]'); grid on;
    title(sprintf('segment %d  (%.0f \\rightarrow %.0f ft/s)', k, UH(k), UH(k+1)), 'FontSize', 9);
    if j == 1, legend({'trim (a=0)', sprintf('a = %.1f', A_REF)}, 'Location','best','FontSize',8); end
end
sgtitle('angle of attack along the segment: the path sweep keeps alpha near trim');
saveas(f1, fullfile(root,'logger','path_sweep_alpha.png'));

f2 = figure('Position',[50 50 1000 520]); hold on;
plot(1:NS, a_star,    'k-o', 'LineWidth', 1.5);
plot(1:NS, a_lin_old, 'b-s', 'LineWidth', 1.5);
plot(1:NS, a_feas,    'r-^', 'LineWidth', 2.0);
xlabel('segment k'); ylabel('a [ft/s^2]'); grid on; xticks(1:NS);
legend({'a\_star (reachable set)','a\_lin\_old (fixed-ubar sweep)', ...
        'a\_feasible (path sweep + model ranges)'}, 'Location','best');
title('three independent acceleration limits');
saveas(f2, fullfile(root,'logger','path_sweep_limits.png'));

f3 = figure('Position',[50 50 1000 440]);
bar(1:NS, clamp_pct, 'FaceColor',[0.90 0.60 0.50]);
xlabel('segment k'); ylabel('points with a clamped variable [%]');
title('aero-model range violations along each segment (at a\_feasible)');
grid on; xticks(1:NS); ylim([0 105]);
saveas(f3, fullfile(root,'logger','path_sweep_clamp.png'));

%% ========================================================================
%  SAVE
%% ========================================================================
meta = struct('WH_IDX',WH_IDX,'WH3',WH3,'UH',UH,'a_sweep',A_SWEEP,'nu',NU, ...
              'ach_tol',ACH_TOL,'rot_act_max_rpm',ROT_ACT_MAX,'a_ref',A_REF, ...
              'vertical_condition','body axis: Z + m*g*cos(theta) = 0 (brief A4)', ...
              'a_horizontal','a_body*cos(theta), since wdot_b = 0 under that condition', ...
              'note','no flight data, schedule or previous total was read');
save(fullfile(dataDir,'mu_theta_along_path.mat'), ...
     'UH','du_seg','th_trim','mu_theta','wb_trim','a_bias','ach1','rat1','wst1', ...
     'vio1','a_feas','fail_why','clamp_pct','branch_use','a_star','a_lin_old', ...
     'a_final','T_k','feas_tab','c1_tab','c2_tab','c3_tab','a_feas_lo','feas_contig','meta');
fprintf('\nsaved data/mu_theta_along_path.mat\n');
fprintf('saved logger/path_sweep_{alpha,limits,clamp}.png\n');

%% ------------------------------------------------------------------------
function R = sweep_segment(k, a, UH, th_trim, mu_theta, XU0, WH3, NU, ...
                           m, G, rho, a_snd, blend, Units, CM, KT, RPM, S_LO, S_HI)
u = linspace(UH(k), UH(k+1), NU)';
f = (u - UH(k)) / (UH(k+1) - UH(k));
tht0 = (1-f)*th_trim(k)  + f*th_trim(k+1);
mut  = (1-f)*mu_theta(k) + f*mu_theta(k+1);
del0 = (1-f)*XU0(13:16,k)' + f*XU0(13:16,k+1)';
omp0 = (1-f)*XU0(17:25,k)' + f*XU0(17:25,k+1)';

R.u = u;
R.a_body = nan(NU,1); R.a_horiz = nan(NU,1); R.a_horiz_hd = nan(NU,1);
R.alpha = nan(NU,1);  R.scale = nan(NU,1);  R.rot_rpm_max = nan(NU,1);
R.anyclamp = false(NU,1);  R.branch = cell(NU,1);
R.clampcount = zeros(1,5);

for i = 1:NU
    th = tht0(i) + mut(i)*a;
    ub = u(i)*cos(th) - WH3*sin(th);
    wb = u(i)*sin(th) + WH3*cos(th);
    del = del0(i,:)';  om0 = omp0(i,:)';

    g = @(s) vert_resid(s, ub, wb, th, del, om0, rho, a_snd, blend, Units, CM, m, G);
    s = bisect(g, S_LO, S_HI);
    R.alpha(i) = atan2(wb, ub);
    [cl, br] = limit_audit(ub, wb, del, om0, KT, RPM);
    R.branch{i} = br;
    if isnan(s)
        R.anyclamp(i) = true;  R.clampcount = R.clampcount + cl;  continue;
    end
    om = om0;  om(1:8) = om0(1:8)*s;
    F = call_aero([ub;0;wb], del, om, rho, a_snd, blend, Units, CM);

    ab = F(1)/m - G*sin(th);
    R.a_body(i)     = ab;
    R.a_horiz(i)    = ab*cos(th);      % wdot_b = 0 under the A4 condition
    R.a_horiz_hd(i) = ab/cos(th);      % heading-frame alternative, for comparison
    R.scale(i)      = s;
    R.rot_rpm_max(i)= max(om(1:8))*RPM;

    [cl, br] = limit_audit(ub, wb, del, om, KT, RPM);
    R.branch{i} = br;
    R.clampcount = R.clampcount + cl;
    R.anyclamp(i) = any(cl);
end
end

function r = vert_resid(s, ub, wb, th, del, om0, rho, a_snd, blend, Units, CM, m, G)
om = om0; om(1:8) = om0(1:8)*s;
F = call_aero([ub;0;wb], del, om, rho, a_snd, blend, Units, CM);
r = F(3)/m + G*cos(th);
end

function [F, V] = call_aero(vb, del, omp, rho, a_snd, blend, Units, CM)
srf = [del(1)+del(2), del(1)-del(2), del(3), del(3), del(4)];
[X,Y,Z,~,~,~,V] = LpC_aero_p_v2([vb(:); 0; 0; 0]', omp(:)', srf, rho, a_snd, blend, Units, CM);
F = [X; Y; Z];
end

function [cl, br] = limit_audit(ub_fps, wb_fps, del, om, KT, RPM)
% Mirrors the check_fac_limits bounds in LpC_interp_p_v2 so a clamped variable
% can be named. cl = [u, w, surface, lift rotor, pusher].
uk = ub_fps/KT;  wk = wb_fps/KT;
Nl = om(1:8)*RPM;  Np = om(9)*RPM;
srf = rad2deg([del(1)+del(2), del(1)-del(2), del(3), del(4)]);
if Np >= 750
    br = 'transition';   ulim = [-5 95];
elseif uk <= 45*1.02
    br = 'hover/blend';  ulim = [-5 45];
else
    br = 'trans-cruise'; ulim = [-5 95];
end
cl = zeros(1,5);
cl(1) = uk < ulim(1) || uk > ulim(2);
cl(2) = abs(wk) > 10;
cl(3) = any(abs(srf) > 30);
cl(4) = any(Nl > 1550) || any(Nl > 1 & Nl < 550);
cl(5) = Np > 1750 || (Np > 1 && Np < 750 && uk > 45);
end

function s = bisect(g, lo, hi)
a = g(lo); b = g(hi);
if isnan(a) || isnan(b) || a*b > 0, s = NaN; return; end
for i = 1:50
    mm = 0.5*(lo+hi);  v = g(mm);
    if isnan(v), s = NaN; return; end
    if a*v <= 0, hi = mm; b = v; else, lo = mm; a = v; end
end
s = 0.5*(lo+hi);
end

function s = clampstr(cnt, n)
nm = {'u','w','surface','lift rotor','pusher'};
on = {};
for i = 1:5
    if cnt(i) > 0, on{end+1} = sprintf('%s(%d/%d)', nm{i}, cnt(i), n); end %#ok<AGROW>
end
if isempty(on), s = 'none'; else, s = strjoin(on, ' '); end
end

function s = ternstr(c,a,b), if c, s = a; else, s = b; end, end
