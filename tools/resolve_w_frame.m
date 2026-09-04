% resolve_w_frame - which frame is the BRT's w axis in, and where should the
% section be cut?
%
% The trim table's (ubar, vbar, wbar) is heading-frame. The body-frame value
% follows from
%     u_body = ubar*cos(theta) - wbar*sin(theta)
%     w_body = ubar*sin(theta) + wbar*cos(theta)
%
% What is NOT known is which of the two the reachable sets were built in. If
% they are heading-frame the aircraft sits at w deviation 0 while accelerating
% and the section is cut at 0. If they are body-frame but centred on the trim
% table's heading-frame number, the aircraft sits at a standing offset and the
% section must be cut there.
%
% PART A  states the two hypotheses and computes H2's offset from the trim
%         table alone - it is arithmetic, not a fitted number
% PART B  tries to tell the two apart from the SHAPE of the stored sets
% PART C  runs the transition calculation once per hypothesis
% PART D  separates out what neither hypothesis can excuse: accelerations that
%         put the pitch outside the range the linearisation was taken in
%
% Nothing here reads a flight log, a schedule, or any previously measured
% deviation. The only inputs are the trim table and the BRT files.
%
% *** Distances scale with the grid extent assumed by
% *** FilterConfig.channelSpec('lon'), which the BRT files do not record.
%
% Reads only. Output: data/w_frame_resolution.mat, logger/w_frame_*.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

K       = 20;
WH_IDX  = 3;
IQ0     = 31;              % q = 0
ITH0    = 16;              % theta = 0
DU_SCAN = 0.05;
A_SCAN  = 0:0.05:20;
ETAS    = [1.0 0.9 0.8 0.7];
LIN_LIM = 0.20;            % rad; edge of the range the linearisation covers
G       = 32.174;
CORR_HI = 0.7;             % above this -> supports H2
CORR_LO = 0.3;             % below this -> supports H1

line = @(c) fprintf('%s\n', repmat(c,1,104));
line('='); fprintf('GRID EXTENT IS AN ASSUMPTION (FilterConfig), NOT VERIFIED FROM THE BRT FILES\n'); line('=');

spec = FilterConfig.channelSpec('lon');
gu  = linspace(spec.grid_min(1), spec.grid_max(1), spec.grid_num(1));
gw  = linspace(spec.grid_min(2), spec.grid_max(2), spec.grid_num(2));
gth = linspace(spec.grid_min(4), spec.grid_max(4), spec.grid_num(4));

T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:K);
WH3 = T.WH(WH_IDX);
th_trim = squeeze(T.XU0_interp(11, 1:K, WH_IDX))';
du_seg = diff(UH);  NS = K-1;

F3 = cell(1,K);  UW = cell(1,K);
for k = 1:K
    S = load(fullfile(dataDir, sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    D = squeeze(S.data(:, :, IQ0, :));            % 21 x 41 x 31  (u, w, theta)
    F3{k} = griddedInterpolant({gu, gw, gth}, D, 'linear', 'none');
    UW{k} = D(:, :, ITH0);                        % 21 x 41  (u, w) at theta = 0
end

mu_theta = zeros(K,1);
for k = 1:K
    W = T.W_lon_interp(:,:,k,WH_IDX);  Bl = T.B_lon_interp(:,:,k,WH_IDX);
    M = (W \ Bl') / (Bl * (W \ Bl'));
    mu_theta(k) = M(12,1);
end

%% ========================================================================
%  PART A - the two hypotheses
%% ========================================================================
line('='); fprintf('PART A - hypotheses\n'); line('=');
fprintf(['H1: the BRT (u,w) axes are heading-frame. Accelerating trim holds wbar,\n' ...
         '    so the w deviation is 0 and the section is cut at 0.\n\n' ...
         'H2: the BRT (u,w) axes are body-frame but the sets are centred on the\n' ...
         '    trim table''s heading-frame wbar. The aircraft then stands off by\n' ...
         '      w_slice(k) = UH(k)*sin(theta) + WH3*(cos(theta) - 1)\n' ...
         '    which is arithmetic from the trim table, not a tuned quantity.\n\n']);

w_slice_H2 = UH .* sin(th_trim) + WH3 * (cos(th_trim) - 1);

%% ========================================================================
%  PART B - can the stored shapes tell the two apart?
%% ========================================================================
cen_w = nan(K,1); cen_u = nan(K,1);          % centroid of V<=0
inc_w = nan(K,1); inc_u = nan(K,1); inc_r = nan(K,1);   % largest inscribed circle
dmin_w = nan(K,1);                            % deepest point (argmin V)
ax_ang = nan(K,1);                            % major-axis angle from the u axis
r_q = nan(K,1);                               % B3 control

[UU, WW] = ndgrid(gu, gw);
for k = 1:K
    Z = UW{k};  in = Z <= 0;
    if ~any(in(:)), continue; end

    cen_u(k) = mean(UU(in));  cen_w(k) = mean(WW(in));

    % largest inscribed circle: both axes are ft/s, so plain Euclidean
    % distance is physical here - no normalisation choice to defend.
    ui = UU(in);  wi = WW(in);
    uo = UU(~in); wo = WW(~in);
    d  = min(sqrt((ui - uo').^2 + (wi - wo').^2), [], 2);
    dedge = min([ui - gu(1), gu(end) - ui, wi - gw(1), gw(end) - wi], [], 2);
    d = min(d, dedge);
    [inc_r(k), im] = max(d);
    inc_u(k) = ui(im);  inc_w(k) = wi(im);

    [~, iz] = min(Z(:));  dmin_w(k) = WW(iz);

    % principal axis of the V<=0 region
    P = [ui - cen_u(k), wi - cen_w(k)];
    C = (P' * P) / size(P,1);
    [Vc, Dc] = eig(C);
    [~, imx] = max(diag(Dc));
    ax_ang(k) = atan2d(Vc(2,imx), Vc(1,imx));
    if ax_ang(k) < 0, ax_ang(k) = ax_ang(k) + 180; end

    % B3 control: q half-width along the q axis (u=w=theta=0)
    S = load(fullfile(dataDir, sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    gq = linspace(spec.grid_min(3), spec.grid_max(3), spec.grid_num(3));
    lq = squeeze(S.data(11, 21, :, ITH0));
    Fq = griddedInterpolant(gq, lq, 'linear', 'none');
    r_q(k) = cross0(@(x) Fq(x), 0.001, gq(end));
end
ax_tilt = ax_ang - 90;                       % w is the long axis; tilt off it

line('='); fprintf('TABLE 1 - trim geometry vs measured section centres\n'); line('=');
fprintf('%3s | %7s | %11s | %11s | %10s %10s %10s | %9s | %8s\n', ...
        'k','UH','th_trim[d]','w_slice_H2','centroid_w','incirc_w','deepest_w', ...
        'axis[deg]','r_q');
line('-');
for k = 1:K
    fprintf('%3d | %7.1f | %11.3f | %11.3f | %10.3f %10.3f %10.3f | %9.2f | %8.3f\n', ...
            k, UH(k), rad2deg(th_trim(k)), w_slice_H2(k), cen_w(k), inc_w(k), dmin_w(k), ...
            ax_ang(k), r_q(k));
end

cc_cen  = corr_(w_slice_H2, cen_w);
cc_inc  = corr_(w_slice_H2, inc_w);
cc_deep = corr_(w_slice_H2, dmin_w);
cc_ax   = corr_(rad2deg(th_trim), ax_tilt);

fprintf('\n--- PART B correlations ---\n');
fprintf('  w_slice_H2  vs  centroid w      : r = %+.3f\n', cc_cen);
fprintf('  w_slice_H2  vs  inscribed-circle w : r = %+.3f\n', cc_inc);
fprintf('  w_slice_H2  vs  deepest-point w : r = %+.3f\n', cc_deep);
fprintf('  theta_trim  vs  major-axis tilt : r = %+.3f\n', cc_ax);
fprintf('\n  magnitudes: max |w_slice_H2| = %.2f ft/s, max |centroid w| = %.2f ft/s,\n', ...
        max(abs(w_slice_H2)), max(abs(cen_w)));
fprintf('              w grid step = %.2f ft/s\n', gw(2)-gw(1));
fprintf(['\n  B3 control: q is the same quantity in both frames, so r_q cannot\n' ...
         '  discriminate. Range %.3f .. %.3f rad/s, spread %.1f %% - %s\n'], ...
        min(r_q), max(r_q), 100*(max(r_q)/min(r_q)-1), ...
        ternary(abs(corr_(rad2deg(th_trim), r_q)) < 0.5, ...
                'no strong tie to theta_trim, as expected', ...
                'unexpectedly tied to theta_trim - investigate'));

%% ========================================================================
%  PART C - the transition calculation under each hypothesis
%% ========================================================================
HYP = {'H1 (w slice = 0)', 'H2 (w slice = computed)'};
a_star = zeros(NS,2);  pitch_at = zeros(NS,2);
Tseg = nan(NS,numel(ETAS),2);  Ttot = nan(numel(ETAS),2);

for h = 1:2
    for k = 1:NS
        ua = (UH(k):DU_SCAN:UH(k+1))';
        if ua(end) < UH(k+1)-1e-9, ua(end+1) = UH(k+1); end %#ok<SAGROW>
        tht = interp1(UH([k k+1]), th_trim([k k+1]),  ua);
        mut = interp1(UH([k k+1]), mu_theta([k k+1]), ua);
        if h == 1
            wsl = zeros(size(ua));
        else
            % same formula as PART A, evaluated continuously along the ramp;
            % it reduces to w_slice_H2(k) at the anchors
            wsl = ua .* sin(tht) + WH3 * (cos(tht) - 1);
        end
        dk = ua - UH(k);  dk1 = ua - UH(k+1);
        shk = tht - th_trim(k);  shk1 = tht - th_trim(k+1);

        best = 0;  bp = 0;
        for a = A_SCAN
            vk  = F3{k  }(dk,  wsl, shk  + mut*a);
            vk1 = F3{k+1}(dk1, wsl, shk1 + mut*a);
            ik  = ~isnan(vk)  & vk  <= 0;
            ik1 = ~isnan(vk1) & vk1 <= 0;
            if ~ik(1) || ~ik1(end), continue; end
            b  = find(~ik, 1, 'first');   isup = ternnum(isempty(b),  numel(ua), b-1);
            b2 = find(~ik1, 1, 'last');   iinf = ternnum(isempty(b2), 1,        b2+1);
            if ua(isup) - ua(iinf) >= 0
                best = a;  bp = max(abs(shk + mut*a));
            end
        end
        a_star(k,h) = best;  pitch_at(k,h) = bp;
    end
    for e = 1:numel(ETAS)
        Tseg(:,e,h) = du_seg ./ max(ETAS(e)*a_star(:,h), eps);
        Ttot(e,h) = sum(Tseg(:,e,h));
    end
end

line('='); fprintf('TABLE 2 - a_star and the pitch it demands, per hypothesis\n'); line('=');
fprintf('%3s | %7s | %11s %11s | %14s %14s | %s\n', ...
        'k','UH','a*(H1)','a*(H2)','pitch H1[deg]','pitch H2[deg]','beyond linear range');
line('-');
for k = 1:NS
    fl = '';
    if pitch_at(k,1) > LIN_LIM, fl = [fl 'H1 ']; end %#ok<AGROW>
    if pitch_at(k,2) > LIN_LIM, fl = [fl 'H2 ']; end %#ok<AGROW>
    fprintf('%3d | %7.1f | %11.2f %11.2f | %14.2f %14.2f | %s\n', ...
            k, UH(k), a_star(k,1), a_star(k,2), ...
            rad2deg(pitch_at(k,1)), rad2deg(pitch_at(k,2)), fl);
end

line('='); fprintf('TABLE 3 - total transition time per hypothesis\n'); line('=');
fprintf('%-26s | %6s | %12s | %s\n', 'hypothesis','eta','T_total [s]','segments beyond linear range');
line('-');
for h = 1:2
    nbad = sum(pitch_at(:,h) > LIN_LIM);
    for e = 1:numel(ETAS)
        fprintf('%-26s | %6.1f | %12.2f | %d of %d\n', ...
                ternstr(e==1, HYP{h}, ''), ETAS(e), Ttot(e,h), nbad, NS);
    end
end

%% ========================================================================
%  PART D - the acceleration at which the linearisation stops applying
%% ========================================================================
a_lin = nan(NS,1);
for k = 1:NS
    ua = (UH(k):DU_SCAN:UH(k+1))';
    tht = interp1(UH([k k+1]), th_trim([k k+1]),  ua);
    mut = interp1(UH([k k+1]), mu_theta([k k+1]), ua);
    shk = tht - th_trim(k);
    a_lin(k) = A_SCAN(end);
    for a = A_SCAN
        if max(abs(shk + mut*a)) >= LIN_LIM, a_lin(k) = a; break; end
    end
end

line('='); fprintf('PART D - restricting to where the linearisation is valid\n'); line('=');
fprintf('D1. segments whose a_star demands more than %.2f rad of pitch:\n', LIN_LIM);
for h = 1:2
    bad = find(pitch_at(:,h) > LIN_LIM);
    fprintf('   %-26s : %s\n', HYP{h}, ternstr(isempty(bad),'none',mat2str(bad')));
end

fprintf('\nTABLE 4 - a_safe = min(a_star, a_linear_limit)\n');
fprintf('%3s | %10s %10s | %10s | %10s %10s | %8s %8s\n', ...
        'k','a*(H1)','a*(H2)','a_lin','a_safe(H1)','a_safe(H2)','T(H1)','T(H2)');
line('-');
a_safe = min(a_star, a_lin);
Tsafe = du_seg ./ max(a_safe, eps);
for k = 1:NS
    fprintf('%3d | %10.2f %10.2f | %10.2f | %10.2f %10.2f | %8.2f %8.2f\n', ...
            k, a_star(k,1), a_star(k,2), a_lin(k), a_safe(k,1), a_safe(k,2), ...
            Tsafe(k,1), Tsafe(k,2));
end
fprintf('%3s | %10s %10s | %10s | %10s %10s | %8.2f %8.2f\n', 'sum','','','','','', ...
        sum(Tsafe(:,1)), sum(Tsafe(:,2)));
fprintf(['\nThose two totals are the minimum transition time inside the range the\n' ...
         'linearisation was taken in, at eta = 1. They still carry no actuator limit.\n']);

%% ========================================================================
%  CHECKS
%% ========================================================================
line('='); fprintf('CHECK 1 - H2 offset must vanish at rest\n'); line('=');
fprintf('  UH1: UH = %.1f ft/s, theta_trim = %.3f deg  ->  w_slice_H2 = %+.5f ft/s\n', ...
        UH(1), rad2deg(th_trim(1)), w_slice_H2(1));
fprintf('  %s\n', ternstr(abs(w_slice_H2(1)) < 0.01, ...
        'PASS - the two frames coincide at rest, as they must.', ...
        'FAIL - the frames should coincide at rest; check the formula.'));

line('='); fprintf('CHECK 2 - does the shape evidence pick a hypothesis?\n'); line('=');
cc_best = max(abs([cc_cen cc_inc cc_deep]));
fprintf('  strongest centre correlation with w_slice_H2 : |r| = %.3f\n', cc_best);
fprintf('  theta_trim vs major-axis tilt                : |r| = %.3f\n', abs(cc_ax));
if cc_best > CORR_HI
    verdict = 'supports H2 (body frame)';
elseif cc_best < CORR_LO
    verdict = 'supports H1 (heading frame)';
else
    verdict = 'INCONCLUSIVE';
end
fprintf('  verdict from the centre test: %s\n', verdict);

line('='); fprintf('CHECK 3 - conclusion\n'); line('=');
if strcmp(verdict, 'INCONCLUSIVE')
    fprintf(['  The shape evidence does not separate the two hypotheses. No choice is\n' ...
             '  forced here. Both PART C results are kept and both are reported.\n' ...
             '  Settling it needs the BRT generating script, which is not in this\n' ...
             '  repository (tools/brt_inspect.m item 7 found consumers only).\n']);
else
    fprintf('  The centre test %s. Note this is indirect evidence about the\n', verdict);
    fprintf(['  generator''s frame, not a reading of it. The generating script would\n' ...
             '  settle it directly and is not in this repository.\n']);
end
fprintf('\n  Both hypotheses are saved. Neither PART C column is preferred here.\n');

%% ========================================================================
%  FIGURES
%% ========================================================================
f1 = figure('Position',[30 30 1500 950]);
for k = 1:K
    subplot(4,5,k);
    contourf(gu, gw, UW{k}', 18, 'LineColor','none'); hold on;
    contour(gu, gw, UW{k}', [0 0], 'k-', 'LineWidth', 1.6);
    plot(0, 0, 'wp', 'MarkerSize', 9, 'MarkerFaceColor','w');
    yline(w_slice_H2(k), 'r--', 'LineWidth', 1.4);
    plot(cen_u(k), cen_w(k), 'g+', 'MarkerSize', 9, 'LineWidth', 1.5);
    title(sprintf('UH%d  u=%.0f', k, UH(k)), 'FontSize', 8);
    xlim([-16 16]); ylim([-33 33]); set(gca,'FontSize',7);
    if k > 15, xlabel('du [ft/s]','FontSize',7); end
    if mod(k,5)==1, ylabel('dw [ft/s]','FontSize',7); end
end
sgtitle('(u,w) sections at q=0, theta=0.  white = origin,  red = H2 slice,  green = centroid');
saveas(f1, fullfile(root,'logger','w_frame_sections.png'));

f2 = figure('Position',[60 60 900 520]); hold on;
plot(rad2deg(th_trim), cen_w,  'bo', 'MarkerSize', 7, 'LineWidth', 1.3);
plot(rad2deg(th_trim), inc_w,  'rs', 'MarkerSize', 7, 'LineWidth', 1.3);
plot(rad2deg(th_trim), dmin_w, 'g^', 'MarkerSize', 7, 'LineWidth', 1.3);
plot(rad2deg(th_trim), w_slice_H2, 'k-', 'LineWidth', 2);
xlabel('\theta_{trim} [deg]'); ylabel('w [ft/s]');
title(sprintf(['section centre vs trim pitch  (r: centroid %+.2f, incircle %+.2f, ' ...
               'deepest %+.2f)'], cc_cen, cc_inc, cc_deep));
legend({'centroid','inscribed-circle centre','deepest point','H2 predicted slice'}, ...
       'Location','northwest');
grid on;
saveas(f2, fullfile(root,'logger','w_frame_correlation.png'));

f3 = figure('Position',[60 60 1000 520]); hold on;
plot(1:NS, a_star(:,1), 'b-o', 'LineWidth', 1.6);
plot(1:NS, a_star(:,2), 'r-s', 'LineWidth', 1.6);
plot(1:NS, a_lin, 'k--', 'LineWidth', 1.4);
xlabel('segment k'); ylabel('a [ft/s^2]');
title('largest acceleration the sets allow, per hypothesis (dashed = linear-range limit)');
legend({HYP{1}, HYP{2}, 'linear-range limit'}, 'Location','best');
grid on; xticks(1:NS);
saveas(f3, fullfile(root,'logger','w_frame_astar.png'));

%% ========================================================================
%  SAVE
%% ========================================================================
meta = struct('WH_IDX', WH_IDX, 'WH3', WH3, 'UH', UH, 'etas', ETAS, ...
              'lin_lim', LIN_LIM, 'a_scan', A_SCAN, 'du_scan', DU_SCAN, ...
              'hypotheses', {HYP}, 'verdict', verdict, ...
              'grid_extent_assumed', true, 'actuator_limits_included', false, ...
              'note', ['no flight data, schedule or previously measured deviation ' ...
                       'was used; w_slice_H2 is arithmetic from the trim table']);
out = fullfile(dataDir,'w_frame_resolution.mat');
save(out, 'UH','th_trim','w_slice_H2','cen_u','cen_w','inc_u','inc_w','inc_r', ...
     'dmin_w','ax_ang','ax_tilt','r_q','cc_cen','cc_inc','cc_deep','cc_ax', ...
     'a_star','pitch_at','a_lin','a_safe','Tseg','Ttot','Tsafe','mu_theta','meta');

fprintf('\nsaved %s\n', out);
fprintf('saved logger/w_frame_{sections,correlation,astar}.png\n');
line('='); fprintf('No actuator limits here. Distances scale with the assumed grid extent.\n'); line('=');

%% ------------------------------------------------------------------------
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
function s = ternstr(c,a,b), if c, s = a; else, s = b; end, end
function v = ternnum(c,a,b), if c, v = a; else, v = b; end, end

function r = corr_(x, y)
x = x(:) - mean(x(:));  y = y(:) - mean(y(:));
d = sqrt(sum(x.^2) * sum(y.^2));
if d == 0, r = NaN; else, r = sum(x.*y)/d; end
end

function r = cross0(fh, step, lim)
if fh(0) > 0, r = 0; return; end
x = 0;
while x + step <= lim
    if fh(x+step) > 0
        lo = x; hi = x+step;
        for i = 1:40, m = (lo+hi)/2; if fh(m) <= 0, lo = m; else, hi = m; end, end
        r = lo; return;
    end
    x = x + step;
end
r = lim;
end
