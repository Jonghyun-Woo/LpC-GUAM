% extract_brt_section - (forward speed, pitch attitude) sections of the 20
% longitudinal reachable sets, and the handover gate read straight off the
% grid.
%
% PART A  slice each set at w = 0, q = 0 and measure its extent by
%         interpolation (not by node snapping)
% PART B  THE RESULT: place two neighbouring sets at their true trim speeds
%         and measure, on the grid itself, how wide the region is that
%         belongs to both. No ellipse, no fit, no linearisation.
% PART C  inscribed ellipse, for reference only, plus how much it costs
%
% PART B is the conclusion. PART C never decides a gate.
%
% Assumptions carried over from tools/brt_inspect.m:
%   array order [u, w, q, theta], size [21 41 61 31]; V <= 0 is inside;
%   coordinates are deviations from the trim point of that file.
%
%   *** The physical grid extent is taken from FilterConfig.channelSpec('lon')
%   *** and is NOT verifiable from the data files (they store only the value
%   *** array). Every distance below scales with that assumption.
%
% Reads only. Output: data/brt_section_WH3.mat, logger/brt_section_*.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

K        = 20;
WH_IDX   = 3;
IW0      = 21;          % w = 0   (41 nodes over [-33, 33])
IQ0      = 31;          % q = 0   (61 nodes over [-1.5, 1.5])
DU_SCAN  = 0.05;        % [ft/s]  sweep step for the gate
TH_SCAN  = -0.20:0.01:0.20;   % 41 pitch values - the requested table
% The sets reach +-0.23..0.35 rad in pitch, so a +-0.20 sweep cannot find
% where the overlap actually ends: the band saturates at the sweep limit and
% reports the sweep, not the data. The wide sweep below runs past the widest
% section so theta_lo/theta_hi are set by the sets themselves. TH_SCAN is
% still evaluated and stored unchanged, as the 19x41 table.
TH_WIDE  = -0.40:0.01:0.40;   % 81 pitch values - used for the reported band
NELL     = 200;         % ellipse boundary samples
RES_MULT = 3;           % half-width below this many grid steps -> suspect
AREA_MIN = 0.6;         % ellipse/section cell ratio below this -> over-conservative
JUMP_REL = 0.30;        % neighbour-to-neighbour change above this -> warn

line = @(c) fprintf('%s\n', repmat(c,1,100));
line('='); fprintf('GRID EXTENT IS AN ASSUMPTION (FilterConfig.channelSpec), NOT VERIFIED FROM DATA\n'); line('=');

%% ========================================================================
%  SETUP - grid vectors, trim speeds, the 20 sections
%% ========================================================================
spec = FilterConfig.channelSpec('lon');
gu  = linspace(spec.grid_min(1), spec.grid_max(1), spec.grid_num(1));   % 21
gth = linspace(spec.grid_min(4), spec.grid_max(4), spec.grid_num(4));   % 31
hu  = gu(2)-gu(1);      % 1.60 ft/s
hth = gth(2)-gth(1);    % 0.05 rad
IU0  = find(abs(gu)  < 1e-12);
ITH0 = find(abs(gth) < 1e-12);

T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:K);

Z = cell(1,K);  F = cell(1,K);
for k = 1:K
    S = load(fullfile(dataDir, sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    Z{k} = squeeze(S.data(:, IW0, IQ0, :));            % 21 x 31, (u, theta)
    F{k} = griddedInterpolant({gu, gth}, Z{k}, 'linear', 'none');
end
fprintf('loaded %d sections, each %s at w=0, q=0\n', K, mat2str(size(Z{1})));
fprintf('grid step: u %.2f ft/s, theta %.3f rad (%.2f deg)\n\n', hu, hth, rad2deg(hth));

%% ========================================================================
%  PART A - extent of each section, connectivity, resolution
%% ========================================================================
r_u_pos = nan(K,1); r_u_neg = nan(K,1);
r_th_pos = nan(K,1); r_th_neg = nan(K,1);
n_cells = zeros(K,1);  n_blobs = zeros(K,1);  centre_val = nan(K,1);

for k = 1:K
    centre_val(k) = F{k}(0,0);

    % A3 - zero crossing along each axis through the centre, by bisection
    r_u_pos(k)  = cross_out(@(x) F{k}(x, 0), 0, +1, gu(end),  DU_SCAN/4);
    r_u_neg(k)  = cross_out(@(x) F{k}(x, 0), 0, -1, gu(1),    DU_SCAN/4);
    r_th_pos(k) = cross_out(@(x) F{k}(0, x), 0, +1, gth(end), hth/100);
    r_th_neg(k) = cross_out(@(x) F{k}(0, x), 0, -1, gth(1),   hth/100);

    % A4 - is V<=0 one blob containing the centre?
    [lab, nb] = blobs(Z{k} <= 0);
    n_blobs(k) = nb;
    if nb == 0
        n_cells(k) = 0;
    else
        cl = lab(IU0, ITH0);
        if cl == 0
            fprintf('  ** UH%d: the grid centre is NOT inside V<=0. Section unusable.\n', k);
            n_cells(k) = 0;
        else
            n_cells(k) = sum(lab(:) == cl);
            if nb > 1
                fprintf(['  ** WARNING UH%d: V<=0 splits into %d pieces. Only the piece\n' ...
                         '     containing the trim point is counted (%d of %d cells).\n'], ...
                         k, nb, n_cells(k), sum(Z{k}(:) <= 0));
            end
        end
    end
end

r_u  = min(r_u_pos,  abs(r_u_neg));      % tighter side
r_th = min(r_th_pos, abs(r_th_neg));
res_u  = r_u  / hu;
res_th = r_th / hth;
resolution_flag = [res_u < RES_MULT, res_th < RES_MULT];

fprintf('\n');
line('='); fprintf('TABLE 1 - section extent per trim point (interpolated, not node-snapped)\n'); line('=');
fprintf('%3s | %8s | %8s %8s | %9s %9s | %7s %7s | %s\n', ...
        'k','UH[ft/s]','r_u_pos','r_u_neg','r_th_pos','r_th_neg','u/step','th/step','flag');
line('-');
for k = 1:K
    fl = '';
    if resolution_flag(k,1), fl = [fl 'u:RES ']; end %#ok<AGROW>
    if resolution_flag(k,2), fl = [fl 'theta:RES ']; end %#ok<AGROW>
    if n_blobs(k) > 1, fl = [fl 'split ']; end %#ok<AGROW>
    fprintf('%3d | %8.1f | %8.3f %8.3f | %9.4f %9.4f | %7.2f %7.2f | %s\n', ...
            k, UH(k), r_u_pos(k), r_u_neg(k), r_th_pos(k), r_th_neg(k), ...
            res_u(k), res_th(k), fl);
end
fprintf(['\n"u/step" and "th/step" are the tighter half-width divided by the grid step.\n' ...
         'Below %d the set is only a few nodes wide and the interpolated boundary is\n' ...
         'as much a property of the grid as of the vehicle - the set may be LARGER\n' ...
         'than it looks here.\n'], RES_MULT);

%% ========================================================================
%  PART B - THE GATE, read directly off the grid (no ellipse)
%% ========================================================================
NS = K-1;  NT = numel(TH_SCAN);  NW = numel(TH_WIDE);
overlap_table      = zeros(NS, NT);   % over TH_SCAN, as requested
gap_table          = zeros(NS, NT);
overlap_table_wide = zeros(NS, NW);   % over TH_WIDE, sets the reported band
overlap_width_theta0 = zeros(NS,1);
gap_width_theta0     = zeros(NS,1);
theta_lo = nan(NS,1);  theta_hi = nan(NS,1);
theta_contig = true(NS,1);  theta_saturated = false(NS,1);
gap_u_range  = nan(NS,2);
du_seg = diff(UH);

for k = 1:NS
    % Midpoint sampling: N cells of equal width across the segment, so the
    % measured lengths add up to exactly du (a trailing partial step would
    % otherwise let "overlap" exceed the segment it is measured over).
    du = du_seg(k);
    N  = max(round(du/DU_SCAN), 1);
    ed = linspace(UH(k), UH(k+1), N+1);
    ua = 0.5*(ed(1:end-1) + ed(2:end));
    wcell = du / N;

    for j = 1:NT+NW
        if j <= NT, th = TH_SCAN(j); else, th = TH_WIDE(j-NT); end
        vk  = F{k  }(ua(:) - UH(k),   th*ones(N,1));
        vk1 = F{k+1}(ua(:) - UH(k+1), th*ones(N,1));
        ink  = ~isnan(vk)  & vk  <= 0;
        ink1 = ~isnan(vk1) & vk1 <= 0;
        both    = ink & ink1;
        neither = ~ink & ~ink1;

        if j <= NT
            overlap_table(k,j) = wcell * sum(both);
            gap_table(k,j)     = wcell * sum(neither);
            if abs(th) < 1e-12
                overlap_width_theta0(k) = overlap_table(k,j);
                gap_width_theta0(k)     = gap_table(k,j);
                if any(neither)
                    gap_u_range(k,:) = [ua(find(neither,1,'first')) - 0.5*wcell, ...
                                        ua(find(neither,1,'last'))  + 0.5*wcell];
                end
            end
        else
            overlap_table_wide(k,j-NT) = wcell * sum(both);
        end
    end

    ok = overlap_table_wide(k,:) > 0;
    if any(ok)
        i1 = find(ok,1,'first');  i2 = find(ok,1,'last');
        theta_lo(k) = TH_WIDE(i1);
        theta_hi(k) = TH_WIDE(i2);
        theta_contig(k)   = all(ok(i1:i2));
        theta_saturated(k) = (i1 == 1) || (i2 == NW);
    end
end
if any(theta_saturated)
    fprintf(['\n  ** WARNING: the pitch band still reaches the sweep limit (%+.2f rad)\n' ...
             '     on segments %s. Those bands are lower bounds set by the sweep,\n' ...
             '     not by the data. Widen TH_WIDE to resolve them.\n'], ...
             TH_WIDE(end), mat2str(find(theta_saturated)'));
end

line('='); fprintf('TABLE 2 - handover gate per segment (grid-direct, PART B - THIS IS THE RESULT)\n'); line('=');
fprintf('%3s | %7s | %14s | %10s | %17s | %s\n', ...
        'k','du','overlap(th=0)','gap(th=0)','theta band [rad]','verdict');
line('-');
for k = 1:NS
    if isnan(theta_lo(k))
        band = sprintf('%17s','none');
    else
        band = sprintf('%+7.3f .. %+7.3f', theta_lo(k), theta_hi(k));
    end
    if gap_width_theta0(k) > 0
        vd = sprintf('GATE FAIL (gap %.2f ft/s)', gap_width_theta0(k));
    elseif overlap_width_theta0(k) > 0
        vd = 'pass';
    else
        vd = 'GATE FAIL (touching, no overlap)';
    end
    if ~theta_contig(k), vd = [vd ' [theta band not contiguous]']; end %#ok<AGROW>
    if theta_saturated(k), vd = [vd ' [band hit sweep limit]']; end %#ok<AGROW>
    fprintf('%3d | %7.3f | %14.3f | %10.3f | %17s | %s\n', ...
            k, du_seg(k), overlap_width_theta0(k), gap_width_theta0(k), band, vd);
end
fprintf('\noverlap can never exceed du; %s\n', ...
        ternary(all(overlap_width_theta0 <= du_seg + 1e-9), ...
                'checked, it does not.', '*** IT DOES - measurement bug. ***'));
fprintf(['\noverlap = length of the u interval (at that pitch) that lies inside BOTH\n' ...
         'sets at once. gap = length that lies inside NEITHER. The theta band is the\n' ...
         'range of pitch deviation over which some overlap still exists - that is the\n' ...
         'pitch budget the acceleration calculation may spend on this segment.\n']);

%% ========================================================================
%  PART C - inscribed ellipse, reference only
%% ========================================================================
s_star = nan(K,1);  rho_u = nan(K,1);  rho_theta = nan(K,1);
area_ratio = nan(K,1);  ell_maxV = nan(K,1);
tt = linspace(0, 2*pi, NELL+1);  tt(end) = [];

for k = 1:K
    if r_u(k) <= 0 || r_th(k) <= 0, continue; end
    lo = 0.05;  hi = 1.0;
    if ~ell_ok(F{k}, r_u(k), r_th(k), lo, tt), s_star(k) = NaN; continue; end
    if  ell_ok(F{k}, r_u(k), r_th(k), hi, tt)
        s_star(k) = hi;
    else
        for it = 1:30
            mid = 0.5*(lo+hi);
            if ell_ok(F{k}, r_u(k), r_th(k), mid, tt), lo = mid; else, hi = mid; end
        end
        s_star(k) = lo;
    end
    rho_u(k)     = s_star(k) * r_u(k);
    rho_theta(k) = s_star(k) * r_th(k);
    ell_maxV(k)  = max(F{k}(s_star(k)*r_u(k)*cos(tt(:)), s_star(k)*r_th(k)*sin(tt(:))));

    [UU, TT2] = ndgrid(gu, gth);
    in_ell = (UU/rho_u(k)).^2 + (TT2/rho_theta(k)).^2 <= 1;
    area_ratio(k) = sum(in_ell(:)) / max(n_cells(k), 1);
end

line('='); fprintf('TABLE 3 - grid-direct gate vs ellipse prediction (PART C, reference only)\n'); line('=');
fprintf('%3s | %8s | %8s | %14s | %14s | %9s | %9s | %s\n', ...
        'k','rho_u(k)','rho_u(k+1)','grid overlap','ellipse pred','diff','area ratio','note');
line('-');
for k = 1:NS
    pred = rho_u(k) + rho_u(k+1) - du_seg(k);
    nt = '';
    if area_ratio(k) < AREA_MIN, nt = [nt sprintf('ellipse over-conservative (%.2f) ', area_ratio(k))]; end %#ok<AGROW>
    fprintf('%3d | %8.3f | %8.3f | %14.3f | %14.3f | %+9.3f | %9.2f | %s\n', ...
            k, rho_u(k), rho_u(k+1), overlap_width_theta0(k), pred, ...
            pred - overlap_width_theta0(k), area_ratio(k), nt);
end
fprintf(['\n"area ratio" = ellipse cells / section cells. Below %.1f the ellipse is\n' ...
         'throwing away more than %.0f%% of the certified region.\n'], AREA_MIN, 100*(1-AREA_MIN));

%% ========================================================================
%  CHECKS
%% ========================================================================
line('='); fprintf('CHECK A - segments that fail the gate\n'); line('=');
failk = find(gap_width_theta0 > 0 | overlap_width_theta0 <= 0);
if isempty(failk)
    fprintf('  none. every segment has a strictly positive overlap at theta = 0.\n');
else
    for k = failk(:)'
        fprintf('  segment %2d (UH %.1f -> %.1f):  gap %.3f ft/s', ...
                k, UH(k), UH(k+1), gap_width_theta0(k));
        if ~any(isnan(gap_u_range(k,:)))
            fprintf('   over u = %.2f .. %.2f ft/s', gap_u_range(k,1), gap_u_range(k,2));
        end
        fprintf('\n');
    end
end

line('='); fprintf('CHECK B - resolution-limited trim points, and the segments they touch\n'); line('=');
resk = find(any(resolution_flag,2));
if isempty(resk)
    fprintf('  none. every half-width is at least %d grid steps.\n', RES_MULT);
else
    fprintf('%5s | %8s | %9s | %9s | %s\n','k','UH','u/step','th/step','segments it takes part in');
    line('-');
    for k = resk(:)'
        segs = intersect([k-1, k], 1:NS);
        tag = '';
        if any(ismember(segs, failk))
            tag = '  <- ALSO A GATE FAILURE: recomputing this BRT on a finer grid may fix it';
        end
        fprintf('%5d | %8.1f | %9.2f | %9.2f | %s%s\n', k, UH(k), res_u(k), res_th(k), ...
                mat2str(segs), tag);
    end
end

line('='); fprintf('CHECK C - value on the inscribed-ellipse boundary (must be <= 0)\n'); line('=');
fprintf('  max over all %d boundary points, all %d trim points : %+.3e   -> %s\n', ...
        NELL, K, max(ell_maxV), ternary(max(ell_maxV) <= 0, 'PASS', 'FAIL'));
bad = find(ell_maxV > 0);
if ~isempty(bad)
    fprintf('  violating trim points: %s\n', mat2str(bad'));
end

line('='); fprintf('CHECK D - smoothness of r_u and r_theta along UH\n'); line('=');
for nm = {'r_u','r_th'}
    v = eval(nm{1});
    d = abs(diff(v)) ./ max(abs(v(1:end-1)), eps);
    jk = find(d > JUMP_REL);
    fprintf('  %-5s : ', nm{1});
    if isempty(jk)
        fprintf('no neighbour-to-neighbour change above %.0f%%\n', 100*JUMP_REL);
    else
        fprintf('\n');
        for j = jk(:)'
            fprintf('     UH%2d -> UH%2d (%.1f -> %.1f ft/s) : %.4f -> %.4f  (%+.0f%%)\n', ...
                    j, j+1, UH(j), UH(j+1), v(j), v(j+1), 100*(v(j+1)/v(j)-1));
        end
    end
end
fprintf(['  reported as found. Nothing is smoothed. Jumps in UH 17..20 are expected\n' ...
         '  (theta_trim peaks at UH17 and the trim table is spliced at UH20/21).\n']);

%% ========================================================================
%  FIGURES
%% ========================================================================
f1 = figure('Position',[40 40 1500 950]);
for k = 1:K
    subplot(4,5,k);
    contourf(gu, gth, Z{k}', 20, 'LineColor','none'); hold on;
    contour(gu, gth, Z{k}', [0 0], 'k-', 'LineWidth', 1.8);
    if ~isnan(s_star(k))
        plot(rho_u(k)*cos(tt), rho_theta(k)*sin(tt), 'r-', 'LineWidth', 1.4);
    end
    plot(0,0,'rp','MarkerSize',8,'MarkerFaceColor','r');
    title(sprintf('UH%d  u=%.0f', k, UH(k)), 'FontSize', 8);
    xlim([-16 16]); ylim([-0.75 0.75]);
    if k > 15, xlabel('du [ft/s]','FontSize',7); end
    if mod(k,5)==1, ylabel('d\theta [rad]','FontSize',7); end
    set(gca,'FontSize',7);
end
sgtitle('BRT sections at w=0, q=0 (black = boundary, red = inscribed ellipse)');
saveas(f1, fullfile(root,'logger','brt_section_a_sections.png'));

f2 = figure('Position',[40 40 1500 950]);
for k = 1:NS
    subplot(4,5,k); hold on;
    uA = UH(k)   + gu;  uB = UH(k+1) + gu;
    contour(uA, gth, Z{k}',   [0 0], 'b-', 'LineWidth', 1.5);
    contour(uB, gth, Z{k+1}', [0 0], 'r-', 'LineWidth', 1.5);
    yline(0,'k:');
    ua = UH(k):DU_SCAN:UH(k+1);
    ink  = F{k  }(ua(:)-UH(k),   zeros(numel(ua),1)) <= 0;
    ink1 = F{k+1}(ua(:)-UH(k+1), zeros(numel(ua),1)) <= 0;
    b = ink & ink1;
    if any(b)
        plot([ua(find(b,1,'first')) ua(find(b,1,'last'))], [0 0], 'g-', 'LineWidth', 4);
    end
    g = ~ink & ~ink1;
    if any(g)
        plot([ua(find(g,1,'first')) ua(find(g,1,'last'))], [0 0], 'm-', 'LineWidth', 4);
    end
    xlim([UH(k)-6, UH(k+1)+6]); ylim([-0.4 0.4]);
    title(sprintf('seg %d: %.0f\\rightarrow%.0f', k, UH(k), UH(k+1)), 'FontSize', 8);
    set(gca,'FontSize',7);
end
sgtitle('neighbouring sets in absolute speed (blue = k, red = k+1, green = overlap, magenta = gap)');
saveas(f2, fullfile(root,'logger','brt_section_b_gate.png'));

f3 = figure('Position',[60 60 980 560]);
imagesc(1:NS, TH_WIDE, overlap_table_wide'); set(gca,'YDir','normal'); hold on;
plot(1:NS, theta_lo, 'w-o', 'LineWidth', 1.6, 'MarkerSize', 4);
plot(1:NS, theta_hi, 'w-o', 'LineWidth', 1.6, 'MarkerSize', 4);
colorbar; xlabel('segment k'); ylabel('\theta deviation [rad]');
title('overlap width [ft/s] vs pitch deviation (white = edge of the usable pitch band)');
saveas(f3, fullfile(root,'logger','brt_section_c_overlap_map.png'));

%% ========================================================================
%  SAVE
%% ========================================================================
meta = struct('WH_IDX', WH_IDX, 'grid_u', gu, 'grid_theta', gth, ...
              'theta_scan', TH_SCAN, 'theta_wide', TH_WIDE, ...
              'du_scan', DU_SCAN, 'UH', UH, 'grid_extent_assumed', true, ...
              'note', ['distances scale with FilterConfig.channelSpec(''lon'') ' ...
                       'grid bounds, which are not verifiable from the BRT files; ' ...
                       'theta_lo/theta_hi come from overlap_table_wide, not from ' ...
                       'the +-0.2 overlap_table']);
out = fullfile(dataDir,'brt_section_WH3.mat');
save(out, 'r_u_pos','r_u_neg','r_th_pos','r_th_neg','r_u','r_th', ...
     'rho_u','rho_theta','s_star','area_ratio','resolution_flag', ...
     'overlap_width_theta0','gap_width_theta0','overlap_table','gap_table', ...
     'overlap_table_wide','theta_lo','theta_hi','theta_contig', ...
     'theta_saturated','n_cells','n_blobs','du_seg','meta');

fprintf('\nsaved %s\n', out);
fprintf('saved logger/brt_section_{a_sections,b_gate,c_overlap_map}.png\n');
line('='); fprintf('REMINDER: every distance above scales with the assumed grid extent.\n'); line('=');

%% ------------------------------------------------------------------------
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end

function r = cross_out(fh, x0, dir, xlim_, step)
% March out from x0 until the value crosses zero, then bisect. Returns the
% crossing coordinate (signed). NaN if the centre is already outside; the
% grid limit if no crossing is found inside the grid.
if fh(x0) > 0, r = NaN; return; end
x = x0;
while true
    xn = x + dir*step;
    if (dir > 0 && xn > xlim_) || (dir < 0 && xn < xlim_)
        r = xlim_; return;                      % never crosses inside the grid
    end
    v = fh(xn);
    if isnan(v), r = x; return; end
    if v > 0
        lo = x; hi = xn;                        % lo inside, hi outside
        for i = 1:40
            mid = 0.5*(lo+hi);
            if fh(mid) <= 0, lo = mid; else, hi = mid; end
        end
        r = lo; return;
    end
    x = xn;
end
end

function [lab, n] = blobs(B)
% 4-connected labelling, written out so no toolbox is needed.
[m, p] = size(B);  lab = zeros(m,p);  n = 0;
for i = 1:m
    for j = 1:p
        if B(i,j) && lab(i,j) == 0
            n = n + 1;  st = [i j];  lab(i,j) = n;
            while ~isempty(st)
                c = st(end,:);  st(end,:) = [];
                for d = [1 0; -1 0; 0 1; 0 -1]'
                    a = c(1)+d(1);  b = c(2)+d(2);
                    if a>=1 && a<=m && b>=1 && b<=p && B(a,b) && lab(a,b)==0
                        lab(a,b) = n;  st(end+1,:) = [a b]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

function ok = ell_ok(Fk, ru, rth, s, tt)
v = Fk(s*ru*cos(tt(:)), s*rth*sin(tt(:)));
ok = all(~isnan(v)) && all(v <= 0);
end
