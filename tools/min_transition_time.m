% min_transition_time - fastest safe transition from the reachable sets alone.
%
% Combines what the previous two tools established:
%   tools/brt_inspect.m         array layout, sign convention, deviation frame
%   tools/extract_brt_section.m the (u, theta) sections and the handover gate
%
% and adds the missing piece: as the command ramps through a segment the
% aircraft does not sit at the trim point, it sits at a pitch offset that is
% proportional to the acceleration. Push the acceleration up and the pitch
% offset walks the aircraft out of the sets. The largest acceleration that
% keeps a valid handover is read straight off the grid.
%
%   theta_dev(u) = [theta_trim(u) - theta_trim(u_k)] + mu_theta(u) * a
%                   \___ the trim schedule moves ___/  \_ accel offset _/
%
% ACTUATOR LIMITS ARE DELIBERATELY ABSENT. This is the reachable-set answer
% only. When the actuator branch produces its own a_max(u), combine with
% min(a_star_BRT, a_max_actuator) - do not fold it in here.
%
% *** All distances scale with the grid extent assumed by
% *** FilterConfig.channelSpec('lon'), which the BRT files do not record.
%
% Reads only. Output: data/min_transition_time_WH3.mat, logger/min_transition_*.png

clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));
dataDir = fullfile(root,'data');

K       = 20;
WH_IDX  = 3;
IW0     = 21;                 % w = 0
IQ0     = 31;                 % q = 0
DU_SCAN = 0.05;               % [ft/s] speed sweep
A_SCAN  = 0:0.05:20;          % [ft/s^2] exhaustive acceleration sweep
ETAS    = [1.0 0.9 0.8 0.7];
ETA_DEF = 0.8;
A_REF   = [3.38 4.69 5.63];   % the band-schedule accelerations
G       = 32.174;
LIN_LIM = 0.20;               % |pitch offset| beyond this leaves linear validity
T_BAND  = [repmat(2.5,1,4) repmat(1.8,1,4) repmat(1.5,1,11)];  % the 33.7 s schedule

line = @(c) fprintf('%s\n', repmat(c,1,96));
line('='); fprintf('GRID EXTENT IS AN ASSUMPTION (FilterConfig), NOT VERIFIED FROM THE BRT FILES\n');
fprintf('ACTUATOR LIMITS ARE NOT INCLUDED - this is the reachable-set answer only\n'); line('=');

%% ========================================================================
%  SETUP - sections, interpolants, trim pitch schedule
%  The interpolants cannot be stored in a .mat, so they are rebuilt exactly
%  as extract_brt_section.m builds them.
%% ========================================================================
spec = FilterConfig.channelSpec('lon');
gu  = linspace(spec.grid_min(1), spec.grid_max(1), spec.grid_num(1));
gth = linspace(spec.grid_min(4), spec.grid_max(4), spec.grid_num(4));

T  = load('trim_table_Poly_ConcatVer4p0.mat');
UH = T.UH(1:K);
th_trim = squeeze(T.XU0_interp(11, 1:K, WH_IDX))';    % 20x1 [rad]
du_seg  = diff(UH);
NS = K-1;

F = cell(1,K);
for k = 1:K
    S = load(fullfile(dataDir, sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    F{k} = griddedInterpolant({gu, gth}, squeeze(S.data(:, IW0, IQ0, :)), 'linear', 'none');
end

B = load(fullfile(dataDir,'brt_section_WH3.mat'));    % PART B results
fprintf('loaded %d sections and data/brt_section_WH3.mat\n\n', K);

%% ========================================================================
%  STEP 1 - pitch sensitivity of the shipped allocation
%  M = W^-1 B' (B W^-1 B')^-1; row 12 is the virtual pitch effector, so
%  M(12,1) is the pitch the allocation spends per unit forward acceleration.
%% ========================================================================
mu_theta = zeros(K,1);
for k = 1:K
    W = T.W_lon_interp(:,:,k,WH_IDX);
    Bl = T.B_lon_interp(:,:,k,WH_IDX);
    M = (W \ Bl') / (Bl * (W \ Bl'));
    mu_theta(k) = M(12,1);
end
fprintf('--- STEP 1: mu_theta ---\n');
fprintf('  hover (UH1)   : %+.6f rad per ft/s^2\n', mu_theta(1));
fprintf('  -1/g          : %+.6f\n', -1/G);
fprintf('  difference    : %+.6f  (%.2f %%)  -> %s\n', ...
        mu_theta(1)+1/G, 100*abs(mu_theta(1)+1/G)*G, ...
        ternary(abs(mu_theta(1)+1/G) < 0.003, 'consistent with -1/g', 'NOT close to -1/g'));
fprintf('  cruise (UH20) : %+.6f   (implied tilt share lambda = %.3f)\n\n', ...
        mu_theta(K), abs(mu_theta(K))*G);

%% ========================================================================
%  STEP 2-4 - S_k(a) by exhaustive sweep, no monotonicity assumed
%% ========================================================================
NA = numel(A_SCAN);
S_curve = -inf(NS, NA);
sup_curve = nan(NS, NA);  inf_curve = nan(NS, NA);
a_star = zeros(NS,1);
nonmono = false(NS,1);
hit_ceiling = false(NS,1);

for k = 1:NS
    ua = UH(k) : DU_SCAN : UH(k+1);
    if ua(end) < UH(k+1) - 1e-9, ua(end+1) = UH(k+1); end %#ok<SAGROW>
    n = numel(ua);
    % trim pitch and pitch sensitivity, linearly interpolated across the segment
    tht = interp1(UH([k k+1]), th_trim([k k+1]),  ua(:));
    mut = interp1(UH([k k+1]), mu_theta([k k+1]), ua(:));
    d_k   = ua(:) - UH(k);           % 0 .. du   in set k's frame
    d_kp1 = ua(:) - UH(k+1);         % -du .. 0  in set k+1's frame
    sh_k   = tht - th_trim(k);       % trim-schedule shift, set k's frame
    sh_kp1 = tht - th_trim(k+1);

    for ia = 1:NA
        a = A_SCAN(ia);
        vk  = F{k  }(d_k,   sh_k   + mut*a);
        vk1 = F{k+1}(d_kp1, sh_kp1 + mut*a);
        in_k   = ~isnan(vk)  & vk  <= 0;
        in_kp1 = ~isnan(vk1) & vk1 <= 0;

        % sup of the run that starts at u_k and never breaks
        if ~in_k(1)
            continue;                        % not inside at the start -> S = -inf
        end
        b = find(~in_k, 1, 'first');
        if isempty(b), i_sup = n; else, i_sup = b-1; end

        % inf of the run that ends at u_{k+1} and never breaks going back
        if ~in_kp1(end)
            continue;
        end
        b2 = find(~in_kp1, 1, 'last');
        if isempty(b2), i_inf = 1; else, i_inf = b2+1; end

        sup_curve(k,ia) = ua(i_sup);
        inf_curve(k,ia) = ua(i_inf);
        S_curve(k,ia)   = ua(i_sup) - ua(i_inf);
    end

    ok = S_curve(k,:) >= 0;
    if any(ok)
        ia_last = find(ok, 1, 'last');
        a_star(k) = A_SCAN(ia_last);
        hit_ceiling(k) = (ia_last == NA);
        % non-monotone if feasibility returns after a gap
        nonmono(k) = ~all(ok(1:ia_last));
    else
        a_star(k) = 0;
    end
end

if any(nonmono)
    fprintf(['** WARNING: S_k(a) is NOT monotone on segments %s - feasibility comes\n' ...
             '   back after being lost. a_star is the LARGEST feasible a, so there is\n' ...
             '   an infeasible pocket below it. Treat those segments with care.\n\n'], ...
             mat2str(find(nonmono)'));
end
if any(hit_ceiling)
    fprintf('** WARNING: a_star hit the sweep ceiling (%.1f ft/s^2) on segments %s.\n\n', ...
            A_SCAN(end), mat2str(find(hit_ceiling)'));
end

%% ========================================================================
%  STEP 5 - assemble the schedule
%% ========================================================================
NE = numel(ETAS);
Tseg = nan(NS, NE);  Ttot = nan(1, NE);
for e = 1:NE
    a = ETAS(e) * a_star;
    Tseg(:,e) = du_seg ./ max(a, eps);
    Ttot(e) = sum(Tseg(:,e));
end
ie = find(ETAS == ETA_DEF, 1);
t_arrive = [0; cumsum(Tseg(:,ie))];

%% ========================================================================
%  TABLE 1
%% ========================================================================
line('='); fprintf('TABLE 1 - per segment (eta = %.1f)\n', ETA_DEF); line('=');
fprintf('%3s | %7s | %8s | %10s | %8s | %9s | %s\n', ...
        'k','du','a_star','a(eta)','T [s]','t [s]','note');
line('-');
for k = 1:NS
    nt = '';
    if nonmono(k),     nt = [nt 'non-monotone ']; end %#ok<AGROW>
    if hit_ceiling(k), nt = [nt 'ceiling ']; end %#ok<AGROW>
    if a_star(k) == 0, nt = [nt 'INFEASIBLE at any a ']; end %#ok<AGROW>
    fprintf('%3d | %7.3f | %8.2f | %10.2f | %8.2f | %9.2f | %s\n', ...
            k, du_seg(k), a_star(k), ETAS(ie)*a_star(k), Tseg(k,ie), t_arrive(k+1), nt);
end
fprintf('%3s | %7.3f | %8s | %10s | %8.2f |\n', 'sum', sum(du_seg), '', '', Ttot(ie));

%% ========================================================================
%  TABLE 2
%% ========================================================================
line('='); fprintf('TABLE 2 - total transition time vs safety factor\n'); line('=');
fprintf('%6s | %12s | %s\n', 'eta', 'T_total [s]', 'vs the 33.7 s band schedule');
line('-');
for e = 1:NE
    fprintf('%6.1f | %12.2f | %+.2f s (%+.1f %%)\n', ETAS(e), Ttot(e), ...
            Ttot(e)-33.7, 100*(Ttot(e)/33.7 - 1));
end

%% ========================================================================
%  TABLE 3
%% ========================================================================
line('='); fprintf('TABLE 3 - pitch sensitivity and headroom against the band accelerations\n'); line('=');
fprintf('%3s | %8s | %11s | %10s | %10s | %10s | %s\n', ...
        'k','UH','mu_theta','a*/3.38','a*/4.69','a*/5.63','pitch at a* [rad]');
line('-');
for k = 1:NS
    fprintf('%3d | %8.1f | %+11.6f | %10.2f | %10.2f | %10.2f | %+8.4f\n', ...
            k, UH(k), mu_theta(k), a_star(k)/A_REF(1), a_star(k)/A_REF(2), ...
            a_star(k)/A_REF(3), mu_theta(k)*a_star(k));
end
fprintf('\nratios below 1.0 mean the reachable sets forbid that band acceleration.\n');

%% ========================================================================
%  CHECK 1 - S_k(0) against the PART B overlap
%% ========================================================================
line('='); fprintf('CHECK 1 - S_k(a=0) vs PART B overlap_width(theta=0)\n'); line('=');
fprintf('%3s | %14s | %18s | %9s | %s\n', 'k','S_k(0)','PART B overlap','diff','note');
line('-');
d1 = nan(NS,1);
for k = 1:NS
    d1(k) = S_curve(k,1) - B.overlap_width_theta0(k);
    fprintf('%3d | %14.3f | %18.3f | %+9.3f |\n', k, S_curve(k,1), B.overlap_width_theta0(k), d1(k));
end
fprintf(['\nmax |diff| = %.3f ft/s^2-free units -> %s\n'], max(abs(d1)), ...
        ternary(max(abs(d1)) < 0.3, 'consistent', 'LARGE - check the frame conversion'));
fprintf(['The two are not identical by construction: PART B probes at exactly\n' ...
         'theta = 0 in both frames, while S_k(0) carries the trim-schedule shift\n' ...
         '(theta_trim(u) - theta_trim(u_k)), which is up to %.4f rad on a segment.\n'], ...
        max(abs(diff(th_trim))));

%% ========================================================================
%  CHECK 2 - pitch offset at a_star against the PART B band
%% ========================================================================
line('='); fprintf('CHECK 2 - pitch offset at a_star inside the PART B band?\n'); line('=');
fprintf('%3s | %14s | %19s | %s\n', 'k','mu*a_star','PART B band','verdict');
line('-');
n_out = 0;
for k = 1:NS
    p = mu_theta(k)*a_star(k);
    inb = p >= B.theta_lo(k) && p <= B.theta_hi(k);
    if ~inb, n_out = n_out + 1; end
    fprintf('%3d | %14.4f | %+8.3f .. %+7.3f | %s\n', k, p, ...
            B.theta_lo(k), B.theta_hi(k), ternary(inb,'inside','OUTSIDE'));
end
fprintf('\noutside the band: %d of %d\n', n_out, NS);
fprintf(['The band is a theta-sweep at a fixed pitch, so it is a reference, not a\n' ...
         'hard bound: the trajectory here varies pitch along the segment.\n']);

%% ========================================================================
%  CHECK 3 - where the time is won or lost against 33.7 s
%% ========================================================================
line('='); fprintf('CHECK 3 - segment-by-segment contribution vs the 33.7 s band schedule\n'); line('=');
fprintf('%3s | %10s | %12s | %10s | %s\n', 'k','ours [s]','band [s]','diff [s]','');
line('-');
dT = Tseg(:,ie) - T_BAND(:);
for k = 1:NS
    mark = '';
    if dT(k) >  0.3, mark = '<- slower (bottleneck)'; end
    if dT(k) < -0.3, mark = '<- faster (gain)'; end
    fprintf('%3d | %10.2f | %12.2f | %+10.2f | %s\n', k, Tseg(k,ie), T_BAND(k), dT(k), mark);
end
fprintf('%3s | %10.2f | %12.2f | %+10.2f |\n', 'sum', Ttot(ie), sum(T_BAND), sum(dT));
[~, iw] = sort(dT, 'descend');
fprintf('\nlargest three losses : segments %s\n', mat2str(iw(1:3)'));
fprintf('largest three gains  : segments %s\n', mat2str(iw(end:-1:end-2)'));

%% ========================================================================
%  CHECK 4 - linear-validity flag
%% ========================================================================
line('='); fprintf('CHECK 4 - pitch offset beyond the linear range (|mu*a*| > %.2f rad)\n', LIN_LIM); line('=');
bad = find(abs(mu_theta(1:NS).*a_star) > LIN_LIM);
if isempty(bad)
    fprintf('  none. every segment stays inside the linear range.\n');
else
    for k = bad(:)'
        fprintf('  segment %2d : mu*a* = %+.4f rad (%.1f deg) - a_star here rests on\n', ...
                k, mu_theta(k)*a_star(k), rad2deg(mu_theta(k)*a_star(k)));
        fprintf('               a pitch the linearised mu_theta was never checked at.\n');
    end
end

%% ========================================================================
%  FIGURES
%% ========================================================================
f1 = figure('Position',[60 60 1000 620]); hold on;
cmap = parula(NS);
for k = 1:NS
    s = S_curve(k,:);  s(isinf(s)) = NaN;
    plot(A_SCAN, s, '-', 'Color', cmap(k,:), 'LineWidth', 1.1);
end
yline(0,'k-','LineWidth',1.8);
for p = 1:numel(A_REF), xline(A_REF(p), 'k--', sprintf('%.2f', A_REF(p))); end
xlabel('acceleration a [ft/s^2]'); ylabel('S_k(a) = sup I_k - inf I_{k+1}  [ft/s]');
title('handover margin vs acceleration (S \geq 0 required); colour = segment 1..19');
colormap(cmap); cb = colorbar; cb.Label.String = 'segment';
caxis([1 NS]); grid on;
saveas(f1, fullfile(root,'logger','min_transition_S_curves.png'));

f2 = figure('Position',[60 60 1000 520]); hold on;
bar(1:NS, a_star, 'FaceColor', [0.75 0.83 0.93]);
plot(1:NS, a_star, 'b-o', 'LineWidth', 1.5);
for p = 1:numel(A_REF)
    yline(A_REF(p), '--', sprintf('band a = %.2f', A_REF(p)), 'LineWidth', 1.2);
end
xlabel('segment k'); ylabel('a^*(k) [ft/s^2]');
title('largest acceleration the reachable sets allow (no actuator limits)');
grid on; xticks(1:NS);
saveas(f2, fullfile(root,'logger','min_transition_astar.png'));

f3 = figure('Position',[60 60 1000 520]); hold on;
tb = [0 cumsum(T_BAND)];
stairs(tb, UH, 'k-', 'LineWidth', 1.6);
plot(tb, UH, 'ko', 'MarkerSize', 4);
for e = 1:NE
    te = [0; cumsum(Tseg(:,e))];
    plot(te, UH, '-', 'LineWidth', 1.4);
end
xlabel('time [s]'); ylabel('commanded speed [ft/s]');
title('speed profile: band schedule vs reachable-set answer at each safety factor');
legend([{sprintf('band schedule (%.1f s)', sum(T_BAND))}, ...
        arrayfun(@(e) sprintf('\\eta = %.1f (%.1f s)', ETAS(e), Ttot(e)), 1:NE, ...
                 'UniformOutput', false)], 'Location','southeast');
grid on;
saveas(f3, fullfile(root,'logger','min_transition_profile.png'));

%% ========================================================================
%  SAVE
%% ========================================================================
meta = struct('WH_IDX', WH_IDX, 'UH', UH, 'a_scan', A_SCAN, 'du_scan', DU_SCAN, ...
              'etas', ETAS, 'eta_default', ETA_DEF, 'T_band', T_BAND, ...
              'actuator_limits_included', false, 'grid_extent_assumed', true, ...
              'note', ['reachable-set answer only; combine with the actuator ' ...
                       'branch as min(a_star, a_max_actuator). Distances scale ' ...
                       'with the assumed FilterConfig grid extent.']);
out = fullfile(dataDir,'min_transition_time_WH3.mat');
save(out, 'UH','du_seg','mu_theta','th_trim','a_star','S_curve','sup_curve', ...
     'inf_curve','A_SCAN','Tseg','Ttot','ETAS','t_arrive','nonmono', ...
     'hit_ceiling','dT','meta');

fprintf('\nsaved %s\n', out);
fprintf('saved logger/min_transition_{S_curves,astar,profile}.png\n');
line('='); fprintf('REMINDER: no actuator limits here. Combine with the actuator branch by min().\n'); line('=');

%% ------------------------------------------------------------------------
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
