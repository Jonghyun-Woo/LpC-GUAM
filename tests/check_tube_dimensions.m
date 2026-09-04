% check_tube_dimensions - how much of the safety tube does each state axis
% actually consume? This is the "2-D or 4-D?" question, answered by measuring.
%
% ---------------------------------------------------------------------------
% THE QUESTION
% ---------------------------------------------------------------------------
% The placement rule treats the margin as one-dimensional: the tube gives
% a(u) ft/s of room along forward speed, the tracking error e eats 2e of it,
% and what is left is the handover width. That is a 1-D reading of a 4-D
% object - the tube also has extent in vertical velocity w, pitch rate q and
% pitch attitude theta, and the vehicle deviates in those too.
%
% The 1-D reading is only legitimate if the other three axes are nowhere near
% their limits. This script measures whether that is true, by flying the
% schedule that was certified safe and recording, at every instant, how far
% each axis has moved from the trim being commanded, against how far that axis
% is allowed to move before leaving the tube.
%
% READING THE RESULT
%   usage well below 100 % on w, q, theta  ->  the 1-D rule is safe to keep;
%       the speed axis is the binding one and the other three are spectators.
%   usage approaching 100 % on any of them ->  the 1-D rule is optimistic. The
%       tube would have to be eroded in that direction as well, which tightens
%       the placement and lengthens the transition.
% ---------------------------------------------------------------------------
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt     = 0.01;
WH_IDX = 3;
N_TRIM = 20;
T_HOLD = 5.0;
AX     = {'u  (forward speed)', 'w  (vertical speed)', ...
          'q  (pitch rate)   ', 'theta (pitch angle)'};
UNIT   = {'ft/s', 'ft/s', 'rad/s', 'rad'};

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4, gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d)); end
ic = [find(gv{1}==0), find(gv{2}==0), find(gv{3}==0), find(gv{4}==0)];

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj0 = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj0.trim_lon(1,:);
brtV = BRTValue(fullfile(root,'data'), traj0.trim_lon, WH_IDX);

%% ---------------------------------------------------------------------
%  1. tube half-width along each of the four axes, at every anchor
%     (walk out from the tube centre along one axis at a time)
%% ---------------------------------------------------------------------
half = zeros(4, N_TRIM);
for k = 1:N_TRIM
    S = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, WH_IDX)), 'data');
    for d = 1:4
        idx = num2cell(ic);  idx{d} = ':';
        line = squeeze(S.data(idx{:}));
        hp = cross_zero(gv{d}, line, ic(d), +1);
        hm = cross_zero(gv{d}, line, ic(d), -1);
        half(d,k) = min(hp, hm);        % the tighter side
    end
end

fprintf('--- tube half-width per axis (tighter side) ---\n');
fprintf('  u [ft/s] |');  fprintf('  %-8s', 'u', 'w', 'q', 'theta');  fprintf('\n');
for k = 1:N_TRIM
    fprintf('%10.1f |', u_anchor(k));
    fprintf('  %8.3f', half(:,k));
    fprintf('\n');
end

%% ---------------------------------------------------------------------
%  2. fly the certified schedule and record the deviation per axis
%     Schedule: the T = 1.90 placement, which verify_optimized_time found to
%     be the fastest one that stays inside the tube.
%% ---------------------------------------------------------------------
O = load(fullfile(root,'logger','segment_time_opt.mat'));
A_f = @(u) interp1(u_anchor, O.a_fwd, clampu(u,u_anchor));
A_b = @(u) interp1(u_anchor, O.a_bwd, clampu(u,u_anchor));
E0  = @(u) interp1(u_anchor, O.e0_u,  clampu(u,u_anchor));
Cf2 = @(u,T) interp2(u_anchor, O.T_GRID(:), O.C, clampu(u,u_anchor), T);
T_USE = 1.90;

u = 0;  us = 0;
while u < u_anchor(end)-1e-6 && numel(us) < 400
    gg = @(D) A_f(u) + A_b(u+D) - D - 2*(E0(u) + Cf2(u,T_USE)*D/T_USE);
    lo = 1e-6; hi = 40;
    for it = 1:60
        mid = (lo+hi)/2;
        if gg(mid) > 0, lo = mid; else, hi = mid; end
    end
    u = u + min(lo, u_anchor(end)-u);  us(end+1) = u; %#ok<SAGROW>
end

uu = repmat(us(1), 1, round(T_HOLD/dt));
for k = 1:numel(us)-1
    n = round(T_USE/dt);
    uu = [uu, us(k) + (us(k+1)-us(k))*(0:n-1)/n]; %#ok<AGROW>
end
uu = [uu, repmat(us(end), 1, round(T_HOLD/dt)+1)];
tt = (0:numel(uu)-1)*dt;
pos = cumtrapz(tt, uu);

guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
N = numel(tt);  dev = zeros(4, N);  V = zeros(1, N);
for k = 1:N
    x   = guam.state([4 6 11 8]);            % u, w, q, theta
    xe  = brtV.trim_at(uu(k));               % the trim being commanded
    dev(:,k) = x(:) - xe(:);
    V(k) = brtV.value(x, uu(k));
    guam.step(struct('pos',[pos(k);0;-80], 'vel',[uu(k);0;0], ...
                     'chi',0, 'chi_dot',0));
end
w_ix = round(T_HOLD/dt)+1 : N-round(T_HOLD/dt);

%% ---------------------------------------------------------------------
%  3. how much of each axis was used
%% ---------------------------------------------------------------------
fprintf('\n--- axis usage on the certified T = %.2f schedule (%.1f s) ---\n', ...
        T_USE, tt(end)-2*T_HOLD);
fprintf('%-20s | %10s | %10s | %8s | %s\n', ...
        'axis', 'peak dev', 'allowed', 'usage %', 'unit');
fprintf('%s\n', repmat('-', 1, 70));
peak = zeros(1,4);  allow = zeros(1,4);
for d = 1:4
    [peak(d), ip] = max(abs(dev(d, w_ix)));
    allow(d) = interp1(u_anchor, half(d,:), clampu(uu(w_ix(ip)), u_anchor));
    fprintf('%-20s | %10.4f | %10.4f | %7.1f%% | %s\n', ...
            AX{d}, peak(d), allow(d), 100*peak(d)/allow(d), UNIT{d});
end
fprintf('\nworst BRT value over the transition: %+.4f  (must stay < 0)\n', ...
        max(V(w_ix)));

fprintf(['\nInterpretation\n' ...
    '  The speed axis is the one the placement rule accounts for. Any other\n' ...
    '  axis whose usage approaches 100%% is margin the 1-D rule is spending\n' ...
    '  without counting it, and would have to be subtracted from the tube\n' ...
    '  before the overlap test - making the transition longer, not shorter.\n']);

%% ---------------------------------------------------------------------
%  4. figure
%% ---------------------------------------------------------------------
f = figure('Position',[100 100 1000 720]);
for d = 1:4
    subplot(4,1,d);
    plot(tt, dev(d,:), 'b-', 'LineWidth', 1.1); hold on;
    hb = interp1(u_anchor, half(d,:), clampu(uu, u_anchor));
    plot(tt,  hb, 'r--', tt, -hb, 'r--', 'LineWidth', 1.0);
    ylabel(sprintf('%s [%s]', strtrim(AX{d}), UNIT{d}), 'FontSize', 8);
    grid on;
    if d == 1, title(sprintf('deviation from the commanded trim vs tube half-width (T = %.2f s)', T_USE)); end
    if d == 4, xlabel('time [s]'); end
end
saveas(f, fullfile(root,'logger','tube_dimensions.png'));
save(fullfile(root,'logger','tube_dimensions.mat'), ...
     'half','dev','V','tt','uu','us','peak','allow','T_USE');
fprintf('\nsaved logger/tube_dimensions.{mat,png}\n');

% =========================================================================
function d = cross_zero(g, line, i0, sgn)
n = numel(g);  i = i0;
while true
    j = i + sgn;
    if j < 1 || j > n, d = abs(g(i)-g(i0)); return; end
    if line(j) >= 0
        f = -line(i)/(line(j)-line(i));
        d = abs(g(i) + f*(g(j)-g(i)) - g(i0));  return;
    end
    i = j;
end
end

function y = clampu(u, ua), y = min(max(u, ua(1)), ua(end)); end
