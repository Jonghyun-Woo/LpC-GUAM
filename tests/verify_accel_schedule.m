% verify_accel_schedule - fly the computed schedule in the real simulator and
% check, at every step, that the aircraft is inside a reachable set that was
% actually computed.
%
% ---------------------------------------------------------------------------
% WHAT IS BEING TESTED
% ---------------------------------------------------------------------------
% optimize_accel_schedule produces 19 accelerations from the LINEAR closed-loop
% model. Everything about that is a prediction: the error propagation, the set
% membership, the margins. This flies the same schedule through the full
% nonlinear simulator - real aero polynomial, real allocation, real actuator
% lags - and asks whether the prediction held.
%
% Three things are checked:
%
%   1. Does the aircraft stay inside a set the whole way? A step counts as
%      covered when it is inside the set of either bracketing anchor, read at
%      the deviation from THAT anchor's own trim point - never a blend of two,
%      which is a set nobody computed.
%
%   2. How far did the simulator end up from the linear prediction? That is
%      the quantity model_eps was supposed to bound.
%
%   3. Does it stay there after the transition ends? The transition is not
%      finished at the last anchor - the loop carries a mode with a 59 s
%      period, so the excursion can still be growing when the ramp stops.
%      T_CRUISE seconds of held cruise are flown and checked too.
%
% The mission descends at WH3 (11.667 ft/s), the condition the trim points and
% their sets were built at. It also starts by SETTLING at that descent: the
% hover trim point is itself a descending condition, so a vehicle at rest is
% 11.7 ft/s away from it and segment 1 would begin already off its anchor.
%
% Output: logger/verify_schedule.mat, logger/verify_schedule.png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

%% 0. setup ------------------------------------------------------------
dt       = 0.01;
WDOT     = 11.667;
T_SETTLE = 45.0;      % acquire the descending hover trim before starting
T_CRUISE = 20.0;      % hold cruise after the last anchor and keep checking
ALT0     = 1500;
SEL      = [4 6 11 8];
LON      = [1 3 5 11];

S  = load(fullfile(root,'logger','accel_schedule.mat'));
a  = S.a;   UH = S.UH;   du_seg = S.du_seg;
T   = load('trim_table_Poly_ConcatVer4p0.mat');
XU0 = T.XU0_interp(:, 1:20, 3);
trim_lon = XU0(LON, :);

Sm  = load(fullfile(root,'logger','closed_loop_model_descending.mat'));
M   = Sm.M;   dXE = Sm.dXE;

spec = FilterConfig.channelSpec('lon');
gv = cell(1,4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
V = cell(1,20);
for k = 1:20
    Sd = load(fullfile(root,'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat',k,3)), 'data');
    V{k} = griddedInterpolant(gv, Sd.data, 'linear', 'nearest');
end

line = @(c) fprintf('%s\n', repmat(c, 1, 92));
fprintf('schedule from logger/accel_schedule.mat, T = %.2f s\n\n', sum(du_seg./a));

%% 1. command profile --------------------------------------------------
v = zeros(1, round(T_SETTLE/dt));                 % settle at the hover trim
seg_of = zeros(size(v));
for k = 1:19
    n  = max(round(du_seg(k)/a(k)/dt), 1);
    v  = [v, UH(k) + (UH(k+1)-UH(k))*(0:n-1)/n];  %#ok<AGROW>
    seg_of = [seg_of, k*ones(1,n)];               %#ok<AGROW>
end
n_cr = round(T_CRUISE/dt);
v  = [v, UH(end)*ones(1, n_cr)];
seg_of = [seg_of, 19*ones(1, n_cr)];
N  = numel(v);
t  = (0:N-1)*dt;
x  = cumtrapz(t, v);
z  = -ALT0 + WDOT*t;

%% 2. fly it -----------------------------------------------------------
params = struct('filter_mode','off');  params.T_seg = 2.0;
guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
st0 = guam.saveState();  st0.state(3) = -ALT0;  guam.restoreState(st0);

XR = zeros(4, N);   ALT = zeros(1, N);
for i = 1:N
    s = guam.state;
    XR(:, i) = [s(4); s(6); s(11); s(8)];
    ALT(i)   = -s(3);
    guam.step(struct('pos',[x(i);0;z(i)], 'vel',[v(i);0;WDOT], 'chi',0,'chi_dot',0));
end

%% 3. the same schedule through the linear model -----------------------
XL = zeros(4, N);   d = zeros(34,1);
for i = 1:N
    [P, xe] = interp_model(M, v(i));
    XL(:, i) = xe(SEL) + d(SEL);
    adot = 0;
    if i < N, adot = (v(i+1)-v(i))/dt; end
    kk = max(min(find(M.u_anchor <= v(i), 1, 'last'), 20), 1);
    d  = P*d - dXE(:, kk)*adot*dt;
end

%% 4. set membership, step by step ------------------------------------
b2h = @(X) [ X(1,:).*cos(X(4,:)) + X(2,:).*sin(X(4,:)) ;
            -X(1,:).*sin(X(4,:)) + X(2,:).*cos(X(4,:)) ;
             X(3,:) ; X(4,:) ];
HR = b2h(XR);   HL = b2h(XL);

Vbest = nan(1, N);   inside = false(1, N);
for i = 1:N
    kb = max(min(find(UH <= v(i), 1, 'last'), 20), 1);
    ks = unique([kb, min(kb+1, 20)]);
    best = inf;
    for kk = ks
        e = HR(:, i) - trim_lon(:, kk);
        if all(e >= spec.grid_min & e <= spec.grid_max)
            best = min(best, V{kk}(e(1), e(2), e(3), e(4)));
        end
    end
    Vbest(i) = best;   inside(i) = best <= 0;
end

%% 5. report -----------------------------------------------------------
tr = seg_of > 0 & t <= T_SETTLE + sum(du_seg./a);
cr = t > T_SETTLE + sum(du_seg./a);

line('='); fprintf('DID IT STAY INSIDE?\n'); line('=');
fprintf('  transition : %6.2f %% of steps inside, worst V = %+.4f\n', ...
        100*mean(inside(tr)), max(Vbest(tr)));
fprintf('  cruise     : %6.2f %% of steps inside, worst V = %+.4f\n', ...
        100*mean(inside(cr)), max(Vbest(cr)));
if all(inside(tr))
    fprintf('  -> covered the whole transition\n');
else
    bad = find(~inside & tr);
    fprintf('  -> LEFT THE SETS: %d steps (%.2f s), first at t = %.2f s, segment %d\n', ...
            numel(bad), numel(bad)*dt, t(bad(1)), seg_of(bad(1)));
    ks = unique(seg_of(bad));
    fprintf('     segments affected: %s\n', mat2str(ks));
end

line('='); fprintf('PREDICTION vs SIMULATOR\n'); line('=');
ER = HR - interp_trim(trim_lon, UH, v);
EL = HL - interp_trim(trim_lon, UH, v);
DIS = abs(EL - ER);
nm = {'u [ft/s]','w [ft/s]','q [deg/s]','th [deg]'};  sc = [1 1 180/pi 180/pi];
fprintf('   channel   | true error max | model error max | model eps allowed\n');
epsm = model_eps(S, mean(a));
for c = 1:4
    fprintf('  %-10s | %13.2f | %14.2f | %14.2f\n', nm{c}, ...
        sc(c)*max(abs(ER(c,tr))), sc(c)*max(DIS(c,tr)), sc(c)*epsm(c));
end

line('='); fprintf('AFTER THE TRANSITION\n'); line('=');
fprintf('  cruise held %.0f s.  final u = %.1f ft/s, w = %.2f, theta = %.2f deg\n', ...
        T_CRUISE, HR(1,end), HR(2,end), rad2deg(HR(4,end)));
fprintf('  altitude %.0f -> %.0f ft (descended %.0f)\n', ALT(1), ALT(end), ALT(1)-ALT(end));
fprintf('  worst V during cruise = %+.4f\n', max(Vbest(cr)));

save(fullfile(root,'logger','verify_schedule.mat'), ...
     't','v','XR','XL','HR','HL','Vbest','inside','seg_of','ALT','a','ER','DIS');

f = figure('Position',[70 70 1150 720], 'Color','w');
subplot(2,2,1); plot(t, v, 'LineWidth',1.6); hold on; plot(t, HR(1,:), 'LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('ft/s'); title('commanded vs flown speed');
legend({'command','flown'},'Location','southeast');
subplot(2,2,2); plot(t, Vbest, 'LineWidth',1.5); hold on; yline(0,'r--','LineWidth',1.4);
grid on; xlabel('t [s]'); ylabel('V'); title('best available set value (V<0 = inside)');
subplot(2,2,3);
plot(t, ER(2,:), 'LineWidth',1.4); hold on; plot(t, EL(2,:), '--', 'LineWidth',1.2);
plot(t, rad2deg(ER(4,:))/10, 'LineWidth',1.4); plot(t, rad2deg(EL(4,:))/10, '--', 'LineWidth',1.2);
grid on; xlabel('t [s]'); title('w [ft/s] and \theta/10 [deg] : solid flown, dashed predicted');
subplot(2,2,4); plot(t, ALT, 'LineWidth',1.6);
grid on; xlabel('t [s]'); ylabel('altitude [ft]'); title('altitude');
saveas(f, fullfile(root,'logger','verify_schedule.png'));
fprintf('\nsaved logger/verify_schedule.mat and logger/verify_schedule.png\n');

%% helpers -------------------------------------------------------------
function [P, xe] = interp_model(M, v)
ua = M.u_anchor;
if v <= ua(1),       k = 1;             al = 0;
elseif v >= ua(end), k = numel(ua)-1;   al = 1;
else
    k = find(ua <= v, 1, 'last');  k = min(k, numel(ua)-1);
    al = (v - ua(k))/(ua(k+1) - ua(k));
end
P  = (1-al)*M.Phi{k}    + al*M.Phi{k+1};
xe = (1-al)*M.xi_e(:,k) + al*M.xi_e(:,k+1);
end

function XT = interp_trim(trim_lon, UH, v)
XT = zeros(4, numel(v));
vv = min(max(v, UH(1)), UH(end));
for c = 1:4
    XT(c,:) = interp1(UH, trim_lon(c,:), vv);
end
end

function e = model_eps(~, a)
A  = [2 3 5 8];
Tb = [0.41 0.71 2.05 5.87; 1.98 3.72 10.45 19.23; ...
      0.0147 0.0227 0.0456 0.0892; 0.0220 0.0545 0.0918 0.1026];
if a <= A(1),       e = Tb(:,1)*(a/A(1));
elseif a >= A(end), e = Tb(:,end);
else,               e = interp1(A, Tb', a)';
end
end
