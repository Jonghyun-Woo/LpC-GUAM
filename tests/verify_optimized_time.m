% verify_optimized_time - fly the candidate schedules and find the fastest
% one that actually stays inside the safety tubes.
%
% ---------------------------------------------------------------------------
% WHY THIS IS NEEDED
% ---------------------------------------------------------------------------
% optimize_segment_time.m drives the segment time T down to whatever the grid
% allows and reports 8.25 s for the whole transition with no dwell. That is
% the linear model talking past the end of its own validity: Phi was measured
% with small perturbations, so it describes a loop with unlimited control
% power. A ramp steep enough saturates the rotors and the elevator, the real
% speed error stops being c*s, and the placement built on c*s stops meaning
% anything.
%
% The model cannot tell us where that happens. Flying can. This script takes
% each candidate schedule, flies the nonlinear closed loop along it, and
% measures the one thing the whole design is about: does the vehicle stay
% inside the safety tube of the trim point currently being commanded?
%
%   V = BRTValue(state, commanded speed)      V < 0 inside, V >= 0 outside
%
% The fastest schedule that keeps V < 0 for the whole transition IS the answer
% to "the fastest transition that can be guaranteed safe".
%
% Reported per candidate:
%   time      how long the transition actually took
%   V>=0      percentage of the transition spent outside the tube
%   worst V   how far outside it got
%   err       peak forward-speed error, against what the placement assumed
% ---------------------------------------------------------------------------
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt      = 0.01;
T_HOLD  = 5.0;                 % hold at each end of the mission [s]
WH_IDX  = 3;

O = load(fullfile(root,'logger','segment_time_opt.mat'));
u_anchor = O.u_anchor;  a_fwd = O.a_fwd;  a_bwd = O.a_bwd;  e0_u = O.e0_u;
T_GRID = O.T_GRID;  C = O.C;  D_dwell = O.D_dwell;

A_f = @(u) interp1(u_anchor, a_fwd, clampu(u,u_anchor));
A_b = @(u) interp1(u_anchor, a_bwd, clampu(u,u_anchor));
E0  = @(u) interp1(u_anchor, e0_u,  clampu(u,u_anchor));

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj0 = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV  = BRTValue(fullfile(root,'data'), traj0.trim_lon, WH_IDX);

%% ---------------------------------------------------------------------
%  candidates
%% ---------------------------------------------------------------------
% c(u, T) by interpolation on the grid the optimiser built
Cf2 = @(u, T) interp2(u_anchor, T_GRID(:), C, ...
                      clampu(u, u_anchor), min(max(T, T_GRID(1)), T_GRID(end)));

cand = {};
% (0) what is flown today
cand{end+1} = struct('name','uniform 8.44, T=2.0 (today)', ...
                     'us', u_anchor, 'Ts', repmat(2.0,1,numel(u_anchor)-1), ...
                     'tds', zeros(1,numel(u_anchor)-1));
% (1..) the placement at a range of uniform segment times
for Tf = [2.0 1.9 1.75 1.5 1.0 0.5]
    [us, Ts] = walk(u_anchor, A_f, A_b, E0, Cf2, @(u) Tf);
    cand{end+1} = struct('name', sprintf('placed, uniform T=%.2f', Tf), ...
                         'us', us, 'Ts', Ts, 'tds', zeros(1,numel(Ts))); %#ok<SAGROW>
end
% hybrid: keep the slow low-speed end, spend the spare high-speed margin on
% going faster instead of on longer steps. This is the thing a single uniform
% T cannot do.
for sp = [60 80]
  for Thi = [1.5 1.0 0.75]
    Tfun = @(u) 2.0*(u < sp) + Thi*(u >= sp);
    [us, Ts] = walk(u_anchor, A_f, A_b, E0, Cf2, Tfun);
    cand{end+1} = struct('name', sprintf('hybrid T=2.0 below %d, %.2f above', sp, Thi), ...
                         'us', us, 'Ts', Ts, 'tds', zeros(1,numel(Ts))); %#ok<SAGROW>
  end
end

%% ---------------------------------------------------------------------
%  fly them
%% ---------------------------------------------------------------------
fprintf('%-34s | %4s | %5s | %6s | %7s | %6s | %6s | %s\n', ...
        'schedule', 'pts', 'time', 'V>=0 %', 'worst V', 'at u', 'err', 'verdict');
fprintf('%s\n', repmat('-', 1, 100));
res = struct([]);
for i = 1:numel(cand)
    c = cand{i};
    [tt, uu] = profile_from(c.us, c.Ts, c.tds, dt, T_HOLD);
    R = fly(root, params, tt, uu, dt, brtV);
    t_trans = tt(end) - 2*T_HOLD;
    if R.viol == 0, vd = 'SAFE'; else, vd = 'unsafe'; end
    fprintf('%-34s | %4d | %5.1f | %6.1f | %+7.3f | %6.1f | %6.2f | %s\n', ...
            c.name, numel(c.us), t_trans, R.viol, R.Vw, R.u_worst, R.err, vd);
    res(i).name = c.name;  res(i).t = t_trans;  res(i).viol = R.viol;
    res(i).Vw = R.Vw;  res(i).err = R.err;  res(i).V = R.V;
    res(i).t_ax = tt;  res(i).u_ref = uu;  res(i).u_act = R.u_act;
end

ok = find([res.viol] == 0);
fprintf('%s\n', repmat('-', 1, 92));
if isempty(ok)
    fprintf('no candidate stayed inside the tube for the whole transition\n');
else
    [~, j] = min([res(ok).t]);
    b = res(ok(j));
    fprintf('FASTEST SAFE: %s  ->  %.1f s  (worst V = %+.3f)\n', b.name, b.t, b.Vw);
end

%% ---------------------------------------------------------------------
%  figure
%% ---------------------------------------------------------------------
f = figure('Position',[100 100 1000 640]);
subplot(2,1,1); hold on;
for i = 1:numel(res)
    plot(res(i).t_ax, res(i).u_ref, 'LineWidth', 1.1);
end
ylabel('commanded speed [ft/s]'); grid on;
legend({res.name}, 'Location','southeast', 'Interpreter','none', 'FontSize', 7);
title('candidate schedules');
subplot(2,1,2); hold on;
for i = 1:numel(res)
    plot(res(i).t_ax(1:numel(res(i).V)), res(i).V, 'LineWidth', 1.1);
end
yline(0, 'k-', 'LineWidth', 1.4);
xlabel('time [s]'); ylabel('BRT value  (V < 0 = inside)'); grid on;
ylim([-0.5 1.0]);
saveas(f, fullfile(root,'logger','optimized_time_verify.png'));
save(fullfile(root,'logger','optimized_time_verify.mat'), 'res');
fprintf('\nsaved logger/optimized_time_verify.{mat,png}\n');

% =========================================================================
function [us, Ts] = walk(ua, A_f, A_b, E0, Cf2, Tfun)
u = 0;  us = 0;  Ts = [];
while u < ua(end)-1e-6 && numel(us) < 400
    T  = Tfun(u);
    gg = @(D) A_f(u) + A_b(u+D) - D - 2*(E0(u) + Cf2(u,T)*D/T);
    if gg(1e-6) <= 0, error('infeasible at u = %.2f', u); end
    lo = 1e-6; hi = 40;
    if gg(hi) > 0, D = hi; else
        for it = 1:60
            mid = (lo+hi)/2;
            if gg(mid) > 0, lo = mid; else, hi = mid; end
        end
        D = lo;
    end
    D = min(D, ua(end)-u);
    u = u + D;  us(end+1) = u;  Ts(end+1) = T; %#ok<AGROW>
end
end

function [t, u] = profile_from(us, Ts, tds, dt, T_hold)
% piecewise-linear speed profile: hold, then ramp/dwell per segment, then hold
u = repmat(us(1), 1, round(T_hold/dt));
for k = 1:numel(Ts)
    n = max(round(Ts(k)/dt), 1);
    u = [u, us(k) + (us(k+1)-us(k)) * (0:n-1)/n]; %#ok<AGROW>
    nd = round(tds(k)/dt);
    if nd > 0, u = [u, repmat(us(k+1), 1, nd)]; end %#ok<AGROW>
end
u = [u, repmat(us(end), 1, round(T_hold/dt) + 1)];
t = (0:numel(u)-1) * dt;
end

function R = fly(root, params, tt, uu, dt, brtV) %#ok<INUSL>
guam = LpC_GUAM(Config('trim_schedule', params));
guam.reset();
N = numel(tt);
pos = cumtrapz(tt, uu);
V = zeros(1, N);  u_act = zeros(1, N);
n_hold = round(5.0/dt);
for k = 1:N
    ref = struct('pos', [pos(k); 0; -80], 'vel', [uu(k); 0; 0], ...
                 'chi', 0, 'chi_dot', 0);
    u_act(k) = guam.state(4);
    V(k) = brtV.value(guam.state([4 6 11 8]), uu(k));
    guam.step(ref);
end
w = n_hold+1 : N-n_hold;                 % the transition itself
R.V    = V;
R.u_act = u_act;
[R.Vw, iw] = max(V(w));
R.u_worst = uu(w(iw));                   % speed at which it is worst
R.viol = 100 * sum(V(w) >= 0) / numel(w);
R.err  = max(abs(u_act(w) - uu(w)));
end

function y = clampu(u, ua), y = min(max(u, ua(1)), ua(end)); end
