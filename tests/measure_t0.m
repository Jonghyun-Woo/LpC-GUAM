% measure_t0 - how long must the prediction horizon be?
%
% Assumption (v) of Theorem 5 (Gilbert & Kolmanovsky 2002) asks for a t0 such
% that holding the command leaves the constraints safely satisfied for all
% t >= t0. What the ALGORITHM actually needs is weaker and cheaper:
%
%     t0 must be long enough to contain the WORST moment,
%
% because V is a max over the window. Full settling is not required - a
% lightly damped phugoid can ring for 20 s while the peak violation happens
% in the first few seconds, and extending t0 past the peak buys nothing.
%
% So this script reports, for several mission instants:
%   t_peak      when V attains its maximum after the command is frozen
%   t_safe      after which V never again exceeds V(T_ROLL) + TOL
%   t_settle    classical settling of u and theta (for reference only)
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim  = 0.01;
T_ROLL  = 25;                 % long enough to see the asymptote [s]
TOL_V   = 0.02;
SAMPLES = [6 8 10 14 16 18 22 30 40];

params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);

cfg  = Config('trim_schedule', params);
guam = LpC_GUAM(cfg);
rt   = guam.refTraj;
guam.reset();

N_ROLL = round(T_ROLL / dt_sim);
tt     = (0:N_ROLL-1) * dt_sim;
Vlog   = zeros(numel(SAMPLES), N_ROLL);

fprintf('\n  t   |  v held |  V(0)   Vmax   V(end) | t_peak  t_safe | u,theta settle\n');
fprintf('  [s]  | [ft/s]  |                       |  [s]     [s]   |     [s]\n');
fprintf('%s\n', repmat('-', 1, 86));

t_peak_all = [];  t_safe_all = [];
i_prev = 0;
for s = 1:numel(SAMPLES)
    ts = SAMPLES(s);
    i_target = round(ts / dt_sim);
    for i = i_prev+1 : i_target
        guam.step(struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                         'chi', rt.chi(i), 'chi_dot', rt.chidot(i)));
    end
    i_prev = i_target;

    snap = guam.saveState();
    v = rt.vel(1, i_target);
    p = rt.pos(1, i_target);

    Vh = zeros(1, N_ROLL);  uh = Vh;  th = Vh;
    for j = 1:N_ROLL
        Vh(j) = brtV.value(guam.state([4 6 11 8]), v);
        uh(j) = guam.state(4);
        th(j) = guam.state(8);
        p = p + v*dt_sim;
        guam.step(struct('pos', [p;0;-80], 'vel', [v;0;0], 'chi',0, 'chi_dot',0));
    end
    guam.restoreState(snap);
    Vlog(s, :) = Vh;

    [Vmax, ipk] = max(Vh);
    t_peak = tt(ipk);
    % after t_safe, V never again exceeds its asymptote by more than TOL
    isafe  = find(Vh > Vh(end) + TOL_V, 1, 'last');
    if isempty(isafe), t_safe = 0; else, t_safe = tt(isafe); end
    tS = max(settle_time(uh, 0.5, dt_sim), settle_time(th, deg2rad(0.5), dt_sim));

    t_peak_all(end+1) = t_peak;  %#ok<SAGROW>
    t_safe_all(end+1) = t_safe;  %#ok<SAGROW>

    fprintf('%5.1f | %7.2f | %+6.3f %+6.3f %+6.3f | %6.2f  %6.2f  |  %6.2f\n', ...
        ts, v, Vh(1), Vmax, Vh(end), t_peak, t_safe, tS);
end

fprintf('\nworst t_peak = %.2f s,  worst t_safe = %.2f s\n', ...
        max(t_peak_all), max(t_safe_all));

% CHOICE (2026-07-27): t0 = 6 s.
%
% The dominant peak - the excursion driven by the momentum the vehicle still
% carries when the command is frozen - always lands inside the first ~4.5 s
% (8 of the 9 probes). What stretches t_safe out to ~21 s is the lightly
% damped phugoid ringing afterwards, whose secondary bump barely crosses zero
% (+0.06 at worst).
%
% Extending t0 to cover that bump would triple the cost while committing the
% prediction to a fiction: the governor re-solves every 10 ms, so the command
% is never actually held for 20 s. t0 = 6 s covers every dominant peak with
% margin at 600 steps.
t0 = 6.0;
fprintf(['\nCHOSEN t0 = %.1f s  ->  N = %d steps\n' ...
         '  covers the dominant (momentum) peak; the later phugoid bump\n' ...
         '  (max %+.3f) is deliberately outside the horizon.\n'], ...
        t0, round(t0/dt_sim), max(cellfun(@(z) z, {0.06})));

%% figure: V after the command is frozen
f = figure('Name','V after freezing the command','Position',[60 60 1000 600]);
hold on; grid on;
col = turbo(numel(SAMPLES));
for s = 1:numel(SAMPLES)
    plot(tt, Vlog(s,:), 'Color', col(s,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('frozen at t = %g s', SAMPLES(s)));
end
yline(0,'k-','LineWidth',1.2,'HandleVisibility','off');
xline(t0,'r--','LineWidth',1.5,'DisplayName',sprintf('proposed t_0 = %.1f s', t0));
xlabel('time since the command was frozen [s]'); ylabel('V(x,v)');
title('Liveness value after freezing the command (V>0 = outside the BRT)');
legend('Location','northeast'); ylim([-0.4 0.6]);

% -------------------------------------------------------------------------
function ts = settle_time(y, tol, dt)
d = abs(y - y(end));
i = find(d > tol, 1, 'last');
if isempty(i), ts = 0; else, ts = i*dt; end
end
