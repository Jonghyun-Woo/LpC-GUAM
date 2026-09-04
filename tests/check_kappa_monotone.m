% check_kappa_monotone - is V(x, v+k(r-v)) monotone in k?
%
% The line search of Eq. (10) looks for the largest admissible kappa. If V is
% monotone increasing in kappa, bisection is valid and cheap. If not, only a
% grid scan is safe, because the admissible set need not be an interval.
%
% This sweeps kappa at several mission instants and reports whether V is
% monotone, plus where the zero crossing sits.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim  = 0.01;
SAMPLES = [7 9 13 16 20 30];       % mission instants to probe [s]
NK      = 21;                      % kappa grid for the sweep
kk      = linspace(0, 1, NK);

params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
cfg  = Config('trim_schedule', params);
guam = LpC_GUAM(cfg);
rt   = guam.refTraj;
guam.reset();

gov = TrimRefGovernor(traj, brtV, guam, dt_sim, struct('T', 6.0));

Vsweep = nan(numel(SAMPLES), NK);
fprintf('\n  t   | v_prev      r    | monotone? | V(k=0)  V(k=1) | zero crossing\n');
fprintf('%s\n', repmat('-', 1, 80));

i_prev = 0;
for s = 1:numel(SAMPLES)
    ts = SAMPLES(s);
    i_target = round(ts / dt_sim);
    for i = i_prev+1 : i_target
        guam.step(struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                         'chi', rt.chi(i), 'chi_dot', rt.chidot(i)));
    end
    i_prev = i_target;

    % v_prev: pretend the governor has fallen a little behind the schedule
    r      = gov.schedule(ts);
    v_prev = max(r - 6.0, traj.trim_lon(1,1));
    gov.pos_n = rt.pos(1, i_target);

    for j = 1:NK
        v_try = v_prev + kk(j)*(r - v_prev);
        Vsweep(s, j) = gov.predict(v_try, false);   % full value, no early exit
    end

    d    = diff(Vsweep(s, :));
    mono = all(d >= -1e-9);
    ineg = find(Vsweep(s, :) <= 0, 1, 'last');
    if isempty(ineg), zc = NaN; else, zc = kk(ineg); end
    fprintf('%5.1f | %6.2f %6.2f  |    %s    | %+6.3f %+6.3f |  k <= %.2f\n', ...
        ts, v_prev, r, ternary(mono,'YES','NO '), ...
        Vsweep(s,1), Vsweep(s,end), zc);
end

fprintf('\nrollouts used: %d\n', gov.n_rollout);
allmono = all(all(diff(Vsweep, 1, 2) >= -1e-9));
if allmono
    fprintf('V is monotone in kappa at every probe -> bisection is valid.\n');
else
    fprintf(['V is NOT monotone at some probe -> keep the grid scan\n' ...
             '(the admissible set may not be an interval).\n']);
end

%% figure
f = figure('Name','V vs kappa','Position',[60 60 950 600]);
hold on; grid on;
col = turbo(numel(SAMPLES));
for s = 1:numel(SAMPLES)
    plot(kk, Vsweep(s,:), '-o', 'Color', col(s,:), 'LineWidth', 1.5, ...
         'MarkerSize', 3, 'DisplayName', sprintf('t = %g s', SAMPLES(s)));
end
yline(0,'k-','LineWidth',1.2,'HandleVisibility','off');
xlabel('\kappa'); ylabel('V at the candidate command');
title('Is the value monotone in \kappa?  (V \leq 0 = admissible)');
legend('Location','northwest');

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
