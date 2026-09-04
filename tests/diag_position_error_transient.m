% diag_position_error_transient - why a command that runs ahead of the
% vehicle produces the overspeed / pitch-down transient.
%
% RSLQR does not track the reference velocity directly. Its outer loop adds
% the body-frame position error to the velocity command with gain 0.1
% (RSLQR.m:227):
%
%   Vb_cmd = ref.vel - trim + 0.1 * e_pos_body,   e_pos_body = R_i2b*(p_ref - p)
%
% So the speed the inner loop actually chases is
%
%   u_effective = ref.vel(1) + 0.1 * e_pos_body(1)
%
% If the commanded point runs ahead, the vehicle lags, the position error
% integrates, and this term silently inflates the speed command.
%
% MEASURED RESULT (2026-07-27) - the proximity constraint does NOT fix this:
%
%   config                       max|e_pos|   max(u-u_ref)   min theta
%   ungoverned                     38.3 ft      +11.30        -9.49 deg
%   option 1 (no proximity)        36.3 ft      +11.69        -8.73 deg
%   option 2 (proximity s_f=6)     35.9 ft      +11.64        -8.60 deg
%
% Proximity caps the INSTANTANEOUS speed gap (8 -> 6 ft/s, visible in
% panel 3 around t = 6..10 s) exactly as designed, but the position error
% is the INTEGRAL of that gap over ~8 s of ramp, so a 2 ft/s cap removes
% only ~6 % of it. The transient's root cause is that the vehicle simply
% cannot match the commanded acceleration out of hover - not that the
% command runs away. Fixing it needs a slower low-speed T_seg or outer-loop
% anti-windup, not a tighter proximity bound.
%
% This script shows the chain for three configurations.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
traj = TrimScheduleTrajectory.build(struct('dt', dt_sim));
spec = FilterConfig.channelSpec('lon');
gv = cell(1, 4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
brt = cell(1, traj.n_trim);
for k = 1:traj.n_trim
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, traj.wh_idx)), 'data');
    brt{k} = S.data;
end
cfg = Config('trim_schedule', struct('filter_mode', 'off'));

R(1) = run_ungoverned(cfg);
R(2) = run_governed(traj, brt, gv, cfg, dt_sim, 'opt1');
R(3) = run_governed(traj, brt, gv, cfg, dt_sim, 'opt2');
names = {'ungoverned (no proximity)', 'option 1 (no proximity)', ...
         'option 2 (proximity s_f=6)'};

%% Numbers at the worst moment
fprintf('\n%-28s | max |e_pos| | that adds | max u over | min theta\n', 'configuration');
fprintf('%-28s |    [ft]     | [ft/s] to  | schedule   |  [deg]\n', '');
fprintf('%-28s |             | the cmd    |  [ft/s]    |\n', '');
fprintf('%s\n', repmat('-', 1, 86));
for c = 1:3
    [me, ie] = max(abs(R(c).e_pos));
    fprintf('%-28s |   %6.1f    |   %+5.2f    |   %+6.2f   |  %+6.2f   (t=%.1f s)\n', ...
        names{c}, me, 0.1 * R(c).e_pos(ie), max(R(c).u - R(c).u_ref), ...
        rad2deg(min(R(c).th)), R(c).t(ie));
end
fprintf(['\n"that adds" = 0.1*e_pos, the speed the outer loop silently adds\n' ...
         'to the schedule command (RSLQR.m:227).\n']);

%% Operational cost, measured in ABSOLUTE terms
% u - u_ref above is measured against a reference that the governor itself
% freezes, so it exaggerates. These are frame-independent quantities.
fprintf('\n%-28s | theta range  | peak |q|  | altitude dev | mission\n', 'configuration');
fprintf('%-28s |    [deg]     |  [deg/s]  |    [ft]      |   [s]\n', '');
fprintf('%s\n', repmat('-', 1, 86));
for c = 1:3
    fprintf('%-28s | %+5.2f .. %+5.2f | %8.2f  |   %6.2f     | %6.2f\n', ...
        names{c}, rad2deg(min(R(c).th)), rad2deg(max(R(c).th)), ...
        rad2deg(max(abs(R(c).q))), max(abs(R(c).alt - 80)), R(c).t(end));
end

%% Figure: the causal chain, top to bottom
f = figure('Name', 'position error -> command inflation -> transient', ...
           'Position', [50 30 1250 950]);
col = [0.85 0.1 0.1; 0.2 0.4 0.9; 0.1 0.6 0.2];

ax(1) = subplot(4, 1, 1); hold on; grid on;
for c = 1:3
    plot(R(c).t, R(c).e_pos, '-', 'Color', col(c, :), 'LineWidth', 1.5, ...
         'DisplayName', names{c});
end
yline(0, 'k-', 'HandleVisibility', 'off');
ylabel('e_{pos} [ft]');
title('1) Position error: how far BEHIND the commanded point the vehicle is');
legend('Location', 'northwest');

ax(2) = subplot(4, 1, 2); hold on; grid on;
for c = 1:3
    plot(R(c).t, 0.1 * R(c).e_pos, '-', 'Color', col(c, :), 'LineWidth', 1.5);
end
yline(0, 'k-');
ylabel('0.1\cdot e_{pos} [ft/s]');
title('2) ... which the outer loop ADDS to the speed command (gain 0.1)');

ax(3) = subplot(4, 1, 3); hold on; grid on;
for c = 1:3
    plot(R(c).t, R(c).u - R(c).u_ref, '-', 'Color', col(c, :), 'LineWidth', 1.5);
end
yline(0, 'k-');
ylabel('u - u_{ref} [ft/s]');
title('3) ... so the vehicle ends up FASTER than the schedule asked for');

ax(4) = subplot(4, 1, 4); hold on; grid on;
for c = 1:3
    plot(R(c).t, rad2deg(R(c).th), '-', 'Color', col(c, :), 'LineWidth', 1.5);
end
yline(0, 'k-');
ylabel('\theta [deg]'); xlabel('time [s]');
title('4) ... and pitches DOWN to make that speed - the (u up, \theta down) BRT exit');
linkaxes(ax, 'x');
xlim(ax(1), [0, 30]);

% -------------------------------------------------------------------------
function S = run_ungoverned(cfg)
guam = LpC_GUAM(cfg);
rt   = guam.refTraj;  M = size(rt.pos, 2);
S = blank(M);
guam.reset();
for i = 1:M
    ref = struct('pos', rt.pos(:, i), 'vel', rt.vel(:, i), ...
                 'chi', rt.chi(i), 'chi_dot', rt.chidot(i));
    S = record(S, i, rt.time(i), guam.state, ref);
    guam.step(ref);
end
S = trim_log(S, M);
end

function S = run_governed(traj, brt, gv, cfg, dt_sim, which)
T_max = 80;  Ng = round(T_max / dt_sim);
guam  = LpC_GUAM(cfg);
if strcmp(which, 'opt1')
    gov = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05);
else
    gov = TrimPointGovernor(traj, brt, gv, dt_sim, 6.0, 1.0, 0.05);
end
guam.reset();  gov.reset();
S = blank(Ng);
for i = 1:Ng
    tt = (i - 1) * dt_sim;
    [ref, gi] = gov.step(guam.state, tt);
    S = record(S, i, tt, guam.state, ref);
    guam.step(ref);
    if gi.progress >= traj.n_trim && tt > traj.t_node(end) + 2
        S = trim_log(S, i);  return;
    end
end
S = trim_log(S, Ng);
end

function S = blank(N)
S = struct('t', zeros(1, N), 'e_pos', zeros(1, N), 'u', zeros(1, N), ...
           'u_ref', zeros(1, N), 'th', zeros(1, N), 'q', zeros(1, N), ...
           'alt', zeros(1, N));
end

function S = record(S, i, t, state, ref)
% body-frame position error, same expression as RSLQR.guidance_error
Rib = RSLQR.rotm_i2b(state(7), state(8), state(9));
e_b = Rib * (ref.pos - state(1:3));
S.t(i)     = t;
S.e_pos(i) = e_b(1);
S.u(i)     = state(4);
S.u_ref(i) = ref.vel(1);
S.th(i)    = state(8);
S.q(i)     = state(11);
S.alt(i)   = -state(3);
end

function S = trim_log(S, n)
fn = fieldnames(S);
for ii = 1:numel(fn), S.(fn{ii}) = S.(fn{ii})(1:n); end
end
