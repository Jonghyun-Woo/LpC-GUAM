% verify_lag_coefficient - does the analytic gap equation predict the real
% speed error of the flown mission?
%
% compute_lag_coefficient.m derives, without simulating any trajectory,
%
%       d(i+1) = Phi(v) d(i) - xi_e'(v) * s(i) * dt                      (*)
%
% where d is the gap between the closed-loop state and the trim the reference
% is currently pointing at, and s is the rate at which the commanded speed is
% rising. This script checks (*) against the thing it claims to predict.
%
% THE CHECK
%   Take the ACTUAL reference profile of the trim-schedule mission, read its
%   slope s(t) off it, and propagate (*) along that profile. That produces a
%   predicted forward-speed error at every instant. Then fly the real
%   nonlinear mission and measure the forward-speed error at every instant.
%   Put the two curves side by side.
%
%   This tests the whole equation, not a summary number: if (*) is right the
%   curves agree in shape, in magnitude, and in where the peaks fall.
%
% WHAT AGREEMENT WOULD MEAN
%   That e can be predicted from Phi and the trim table alone, so the
%   trim-point placement rule can be solved for the step length instead of
%   being searched by flying candidate schedules.
%
% WHAT THE TWO PIECES OF e ARE
%   e0  the error that survives when the command is HELD (s = 0). It shows up
%       here as xi_e(u) - v: where the loop actually settles versus what was
%       commanded. This is the piece that decides whether a speed band is
%       crossable at all, because it does not shrink when the step shrinks.
%   c   the ramp-proportional piece, e = c * s.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt    = 0.01;
T_SEG = 2.0;
U_ROW = 4;

L = load(fullfile(root, 'logger', 'lag_coefficient.mat'));
S = load(fullfile(root, 'logger', 'closed_loop_model_ramped.mat'));  M = S.M;

params = struct('filter_mode', 'off');  params.T_seg = T_SEG;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));

%% ---------------------------------------------------------------------
%  1. e0 : where does the loop settle when the command is simply HELD?
%% ---------------------------------------------------------------------
off = M.xi_e(U_ROW, :) - M.u_anchor;
fprintf('--- e0: settled speed minus commanded speed, command HELD ---\n');
fprintf('  k |    u    | settled u | offset [ft/s]\n');
for k = 1:numel(M.u_anchor)
    fprintf('%3d | %7.2f | %9.2f | %+8.3f\n', k, M.u_anchor(k), ...
            M.xi_e(U_ROW, k), off(k));
end
fprintf('offset: min %+.3f  median %+.3f  max %+.3f ft/s\n\n', ...
        min(off), median(off), max(off));

%% ---------------------------------------------------------------------
%  2. fly the real mission
%% ---------------------------------------------------------------------
guam = LpC_GUAM(Config('trim_schedule', params));
rt   = guam.refTraj;  N = size(rt.pos, 2);
guam.reset();
u_act = zeros(1, N);
for k = 1:N
    ref = struct('pos', rt.pos(:,k), 'vel', rt.vel(:,k), ...
                 'chi', rt.chi(k),   'chi_dot', rt.chidot(k));
    u_act(k) = guam.state(4);
    guam.step(ref);
end
u_ref = rt.vel(1, :);
err_meas = u_act - u_ref;

%% ---------------------------------------------------------------------
%  3. propagate (*) along the SAME reference profile
%     Phi and xi_e' are interpolated in the commanded speed, exactly the way
%     ClosedLoopModel.at() interpolates Phi for the linear rollout.
%% ---------------------------------------------------------------------
dXE = L.dXE;             % 34 x n, dxi_e/dv, position rows already zeroed
                         % by compute_lag_coefficient (see the note there)
s_t = [diff(u_ref) / dt, 0];                   % commanded slope [ft/s^2]

d = zeros(ClosedLoopModel.N_XI, 1);
err_pred = zeros(1, N);
for i = 1:N
    err_pred(i) = d(U_ROW);
    [P, ~, ~] = M.at(u_ref(i));
    ve = interp_cols(dXE, M.u_anchor, u_ref(i));
    d  = P*d - ve*s_t(i)*dt;
end

%% ---------------------------------------------------------------------
%  4. compare, per segment
%% ---------------------------------------------------------------------
fprintf('--- measured vs predicted forward-speed error, per segment ---\n');
fprintf('slope during every transition segment = %.2f ft/s^2\n\n', 8.44/T_SEG);
fprintf('seg |  u range [ft/s]  | measured peak | predicted peak | pred-meas\n');
fprintf('%s\n', repmat('-', 1, 72));
tn = traj.t_node;
pk_m = zeros(1, numel(tn)-1);  pk_p = pk_m;
for k = 1:numel(tn)-1
    i0 = max(round(tn(k)/dt), 1);  i1 = min(round(tn(k+1)/dt), N);
    pk_m(k) = max(abs(err_meas(i0:i1)));
    pk_p(k) = max(abs(err_pred(i0:i1)));
    fprintf('%3d | %6.1f -> %6.1f | %13.2f | %14.2f | %+9.2f\n', ...
            k, traj.trim_lon(1,k), traj.trim_lon(1,k+1), pk_m(k), pk_p(k), ...
            pk_p(k)-pk_m(k));
end
fprintf('%s\n', repmat('-', 1, 72));
fprintf('mission peak: measured %.2f  predicted %.2f ft/s\n', ...
        max(abs(err_meas)), max(abs(err_pred)));
rel = abs(pk_p - pk_m) ./ max(pk_m, eps);
fprintf('per-segment relative error: median %.1f %%, worst %.1f %%\n', ...
        100*median(rel), 100*max(rel));

%% ---------------------------------------------------------------------
%  5. figure
%% ---------------------------------------------------------------------
f = figure('Position', [100 100 1000 620]);
subplot(2,1,1);
plot(rt.time, u_ref, 'k-', 'LineWidth', 1.2); hold on;
plot(rt.time, u_act, 'b-', 'LineWidth', 1.2);
xline(tn, ':', 'Color', [.7 .7 .7]);
ylabel('forward speed [ft/s]'); legend('reference','flown','Location','northwest');
title('trim-schedule mission, T_{seg} = 2 s, no safety filter'); grid on;

subplot(2,1,2);
plot(rt.time, err_meas, 'b-', 'LineWidth', 1.2); hold on;
plot(rt.time, err_pred, 'r--', 'LineWidth', 1.4);
xline(tn, ':', 'Color', [.7 .7 .7]); yline(0, 'k-');
xlabel('time [s]'); ylabel('speed error: flown - reference [ft/s]');
legend('measured (nonlinear)', 'predicted by the gap equation', ...
       'Location','southwest');
grid on;
saveas(f, fullfile(root, 'logger', 'lag_coefficient_verify.png'));
fprintf('\nsaved logger/lag_coefficient_verify.png\n');

save(fullfile(root,'logger','lag_verify.mat'), ...
     'err_meas','err_pred','u_ref','u_act','pk_m','pk_p','off','tn');

% -------------------------------------------------------------------------
function col = interp_cols(A, x, xq)
% linear interpolation of the columns of A over the breakpoints x
if xq <= x(1),   col = A(:,1);    return;  end
if xq >= x(end), col = A(:,end);  return;  end
k = find(x <= xq, 1, 'last');  k = min(k, numel(x)-1);
a = (xq - x(k)) / (x(k+1) - x(k));
col = (1-a)*A(:,k) + a*A(:,k+1);
end
