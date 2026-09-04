% measure_tau_step - identify the velocity-loop time constant directly, by
% step response, instead of inferring it from the ramp-following bias.
%
% The bias argument gives tau = k_pos*e_pos/a = 0.91 s, but that rests on a
% first-order assumption and on the loop reaching steady state inside a 2 s
% segment (it is only ~2.2 tau, so it does not). A step response assumes
% nothing: settle at a trim, step the command, and read off when the speed
% has covered 63.2 % of the step.
%
% Two variants are run at each trim:
%   k_pos = 0    inner velocity loop alone - this is the tau in the formula
%   k_pos = 0.1  as flown, inner + outer position loop
% Their difference shows how much the outer loop is doing.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  DV = 3.0;                 % step size [ft/s]
T_SETTLE = 15;  T_STEP = 12;
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
U_A = traj.trim_lon(1, :);
a_plan = (U_A(2) - U_A(1)) / params.T_seg;

fprintf('step size %.1f ft/s, plan acceleration a = %.2f ft/s^2\n\n', DV, a_plan);
fprintf(['  k |   u    | tau_inner | tau_flown | rise10-90 | overshoot | ' ...
         'settle 2%% | implied bias a*tau\n']);
fprintf('%s\n', repmat('-', 1, 96));

R = struct('tau_in', zeros(1,numel(U_A)), 'tau_fl', zeros(1,numel(U_A)));
for k = 1:numel(U_A)
    v = U_A(k);
    [ti, ~,  ~,  ~ ] = step_test(params, v, DV, 0.0, dt, T_SETTLE, T_STEP);
    [tf, tr, os, ts] = step_test(params, v, DV, 0.1, dt, T_SETTLE, T_STEP);
    R.tau_in(k) = ti;  R.tau_fl(k) = tf;
    fprintf('%3d | %6.1f | %9.3f | %9.3f | %9.3f | %8.1f%% | %8.3f | %8.2f ft/s\n', ...
            k, v, ti, tf, tr, 100*os, ts, a_plan*ti);
end

fprintf('\ninner loop  tau: min %.3f  median %.3f  max %.3f s\n', ...
        min(R.tau_in), median(R.tau_in), max(R.tau_in));
fprintf('as flown    tau: min %.3f  median %.3f  max %.3f s\n', ...
        min(R.tau_fl), median(R.tau_fl), max(R.tau_fl));
fprintf(['\nbias inferred from the ramp argument was 3.83 ft/s (tau = 0.91 s);\n' ...
         'compare against a*tau_inner above, trim by trim.\n']);

% =========================================================================
function [tau, t_rise, overshoot, t_settle] = step_test(params, v, dv, kpos, ...
                                                        dt, T_SETTLE, T_STEP)
guam = LpC_GUAM(Config('trim_schedule', params));
guam.controller.k_pos = kpos;
guam.reset();
mkref = @(p, vv) struct('pos',[p;0;-80],'vel',[vv;0;0],'chi',0,'chi_dot',0);

p = 0;
for i = 1:round(T_SETTLE/dt)
    p = p + v*dt;  guam.step(mkref(p, v));
end
u0 = guam.state(4);

v2 = v + dv;  N = round(T_STEP/dt);  U = zeros(1, N);
for i = 1:N
    p = p + v2*dt;  guam.step(mkref(p, v2));
    U(i) = guam.state(4);
end
t = (1:N)*dt;
uf = mean(U(end-round(0.5/dt):end));        % achieved final value
d  = uf - u0;
if abs(d) < 1e-6
    tau = NaN; t_rise = NaN; overshoot = NaN; t_settle = NaN;  return;
end

cross = @(f) first_cross(t, U, u0 + f*d);
tau       = cross(0.632);
t_rise    = cross(0.90) - cross(0.10);
overshoot = (max((U - u0)/d) - 1);
out       = find(abs(U - uf) > 0.02*abs(d), 1, 'last');
if isempty(out), t_settle = 0; else, t_settle = t(min(out+1, N)); end
end

function tc = first_cross(t, y, lvl)
i = find((y - lvl) .* sign(lvl - y(1)) >= 0, 1, 'first');
if isempty(i) || i == 1, tc = NaN; return; end
% linear interpolation between samples
tc = t(i-1) + (lvl - y(i-1)) / (y(i) - y(i-1)) * (t(i) - t(i-1));
end
