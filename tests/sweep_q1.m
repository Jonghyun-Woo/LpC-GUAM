% sweep_q1 - raise the integral weight on velocity tracking and measure what
% it buys.
%
% Qlon0 = [ q1 q2 q3 | q4 q5 q6 ] weights, in order,
%   the accumulated errors in u, w, q, then the deviations of u, w, q from
%   trim. The design solves a Riccati equation on a plant whose input is the
%   acceleration itself, so per channel the gains come out as
%
%       Ki = sqrt(q1/r)              integral action on the tracking error
%       Kx = sqrt(q4/r + 2*Ki)       proportional action on the tracking error
%
% and the control law is a plain PI:  mdes = Kx*error + Ki*integral(error).
% (Kx multiplies Xlon, the deviation from the trim looked up AT THE COMMANDED
% speed, and that trim's speed IS the command, so Xlon(1) = -error.)
%
% q4 is left at 0 throughout: raising it adds proportional action, but the
% mission is a ramp and the failure we are chasing (38 ft position error ->
% +3.8 ft/s bias -> BRT departure) comes from tracking, so q1 is swept alone.
%
% RSLQRConfig.Qlon0 is a Constant property, so it cannot be set at run time.
% ctrl_lon() takes the weights as arguments, so the tables are regenerated
% directly instead.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q1_LIST = [0.01 0.04 0.16 0.64 2.56];      % Ki = 0.1, 0.2, 0.4, 0.8, 1.6
STEP    = 3.0;                              % step size for the response test
K_STEP  = [5 10 15];                        % trims to run the step test at

params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

fprintf('sweeping q1 (and q2, kept equal); q3 = 1000 and q4..q6 = 0 fixed\n');
fprintf('step test %+.1f ft/s at trims %s\n\n', STEP, mat2str(K_STEP));

R = struct([]);
for a = 1:numel(Q1_LIST)
    q1 = Q1_LIST(a);
    Q  = [q1 q1 1000 0 0 0]';

    % --- regenerate the longitudinal gains at this Q -----------------------
    g0  = LpC_GUAM(Config('trim_schedule', params));
    bc  = g0.controller.baseline_controller;
    cfg = RSLQRConfig;
    N = cfg.N_trim;  M = cfg.M_trim;
    Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);  bad = false;
    for jj = 1:M
        for ii = 1:N
            tp = bc.trim_xu_eq(bc.XU0(:,ii,jj));
            [lon, err] = bc.ctrl_lon(tp, Q, cfg.Rlon0, cfg.Wlon0);
            if err, bad = true; end
            Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
        end
    end
    kI = Ki(1,1,8,2);  kX = Kx(1,1,8,2);       % trim 8, level-flight row

    if bad
        fprintf('q1 = %6.3f | UNSTABLE at design time - skipped\n', q1);
        continue;
    end

    % --- step response ----------------------------------------------------
    t63 = nan(1,numel(K_STEP));  os = nan(1,numel(K_STEP));
    for b = 1:numel(K_STEP)
        [t63(b), os(b)] = step_test(params, traj.trim_lon(1,K_STEP(b)), ...
                                    STEP, dt, Ki, Kx);
    end

    % --- ungoverned mission ----------------------------------------------
    [eu, ep, viol, ad, blown] = mission_test(params, brtV, dt, Ki, Kx);

    R(end+1).q1 = q1;      %#ok<SAGROW>
    R(end).Ki = kI;   R(end).Kx = kX;
    R(end).t63 = mean(t63, 'omitnan');  R(end).os = mean(os, 'omitnan');
    R(end).eu = eu;   R(end).ep = ep;   R(end).viol = viol;
    R(end).ad = ad;   R(end).blown = blown;

    fprintf('q1 = %6.3f done  (Ki %.3f, Kx %.3f)\n', q1, kI, kX);
end

%% Report
fprintf('\n%7s | %6s %6s | %7s %8s | %8s %8s %8s %8s\n', ...
        'q1','Ki','Kx','t63[s]','over[%]','|u-r|','e_pos','BRT>=0','alt dev');
fprintf('%s\n', repmat('-', 1, 84));
for a = 1:numel(R)
    if R(a).blown
        fprintf('%7.3f | %6.3f %6.3f | %7s %8s | %8s %8s %8s %8s  DIVERGED\n', ...
                R(a).q1, R(a).Ki, R(a).Kx, '--','--','--','--','--','--');
    else
        fprintf('%7.3f | %6.3f %6.3f | %7.2f %8.1f | %8.2f %8.1f %8.1f %8.1f\n', ...
                R(a).q1, R(a).Ki, R(a).Kx, R(a).t63, 100*R(a).os, ...
                R(a).eu, R(a).ep, R(a).viol, R(a).ad);
    end
end
fprintf(['\nbaseline is q1 = 0.010 (the shipped design).\n' ...
         't63 = time to cover 63%% of the commanded step, averaged over the\n' ...
         'test trims; over = peak / commanded step - 1.\n']);
save(fullfile(root,'logger','sweep_q1.mat'), 'R', 'Q1_LIST');

% =========================================================================
function [t63, os] = step_test(params, v, dv, dt, Ki, Kx)
guam = LpC_GUAM(Config('trim_schedule', params));
set_gains(guam, Ki, Kx);
guam.reset();
mk = @(p,vv) struct('pos',[p;0;-80],'vel',[vv;0;0],'chi',0,'chi_dot',0);
p = 0;
for i = 1:round(20/dt), p = p + v*dt;  guam.step(mk(p,v)); end
u0 = guam.state(4);
v2 = v + dv;  N = round(12/dt);  U = zeros(1,N);
for i = 1:N
    p = p + v2*dt;  guam.step(mk(p,v2));  U(i) = guam.state(4);
end
t = (1:N)*dt;
lvl = u0 + 0.632*dv;                    % against the COMMANDED step
i63 = find(U >= lvl, 1, 'first');
if isempty(i63), t63 = NaN; else, t63 = t(i63); end
os = max(U - u0)/dv - 1;
end

function [eu, ep, viol, ad, blown] = mission_test(params, brtV, dt, Ki, Kx)
guam = LpC_GUAM(Config('trim_schedule', params));
set_gains(guam, Ki, Kx);
guam.reset();
rt = guam.refTraj;  N = size(rt.pos,2);
eu = 0; ep = 0; nv = 0; ad = 0; blown = false;
for i = 1:N
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, blown = true; break; end
    R = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e = R * (rt.pos(:,i) - s(1:3));
    eu = max(eu, abs(s(4) - rt.vel(1,i)));
    ep = max(ep, abs(e(1)));
    ad = max(ad, abs(-s(3) - 80));
    nv = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
viol = 100*nv/N;
end

function set_gains(guam, Ki, Kx)
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;
bc.LON.Kx = Kx;
end
