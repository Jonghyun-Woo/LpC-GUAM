% build_descending_model - rebuild Phi and xi_e with the loop settled against
% the DESCENDING reference the trim points were actually computed for.
%
% ---------------------------------------------------------------------------
% WHY
% ---------------------------------------------------------------------------
% The 20 trim points are WH3 points: trimmed at an 11.667 ft/s (700 ft/min)
% descent. The stored model was settled against a LEVEL reference, so its
% equilibrium sits a standing -11.667 ft/s off every anchor in heading-frame w
% and +2 to +4 deg in pitch. The reachable sets still contain that point, but
% most of the room to manoeuvre is spent before the vehicle accelerates at all
% - measured, level flight admits a uniform 1.0 ft/s^2 where sitting on the
% anchors admits 3.5.
%
% Measured against the simulator with the descending reference, the settled
% deviations collapse:
%
%     anchor   level: w err / theta err     descending: w err / theta err
%       1       -11.67 ft/s  /  -0.13 deg      +0.13 ft/s  /  +0.01 deg
%      10       -10.65       /  +2.17          +0.39       /  +0.10
%      17       -11.67       /  +2.19           0.00       /   0.00
%      20       -11.67       /  +4.36          -0.02       /  +0.01
%
% (The cruise anchors need a long hold to show this - the loop carries a mode
% with a 59 s period, and a 30 s settle left anchor 20 apparently 10 ft/s off.
% T_SETTLE below is sized for that, not for the eye.)
%
% ---------------------------------------------------------------------------
% HOW
% ---------------------------------------------------------------------------
% Same construction as compute_lag_coefficient: walk the anchors in order,
% RAMPING the command between them rather than stepping it. A step of 8.44
% ft/s is large enough to drive the loop out of the region where it behaves
% and the equilibrium search diverges - the stored closed_loop_model_gs15.mat
% has ||c|| = 5.7e60 from exactly that. Phi is then measured by nudging each
% of the 34 states and taking ONE step, so the actuator lags, the allocation
% and the aero polynomial are all captured without modelling them by hand.
%
% Output: logger/closed_loop_model_descending.mat
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt       = 0.01;
T_RAMP   = 6.0;     % how gently to move the command between anchors [s]
T_SETTLE = 90.0;    % hold at the anchor this long before measuring    [s]
WDOT     = 11.667;  % WH3, the descent the trim points were built at [ft/s]

params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj.trim_lon(1, :);
trim_lon = traj.trim_lon;
n = numel(u_anchor);

% ---------------------------------------------------------------------------
% TANGENT OR SECANT
% ---------------------------------------------------------------------------
% Phi is a finite difference, so the perturbation size decides WHICH slope it
% measures. The shipped sizes (0.1 ft/s, 0.001 rad) measure the tangent at the
% trim point. But a segment carries u across 8.44 ft/s, and lift, drag and
% control-surface authority all go as u^2 over that span, so the tangent is
% not the slope the transition actually travels on.
%
% Measured against the simulator, the tangent model's error reaches 8 %% of the
% reachable set's half-width at a = 2 and 42 %% at a = 5. That uncertainty has
% to be carried as margin, and it costs about 19 s of the 43 s transition.
%
% 'secant' sizes each perturbation like the excursion actually flown, from
% tests/verify_accel_schedule: u 3.9 ft/s, w 16.3, q 6.9 deg/s, theta 16.9 deg.
% The central difference then measures the average slope across the range the
% vehicle uses rather than the slope at one point of it.
%
% This is a trade, not a free win: a large perturbation mixes nonlinearity INTO
% the slope, which is exactly what ClosedLoopModel.tune_h was written to avoid.
% Whether it helps is decided by measure_model_validity, not by argument.
H_MODE = 'secant';

h0 = zeros(ClosedLoopModel.N_XI, 1);
switch H_MODE
    case 'tangent'
        h0(1:3)   = 0.5;    h0(4:6)   = 0.1;   h0(7:9)   = 1e-3;
        h0(10:12) = 1e-3;   h0(13:21) = 0.5;   h0(22:26) = 1e-3;
        h0(27:32) = 1e-2;   h0(33:34) = 1e-3;
        out_name  = 'closed_loop_model_descending.mat';
    case 'secant'
        h0(1:3)   = 5.0;                       % position error [ft]
        h0(4) = 4.0;  h0(5) = 1.0;  h0(6) = 8.0;      % u, v, w [ft/s]
        h0(7) = 0.02; h0(8) = 0.15; h0(9) = 0.02;     % phi, theta, psi [rad]
        h0(10)= 0.01; h0(11)= 0.06; h0(12)= 0.01;     % p, q, r [rad/s]
        h0(13:21) = 5.0;                       % rotors [rad/s]
        h0(22:26) = 0.02;                      % surfaces [rad]
        h0(27:32) = 0.1;                       % integrators
        h0(33:34) = 0.02;                      % virtual attitude [rad]
        out_name  = 'closed_loop_model_secant.mat';
    otherwise
        error('build_descending_model:mode', 'H_MODE must be tangent or secant.');
end

line = @(c) fprintf('%s\n', repmat(c, 1, 88));

%% build ---------------------------------------------------------------
guam = LpC_GUAM(Config('trim_schedule', params));
M = ClosedLoopModel(dt);
M.w_ref    = WDOT;
M.u_anchor = u_anchor(:)';
M.Phi      = cell(1, n);
M.xi_e     = zeros(ClosedLoopModel.N_XI, n);
M.c        = zeros(ClosedLoopModel.N_XI, n);
M.resid    = zeros(1, n);
M.h        = h0(:);

guam.reset();
xi = M.pack(guam, [0; 0; ClosedLoopModel.Z_REF]);
v_prev = u_anchor(1);

line('='); fprintf('BUILDING AT w_ref = %.3f ft/s (%.0f ft/min descent)\n', ...
                   WDOT, WDOT*60); line('=');
fprintf('  ramp %.1f s between anchors, settle %.0f s, then 34-column Jacobian\n\n', ...
        T_RAMP, T_SETTLE);
tic;
for k = 1:n
    v = u_anchor(k);
    for i = 1:round(T_RAMP/dt)
        vv = v_prev + (v - v_prev) * i / round(T_RAMP/dt);
        xi = M.step_map(guam, xi, vv);
    end
    for i = 1:round(T_SETTLE/dt)
        xi = M.step_map(guam, xi, v);
    end
    M.xi_e(:, k) = xi;
    M.c(:, k)    = M.step_map(guam, xi, v) - xi;
    M.resid(k)   = norm(M.c(:, k));
    M.Phi{k}     = M.jacobian(guam, xi, v, h0);
    v_prev       = v;
    fprintf('  anchor %2d  v = %6.1f  drift %.2e  max|eig| %.5f  (%.0f s)\n', ...
            k, v, M.resid(k), max(abs(eig(M.Phi{k}))), toc);
end
fprintf('\n  built in %.0f s\n\n', toc);

%% health checks -------------------------------------------------------
line('='); fprintf('HEALTH CHECKS\n'); line('=');

lam = cellfun(@(P) max(abs(eig(P))), M.Phi);
fprintf('  max|eig(Phi)| over anchors : %.5f .. %.5f  (must be <= 1)\n', ...
        min(lam), max(lam));
assert(max(lam) < 1 + 1e-6, 'build_descending_model:unstable', ...
    'Some anchor linearised to an unstable loop; the settle did not converge.');

dXE = zeros(size(M.xi_e));
for r = 1:ClosedLoopModel.N_XI
    dXE(r, :) = gradient(M.xi_e(r, :), M.u_anchor);
end
fprintf('  du_e/dv                    : %.4f .. %.4f  (must be ~1)\n', ...
        min(dXE(4,:)), max(dXE(4,:)));

% the point of the whole exercise: does the settled state sit on the anchors?
SEL = [4 6 11 8];
XE  = M.xi_e(SEL, :);
th  = XE(4, :);
XH  = [ XE(1,:).*cos(th) + XE(2,:).*sin(th) ;      % body -> heading
       -XE(1,:).*sin(th) + XE(2,:).*cos(th) ;
        XE(3,:) ; th ];
dev = XH - trim_lon;
fprintf('\n  settled state minus trim point, heading frame:\n');
fprintf('   k   u_trim     du [ft/s]    dw [ft/s]    dq [deg/s]   dth [deg]\n');
for k = 1:n
    fprintf('%4d %8.1f %12.3f %12.3f %13.3f %11.3f\n', k, u_anchor(k), ...
            dev(1,k), dev(2,k), rad2deg(dev(3,k)), rad2deg(dev(4,k)));
end
fprintf('\n  worst |dw| = %.3f ft/s, worst |dtheta| = %.3f deg\n', ...
        max(abs(dev(2,:))), rad2deg(max(abs(dev(4,:)))));

save(fullfile(root,'logger','closed_loop_model_descending.mat'), 'M', 'dXE', ...
     'trim_lon', 'u_anchor', 'WDOT', 'T_RAMP', 'T_SETTLE', '-v7.3');
fprintf('\nsaved logger/closed_loop_model_descending.mat\n');
