% build_closed_loop_model - measure Phi(v) and xi_e(v) at the 20 trim points
% and check that the result is usable.
%
% Produces logger/closed_loop_model.mat, which is what the linear rollout
% will read instead of stepping the plant.
%
% Three checks are run afterwards, in the order that isolates faults:
%   (1) did each trim actually settle?      -> residual of xi_e
%   (2) is Phi stable (Schur)?              -> holding a command must settle
%   (3) does the linear model track the nonlinear one from a small nudge?
%       If the two disagree from the FIRST step the measurement is wrong;
%       disagreeing later is just nonlinearity and is expected.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
guam = LpC_GUAM(Config('trim_schedule', params));

%% Perturbation sizes, by physical channel
h0 = zeros(ClosedLoopModel.N_XI, 1);
h0(1:3)   = 0.5;        % e_pos      [ft]
h0(4:6)   = 0.1;        % u, v, w    [ft/s]
h0(7:9)   = 1e-3;       % euler      [rad]
h0(10:12) = 1e-3;       % rates      [rad/s]
h0(13:21) = 0.5;        % propellers [rad/s]
h0(22:26) = 1e-3;       % surfaces   [rad]
h0(27:32) = 1e-2;       % integrators
h0(33:34) = 1e-3;       % phi_v, theta_v

%% Build
M = ClosedLoopModel(dt);
fprintf('building closed-loop model at %d trim points\n\n', traj.n_trim);
tic;
M.build(guam, traj.trim_lon(1, :), h0, 15.0);
fprintf('\nbuild time %.0f s\n', toc);

%% (1) how far xi_e is from a true fixed point, and which state drifts
fprintf('\n--- (1) one-step drift c (absorbed by the affine form) ---\n');
fprintf('max |c| %.2e, median %.2e\n', max(M.resid), median(M.resid));
[~, kw] = max(M.resid);
[~, jw] = max(abs(M.c(:, kw)));
nmxi = [repmat({'e_pos'},1,3) {'u','v','w','phi','theta','psi','p','q','r'} ...
        repmat({'eng'},1,9) repmat({'srf'},1,5) ...
        repmat({'lon_i'},1,3) repmat({'lat_i'},1,3) {'phi_v','theta_v'}];
fprintf('worst trim %d: largest drifting state is #%d (%s) at %.2e per step\n', ...
        kw, jw, nmxi{jw}, M.c(jw, kw));

%% (2) stability
fprintf('\n--- (2) is each Phi Schur? ---\n');
lam = zeros(1, traj.n_trim);
for k = 1:traj.n_trim, lam(k) = max(abs(eig(M.Phi{k}))); end
fprintf('max|eig| over all trims: %.6f  (must be < 1)\n', max(lam));
if max(lam) >= 1
    fprintf('WARNING trims %s are not Schur\n', mat2str(find(lam >= 1)));
end
% slowest mode -> settling time -> the horizon T the governor should use
Tset = -log(0.02) ./ -log(max(lam, eps)) * dt;
fprintf('\n k |   u    | max|eig| | implied settling T [s]\n');
for k = 1:traj.n_trim
    fprintf('%2d | %6.1f | %8.5f | %6.2f\n', k, M.u_anchor(k), lam(k), ...
            -log(0.02)/(-log(lam(k)))*dt);
end
fprintf('\nlongest implied T = %.2f s  (current governor uses 6.0 s)\n', max(Tset));

%% (3) linear vs nonlinear from a small nudge
fprintf('\n--- (3) small-perturbation agreement ---\n');
K_TEST = [3 8 14];  N_CHK = 200;                 % 2 s
fprintf('%4s %8s | %10s %10s %10s\n', 'trim', 'u', 'du@1step', 'du@0.5s', 'du@2s');
for kk = K_TEST
    v  = M.u_anchor(kk);
    xe = M.xi_e(:, kk);
    d0 = zeros(ClosedLoopModel.N_XI, 1);
    d0(4) = 1.0;                                  % 1 ft/s faster than trim
    % nonlinear
    xi = xe + d0;  XN = zeros(1, N_CHK+1);  XN(1) = xi(4);
    for i = 1:N_CHK, xi = M.step_map(guam, xi, v);  XN(i+1) = xi(4); end
    % linear, affine form
    d = d0;  XL = zeros(1, N_CHK+1);  XL(1) = xe(4) + d(4);
    for i = 1:N_CHK
        d = M.c(:,kk) + M.Phi{kk}*d;  XL(i+1) = xe(4) + d(4);
    end
    e = abs(XL - XN);
    fprintf('%4d %8.1f | %10.2e %10.2e %10.2e\n', kk, v, e(2), e(51), e(end));
end

%% Save
out = fullfile(root, 'logger', 'closed_loop_model.mat');
save(out, 'M', 'h0', '-v7.3');
fprintf('\nsaved %s\n', out);
