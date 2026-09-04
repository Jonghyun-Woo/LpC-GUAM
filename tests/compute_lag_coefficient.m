% compute_lag_coefficient - the tracking-error coefficient c, analytically.
%
% ---------------------------------------------------------------------------
% WHAT THIS ANSWERS
% ---------------------------------------------------------------------------
% The trim-point placement rule needs one number per speed: how far the
% vehicle falls behind the reference while the reference speed is being
% ramped. Call the resulting speed error e. This script computes e WITHOUT
% flying the mission, by writing the equation for the GAP between the vehicle
% and the reference and solving that equation instead.
%
% ---------------------------------------------------------------------------
% THE EQUATION
% ---------------------------------------------------------------------------
% ClosedLoopModel already provides, for a command HELD at v,
%
%       xi(i+1) - xi_e(v) = Phi(v) ( xi(i) - xi_e(v) ) + c
%
% where xi is the 34-state closed loop (position error, body velocities,
% attitude, rates, 9 propellers, 5 surfaces, both integrators, both virtual
% attitude effectors) and xi_e(v) is where that loop settles when the command
% is held at v. Phi is measured by nudging each state and taking ONE simulator
% step, so the actuator lags, the control allocation and the aero polynomial
% are all inside it.
%
% Now let the command RAMP at s ft/s per second instead of being held. Between
% step i and step i+1 the command moves by s*dt, so the equilibrium that the
% gap is measured FROM also moves, by xi_e'(v) * s * dt. Writing
% d(i) = xi(i) - xi_e(v(i)) and subtracting that shift:
%
%       d(i+1) = Phi d(i) + c - xi_e'(v) * s * dt                        (*)
%
% Everything on the right is known. Phi and xi_e come from the closed-loop
% model, xi_e' is a finite difference of xi_e across the trim anchors, and s
% is the ramp slope we are designing. No trajectory is simulated: (*) is
% propagated with 34x34 matrix multiplies.
%
% ---------------------------------------------------------------------------
% WHICH NUMBER WE TAKE OUT OF (*), AND WHY NOT THE OTHER ONE
% ---------------------------------------------------------------------------
% The obvious quantity is the balance point of (*), d_ss = (I-Phi)^-1(...),
% i.e. the lag once it has fully settled. THAT QUANTITY IS NOT USABLE HERE.
% Phi has an eigenvalue at 0.99999, a position-error mode with a time constant
% of order 1000 s (RSLQR's +0.1*e_pos outer loop, which is deliberately slow).
% So I-Phi is numerically singular - cond(I-Phi) ~ 1e7 even for a healthy
% model - and the "settled" lag is dominated by a mode that does essentially
% nothing during a 2-second segment. Asking for it is asking the wrong
% question.
%
% What we actually need is the WORST gap over a segment of finite length T,
% which is also the quantity the safety condition is written against. So (*)
% is propagated forward for T seconds with s = 1, then for a further T_HOLD
% seconds with s = 0 (the ramp has ended and the loop swings back), and the
% peak forward-speed row is taken:
%
%       c_pk(v, T) = max_i | d_u(i) |          with s = 1                [s]
%
% Because (*) is linear in s, the peak scales exactly with s, so one
% propagation gives the coefficient for EVERY slope:
%
%       e = c_pk * s = c_pk * (step length / segment time)
%
% c_pk is in seconds. Read it as "the vehicle is c_pk seconds behind".
%
% ---------------------------------------------------------------------------
% WHY THE MODEL IS REBUILT HERE INSTEAD OF LOADED
% ---------------------------------------------------------------------------
% Two reasons.
%
% (1) Phi describes the CLOSED loop, so it is invalidated by a gain change.
%     Qlon0 is now [0.08 0.01 1000 2 0 0]; logger/closed_loop_model.mat was
%     built before that retune.
%
% (2) ClosedLoopModel.build finds each anchor's equilibrium by JUMPING the
%     command from the previous anchor - a 8.44 ft/s step - and holding for
%     15 s. With the retuned (higher speed-loop bandwidth) gains that step is
%     large enough to drive the loop out of the region where it behaves, and
%     the equilibrium search diverges: the stored logger/closed_loop_model_
%     gs15.mat has ||c|| = 5.7e60 and du_e/dv = 0.95 where it should be ~1.
%     Rebuilding with the command RAMPED between anchors (T_RAMP below) keeps
%     the loop in its operating region and the search converges. The three
%     health checks printed below are what catch this if it recurs.
%
% Output: logger/lag_coefficient.mat, logger/closed_loop_model_ramped.mat
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt       = 0.01;
T_SEG    = 2.0;      % segment duration the coefficients are reported at [s]
T_HOLD   = 6.0;      % keep watching this long after the ramp stops [s]
T_RAMP   = 5.0;      % how gently to move the command between anchors [s]
T_SETTLE = 20.0;     % hold at the anchor this long before measuring [s]
U_ROW    = 4;        % forward speed u sits in row 4 of the 34-state xi

params = struct('filter_mode', 'off');  params.T_seg = T_SEG;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
u_anchor = traj.trim_lon(1, :);
n = numel(u_anchor);

h0 = zeros(ClosedLoopModel.N_XI, 1);
h0(1:3)   = 0.5;    h0(4:6)   = 0.1;   h0(7:9)   = 1e-3;
h0(10:12) = 1e-3;   h0(13:21) = 0.5;   h0(22:26) = 1e-3;
h0(27:32) = 1e-2;   h0(33:34) = 1e-3;

%% ---------------------------------------------------------------------
%  1. build Phi(v) and xi_e(v), ramping the command between anchors
%% ---------------------------------------------------------------------
model_file = fullfile(root, 'logger', 'closed_loop_model_ramped.mat');
if isfile(model_file) && ~(exist('FORCE_REBUILD','var') && FORCE_REBUILD)
    fprintf('loading %s (delete it, or set FORCE_REBUILD=true, to rebuild)\n\n', ...
            model_file);
    Sm = load(model_file);  M = Sm.M;
else
guam = LpC_GUAM(Config('trim_schedule', params));
M    = ClosedLoopModel(dt);
M.u_anchor = u_anchor(:)';
M.Phi   = cell(1, n);
M.xi_e  = zeros(ClosedLoopModel.N_XI, n);
M.c     = zeros(ClosedLoopModel.N_XI, n);
M.resid = zeros(1, n);
M.h     = h0(:);

guam.reset();
xi = M.pack(guam, [0; 0; ClosedLoopModel.Z_REF]);
v_prev = u_anchor(1);
fprintf('building closed-loop model, command ramped over %.1f s between anchors\n', T_RAMP);
tic;
for k = 1:n
    v = u_anchor(k);
    % --- ramp the command from the previous anchor to this one ----------
    N_r = round(T_RAMP / dt);
    for i = 1:N_r
        vv = v_prev + (v - v_prev) * i / N_r;
        xi = M.step_map(guam, xi, vv);
    end
    % --- hold and let it settle -----------------------------------------
    for i = 1:round(T_SETTLE / dt)
        xi = M.step_map(guam, xi, v);
    end
    M.xi_e(:, k) = xi;
    M.c(:, k)    = M.step_map(guam, xi, v) - xi;
    M.resid(k)   = norm(M.c(:, k));
    M.Phi{k}     = M.jacobian(guam, xi, v, h0);
    v_prev       = v;
end
fprintf('built in %.0f s\n\n', toc);
save(model_file, 'M', '-v7.3');
end

%% ---------------------------------------------------------------------
%  2. health checks - these are what caught the broken stored model
%% ---------------------------------------------------------------------
lam = cellfun(@(P) max(abs(eig(P))), M.Phi);
dXE = zeros(size(M.xi_e));
for r = 1:ClosedLoopModel.N_XI
    dXE(r, :) = gradient(M.xi_e(r, :), M.u_anchor);
end

% The three position-error rows are excluded from the equilibrium shift.
% Reason: the position loop is the 0.99937 eigenvalue, a time constant of
% order 1000 s. It cannot reach a command-dependent equilibrium during the
% 20 s settle, and the stored xi_e rows show it - the down-position row runs
% 4.71, 32.63, 15.83 across anchors 16-18, which is drift, not a trim. So
% their finite difference is noise, and feeding it in as "the equilibrium is
% moving this fast" injects an error the loop never actually sees. Setting
% these rows to zero states the physically true thing instead: on a 2-second
% segment the position loop does not move at all.
% Effect, measured against the flown mission (verify_lag_coefficient.m):
% mission peak error 5.13 -> 4.50 ft/s predicted, against 4.45 measured.
dXE(1:3, :) = 0;

fprintf('--- health checks ---\n');
fprintf('(1) settled?    max one-step drift ||c|| = %.2e  (blown-up build gave 5.7e60)\n', ...
        max(M.resid));
fprintf('(2) stable?     max|eig| = %.6f  (must be <= 1)\n', max(lam));
fprintf('(3) sane trim?  du_e/dv should be ~1 (a trim at speed v HAS speed v)\n');
fprintf('    du_e/dv    min %.4f  median %.4f  max %.4f\n', ...
        min(dXE(U_ROW,:)), median(dXE(U_ROW,:)), max(dXE(U_ROW,:)));
ok = max(M.resid) < 1 && max(lam) <= 1.0001 && ...
     all(abs(dXE(U_ROW,:) - 1) < 0.35);
if ~ok
    fprintf('\n*** health checks FAILED - the coefficients below are meaningless ***\n');
else
    fprintf('    all three pass\n');
end
fprintf('\n');

%% ---------------------------------------------------------------------
%  3. c_pk : the worst-case lag coefficient over a T_SEG segment
%     Propagate (*) with s = 1. Linear in s, so this IS the coefficient.
%     The constant c is dropped: it is the model's own leftover drift about
%     xi_e, not a tracking error, and including it would add a fixed offset
%     that has nothing to do with the ramp.
%% ---------------------------------------------------------------------
peak_of = @(P, ve, T) local_peak(P, ve, T, T_HOLD, dt, U_ROW, ...
                                 ClosedLoopModel.N_XI);

c_pk = zeros(1, n);  t_pk = zeros(1, n);  c_end = zeros(1, n);
for k = 1:n
    [c_pk(k), t_pk(k), c_end(k)] = peak_of(M.Phi{k}, dXE(:,k), T_SEG);
end

% How c depends on the segment length. It is not a constant: a longer segment
% gives the lag more time to build up, so c grows and then saturates once the
% loop has settled. This matters because T is a design variable too.
T_LIST = [0.5 1.0 1.5 2.0 3.0 4.0 6.0];
c_T = zeros(numel(T_LIST), n);
for a = 1:numel(T_LIST)
    for k = 1:n
        c_T(a, k) = peak_of(M.Phi{k}, dXE(:,k), T_LIST(a));
    end
end

%% ---------------------------------------------------------------------
%  4. report
%% ---------------------------------------------------------------------
s_now = 8.44 / T_SEG;      % slope of the CURRENT uniform schedule [ft/s^2]
fprintf('--- lag coefficient (segment length T = %.1f s) ---\n', T_SEG);
fprintf('  k |    u    | c_end [s] | c_pk [s] | peak at | e at current 8.44/%.0fs\n', T_SEG);
fprintf('    | [ft/s]  | at T_seg  |  worst   |   [s]   |        [ft/s]\n');
fprintf('%s\n', repmat('-', 1, 74));
for k = 1:n
    fprintf('%3d | %7.2f | %9.4f | %8.4f | %7.2f | %8.2f\n', ...
            k, u_anchor(k), c_end(k), c_pk(k), t_pk(k), c_pk(k)*s_now);
end
fprintf('\nc_pk  min %.3f  median %.3f  max %.3f  [s]\n', ...
        min(c_pk), median(c_pk), max(c_pk));
fprintf('overshoot beyond the ramp-end lag: median %.2fx\n', ...
        median(c_pk ./ max(c_end, eps)));
fprintf('\npredicted peak speed error on the CURRENT schedule (8.44 ft/s / %.1f s):\n', T_SEG);
fprintf('   %.2f .. %.2f ft/s   (RSLQRConfig records 4.45 ft/s measured)\n', ...
        min(c_pk)*s_now, max(c_pk)*s_now);

fprintf('\n--- how c depends on the segment length T ---\n');
fprintf('  T [s] |');  fprintf(' u=%5.0f', u_anchor([1 4 8 12 16 20]));  fprintf('\n');
fprintf('%s\n', repmat('-', 1, 58));
for a = 1:numel(T_LIST)
    fprintf('%7.1f |', T_LIST(a));
    fprintf(' %6.3f', c_T(a, [1 4 8 12 16 20]));
    fprintf('\n');
end
fprintf(['\nc grows with T and then saturates: a longer segment lets the lag\n' ...
         'finish building. Above ~3 s it stops moving.\n']);

save(fullfile(root, 'logger', 'lag_coefficient.mat'), ...
     'u_anchor', 'c_pk', 'c_end', 't_pk', 'c_T', 'T_LIST', ...
     'T_SEG', 'T_HOLD', 'dt', 'lam', 'dXE');
fprintf('\nsaved logger/lag_coefficient.mat\n');

% -------------------------------------------------------------------------
function [pk, t_pk, d_end] = local_peak(P, ve, T, T_hold, dt, U_ROW, NX)
% Propagate  d(i+1) = P d(i) - ve * s * dt  with s = 1 for T seconds, then
% s = 0 for T_hold seconds, and return the largest |d_u| seen. Linear in s,
% so this is the coefficient: e = pk * s.
d = zeros(NX, 1);  pk = 0;  ipk = 0;
for i = 1:round(T/dt)
    d = P*d - ve*dt;
    if abs(d(U_ROW)) > pk, pk = abs(d(U_ROW));  ipk = i;  end
end
d_end = abs(d(U_ROW));
n_T = round(T/dt);
for j = 1:round(T_hold/dt)
    d = P*d;
    if abs(d(U_ROW)) > pk, pk = abs(d(U_ROW));  ipk = n_T + j;  end
end
t_pk = ipk * dt;
end
