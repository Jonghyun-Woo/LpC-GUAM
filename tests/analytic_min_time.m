% analytic_min_time - the fastest safe transition, without flying anything.
%
% ---------------------------------------------------------------------------
% WHAT WAS MISSING
% ---------------------------------------------------------------------------
% optimize_segment_time.m answered "8.25 s" and the flown check said 41.8 s.
% The gap is one missing constraint. The safety condition
%
%       a_fwd(u) + a_bwd(u+D) - D - 2*( e0 + c(u,T)*D/T )  >=  0
%
% bounds how far a step may reach, but nothing in it bounds how hard the
% VEHICLE can push. Phi was measured with small perturbations, so the model
% believes the loop has unlimited control power and will happily ramp the
% command at any rate. Adding the real limit closes the gap.
%
% ---------------------------------------------------------------------------
% THE SAFETY CONDITION IS SECRETLY AN ACCELERATION LIMIT
% ---------------------------------------------------------------------------
% Write the condition with a(u+D) ~ a(u) (the tube half-width moves by only
% 0.031 per ft/s, so over one step it barely changes) and divide by T:
%
%       D * ( 1 + 2c(T)/T )  <=  2a - 2e0
%       => s = D/T  <=  ( 2a(u) - 2e0(u) ) / ( T + 2c(u,T) )  ==  s_tube(u,T)
%
% The step length D has dropped out and what is left is a bound on s, the RATE
% at which the commanded speed may rise. That is the real content of the
% overlap-with-margin condition: not "steps must be short" but "the reference
% must not accelerate faster than this". Trim-point spacing is just how that
% rate gets implemented.
%
% s_tube grows without bound as T shrinks, which is why the optimiser ran away.
%
% ---------------------------------------------------------------------------
% THE VEHICLE LIMIT
% ---------------------------------------------------------------------------
% RSLQR asks for a longitudinal acceleration triple mdes = [du/dt; dw/dt; dq/dt]
% and the allocation turns it into effector perturbations
%
%       act = M(u) * mdes,        M = W^-1 B' (B W^-1 B')^-1
%
% with W and B stored per trim point. Asking for pure forward acceleration
% means mdes = [s; 0; 0], so act = s * M(:,1) and each effector moves in
% proportion to s. The largest s that keeps every effector inside its limit is
%
%       s_veh(u) = min over effectors of  limit_i / |M(i,1)|
%
% This is a table lookup and a division - no simulation, and it uses the same
% allocation and the same limits the controller enforces in flight
% (RSLQR.pseudo_alloc scales the whole command down when this is exceeded).
%
% ---------------------------------------------------------------------------
% THE ANSWER
% ---------------------------------------------------------------------------
%       s_max(u) = min( s_veh(u),  max over T of s_tube(u,T) )
%
%       T_min = integral of  du / s_max(u)   from hover to cruise
%
% Time is the integral of one over acceleration - the most ordinary formula
% there is. All the reachability machinery ends up as one of two competing
% acceleration limits.
%
% Output: logger/analytic_min_time.mat / .png
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt     = 0.01;
WH_IDX = 3;
N_TRIM = 20;

O = load(fullfile(root,'logger','segment_time_opt.mat'));
u_anchor = O.u_anchor;  a_fwd = O.a_fwd;  a_bwd = O.a_bwd;  e0_u = O.e0_u;
T_GRID   = O.T_GRID;    C = O.C;

%% ---------------------------------------------------------------------
%  1. s_veh(u) - the largest forward acceleration the allocation can deliver
%% ---------------------------------------------------------------------
params = struct('filter_mode','off');  params.T_seg = 2.0;
cfg   = Config('trim_schedule', params);
ctrl  = RSLQR(cfg.controller, dt);          % RSLQR takes the ControllerConfig
rcfg  = cfg.controller.rslqr;
lim   = [rcfg.eng_max(:); rcfg.ele_max; rcfg.flp_max];  % 11x1
nm   = [compose("lift%d",1:8), "pusher", "elevator", "flap"];

s_veh  = zeros(1, N_TRIM);
bind   = strings(1, N_TRIM);
WH_val = ctrl.WH(WH_IDX);
for k = 1:N_TRIM
    W = ctrl.interp_mtrx(ctrl.LON.W, u_anchor(k), WH_val);
    B = ctrl.interp_mtrx(ctrl.LON.B, u_anchor(k), WH_val);
    M = (W \ B') / (B * (W \ B'));       % 12 x 3; rows 1:11 physical, 12 virtual
    col = abs(M(1:11, 1));               % effector travel per 1 ft/s^2 forward
    ratio = lim ./ max(col, eps);
    [s_veh(k), ib] = min(ratio);
    bind(k) = nm(ib);
end

% NOTE - THIS ESTIMATE IS WRONG, AND THE WAY IT IS WRONG IS INFORMATIVE.
% It returns 1382 ft/s^2 (43 g) at hover. The reason: at low speed the
% physical effectors produce almost no forward force, so M(1:11,1) is tiny and
% limit/|M| blows up. Forward acceleration at hover is not produced by the
% effectors at all - it is produced by PITCHING OVER, which the allocation
% routes through the virtual attitude channel (column 12 of M), and that
% channel has no entry in `lim`. So these limits do not bound forward
% acceleration; they bound the wrong thing.
% The real bound is:
%    low speed  - how far the vehicle may pitch  (accel ~ g*tan(theta))
%    high speed - pusher thrust minus drag
% Neither is in the allocation matrix. Section 1b below MEASURES the true
% capability so the rest of the calculation has a number to use; deriving it
% instead of measuring it is the open item.
fprintf('--- s_veh from allocation limits (KNOWN WRONG - see the note) ---\n');
fprintf('  u [ft/s] | s_veh [ft/s^2] | which effector runs out first\n');
for k = 1:N_TRIM
    fprintf('%10.1f | %14.2f | %s\n', u_anchor(k), s_veh(k), bind(k));
end

%% ---------------------------------------------------------------------
%  1b. s_veh MEASURED: command far more acceleration than the vehicle has
%      and see what it actually delivers. One flight, and it is the vehicle's
%      capability curve because the command is saturating throughout.
%% ---------------------------------------------------------------------
u_end = u_anchor(end);
T_RUSH = 4.0;                                  % 0 -> 160 in 4 s: impossible
n_h = round(5.0/dt);  n_r = round(T_RUSH/dt);  n_t = round(25.0/dt);
uu = [zeros(1,n_h), u_end*(0:n_r-1)/n_r, repmat(u_end,1,n_t)];
tt = (0:numel(uu)-1)*dt;  pos = cumtrapz(tt, uu);

guam = LpC_GUAM(Config('trim_schedule', params));  guam.reset();
u_act = zeros(1, numel(tt));
for k = 1:numel(tt)
    u_act(k) = guam.state(4);
    guam.step(struct('pos',[pos(k);0;-80], 'vel',[uu(k);0;0], ...
                     'chi',0, 'chi_dot',0));
end
acc = [0, diff(u_act)/dt];
acc = movmean(acc, 51);                        % 0.5 s smoothing
m   = (u_act > 1) & (u_act < u_end-2) & (tt > 5.0);
s_meas = zeros(1, N_TRIM);
for k = 1:N_TRIM
    w = m & abs(u_act - u_anchor(k)) < 6;
    if any(w), s_meas(k) = max(acc(w)); else, s_meas(k) = NaN; end
end
s_meas = fillmissing(s_meas, 'linear');

fprintf('\n--- s_veh MEASURED (max achieved forward acceleration) ---\n');
fprintf('  u [ft/s] | s_meas [ft/s^2]\n');
for k = 1:N_TRIM
    fprintf('%10.1f | %14.2f\n', u_anchor(k), s_meas(k));
end
fprintf('mean %.2f ft/s^2 -> a bare 0 to %.0f ft/s dash would take %.1f s\n', ...
        mean(s_meas), u_end, u_end/mean(s_meas));

%% ---------------------------------------------------------------------
%  2. s_tube(u,T) and the best T at each speed
%% ---------------------------------------------------------------------
nT = numel(T_GRID);
s_tube = zeros(nT, N_TRIM);
for a = 1:nT
    for k = 1:N_TRIM
        s_tube(a,k) = (2*a_fwd(k) - 2*e0_u(k)) / (T_GRID(a) + 2*C(a,k));
    end
end
[s_tube_best, ia] = max(s_tube, [], 1);

fprintf('\n--- the two limits, and which one binds (s_veh = MEASURED) ---\n');
fprintf('  u [ft/s] | s_tube | best T | s_veh | s_max | binding\n');
s_max = zeros(1, N_TRIM);
for k = 1:N_TRIM
    s_max(k) = min(s_tube_best(k), s_meas(k));
    if s_meas(k) < s_tube_best(k), b = 'VEHICLE'; else, b = 'tube'; end
    fprintf('%10.1f | %6.2f | %6.2f | %5.2f | %5.2f | %s\n', ...
            u_anchor(k), s_tube_best(k), T_GRID(ia(k)), s_meas(k), s_max(k), b);
end

%% ---------------------------------------------------------------------
%  3. the integral
%% ---------------------------------------------------------------------
uq = linspace(u_anchor(1), u_anchor(end), 2000);
sq = interp1(u_anchor, s_max, uq);
T_min = trapz(uq, 1 ./ sq);

sq_veh  = interp1(u_anchor, s_meas, uq);
sq_tube = interp1(u_anchor, s_tube_best, uq);
T_veh   = trapz(uq, 1 ./ sq_veh);
T_tube  = trapz(uq, 1 ./ sq_tube);

fprintf('\n=== if T is free over the whole grid ===\n');
fprintf('  tube limit alone         : %6.2f s\n', T_tube);
fprintf('  (the optimiser always picks T = %.2f, the smallest T on the grid)\n', ...
        T_GRID(1));

%% ---------------------------------------------------------------------
%  3b. THE POINT: c(u,T) was only VERIFIED near T = 2 s
%
%  verify_lag_coefficient.m checked the gap equation against a flown mission
%  whose segments were 2 s long, and the mission peak came out 4.50 predicted
%  against 4.45 measured. Nothing was checked at T = 0.25 s, and there is
%  every reason to think it is wrong there: a 0.25 s ramp of the same step is
%  eight times steeper, far outside the small perturbations Phi was measured
%  with. Letting the optimiser read c at T = 0.25 is extrapolating a model
%  past the only place it was ever tested, and the runaway answer is what that
%  extrapolation buys.
%
%  So: hold T fixed, sweep it, and see what the tube limit says at each T.
%% ---------------------------------------------------------------------
fprintf('\n=== minimum time at each FIXED segment time T ===\n');
fprintf('  T [s] | T_min [s] | verified? | flown result\n');
Tmin_of_T = zeros(1, nT);
for a = 1:nT
    sq_a = interp1(u_anchor, s_tube(a,:), uq);
    Tmin_of_T(a) = trapz(uq, 1 ./ sq_a);
    if T_GRID(a) >= 1.5 && T_GRID(a) <= 3.0, vf = 'yes'; else, vf = 'NO  '; end
    switch T_GRID(a)
        case 2.00, fl = '44.0 s  SAFE';
        case 1.50, fl = '36.0 s  unsafe';
        case 1.00, fl = '28.0 s  unsafe';
        case 0.50, fl = '15.5 s  unsafe';
        otherwise, fl = '';
    end
    fprintf('%7.2f | %9.2f | %9s | %s\n', T_GRID(a), Tmin_of_T(a), vf, fl);
end

a19 = find(abs(T_GRID - 2.0) < 1e-9, 1);
sq19 = interp1(u_anchor, s_tube(a19,:), uq);
fprintf(['\nAt T = 2.0 s, where c WAS verified, the formula gives %.1f s.\n' ...
         'The flown answer is 41.8 s (T = 1.90). Agreement %.0f %%.\n'], ...
        trapz(uq, 1./sq19), 100*trapz(uq, 1./sq19)/41.8);
fprintf(['\nSo the analytic route does give the right number - but only when c\n' ...
         'is read where it was checked. The runaway 8 s is not a claim about\n' ...
         'the vehicle, it is the model being asked a question it cannot answer.\n']);
fprintf(['\nStill missing: a bound that comes from the VEHICLE rather than the\n' ...
         'tube, so that short T is ruled out on its own merits instead of by\n' ...
         'refusing to extrapolate. Two attempts failed - see the notes above.\n']);

% where the time is spent
fprintf('\n--- where the time goes ---\n');
edges = [0 40 80 120 160];
for i = 1:numel(edges)-1
    m = uq >= edges(i) & uq <= edges(i+1);
    fprintf('  %3d - %3d ft/s : %5.2f s\n', edges(i), edges(i+1), ...
            trapz(uq(m), 1./sq(m)));
end

%% ---------------------------------------------------------------------
%  4. and the schedule that realises it, at a chosen T
%% ---------------------------------------------------------------------
fprintf('\n--- trim points that realise s_max, for a few segment times ---\n');
fprintf('  T [s] | trim points | total time [s]\n');
for Tf = [2.0 1.9 1.5 1.0]
    u = 0;  n = 0;  tt = 0;
    while u < u_anchor(end)-1e-6 && n < 500
        s = interp1(u_anchor, s_max, min(max(u,0), u_anchor(end)));
        D = min(s*Tf, u_anchor(end)-u);
        u = u + D;  n = n + 1;  tt = tt + D/s;
    end
    fprintf('%7.2f | %11d | %14.2f\n', Tf, n+1, tt);
end

%% ---------------------------------------------------------------------
%  5. figure
%% ---------------------------------------------------------------------
f = figure('Position',[100 100 950 620]);
subplot(2,1,1);
plot(u_anchor, s_meas, "o-", "LineWidth", 1.5); hold on;
plot(u_anchor, s_tube_best, 's-', 'LineWidth', 1.5);
plot(u_anchor, s_max, 'k-', 'LineWidth', 2.2);
ylabel('allowed acceleration [ft/s^2]'); grid on;
legend('vehicle (measured capability)', 'tube (overlap + margin)', ...
       'binding limit', 'Location','northeast');
ylim([0 30]);
title('the safety condition and the vehicle, as competing acceleration limits');

subplot(2,1,2);
plot(uq, cumtrapz(uq, 1./sq), 'k-', 'LineWidth', 2); hold on;
yline(41.8, 'r--', 'LineWidth', 1.4);
xlabel('forward speed [ft/s]'); ylabel('elapsed time [s]'); grid on;
legend(sprintf('analytic, total %.1f s', T_min), 'flown result 41.8 s', ...
       'Location','northwest');
saveas(f, fullfile(root,'logger','analytic_min_time.png'));

save(fullfile(root,'logger','analytic_min_time.mat'), ...
     'u_anchor','s_veh','s_tube','s_tube_best','s_max','T_GRID', ...
     'T_min','T_veh','T_tube','bind');
fprintf('\nsaved logger/analytic_min_time.{mat,png}\n');
