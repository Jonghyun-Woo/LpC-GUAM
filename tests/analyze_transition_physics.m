% analyze_transition_physics - how the GUAM Lift+Cruise actually transitions.
%
% At each of the 20 scheduled trim points the aero model is called three times
%   (a) airframe only          : all 9 propellers set to zero
%   (b) airframe + lift rotors : propellers 1..8 at their trim speed
%   (c) everything             : pusher 9 added
% so the vertical force can be split into wing/tail/fuselage, lift rotors and
% pusher. That split IS the transition: the load handing over from thrust-borne
% to wing-borne.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
guam = LpC_GUAM(Config('trim_schedule', params));

W = guam.vehicleConfig.mass * 32.174;          % weight [lbf]
S = guam.vehicleConfig.S;
fprintf('mass %.1f slug -> weight %.0f lbf, wing area %.0f ft^2, span %.1f ft\n\n', ...
        guam.vehicleConfig.mass, W, S, guam.vehicleConfig.b);

T = traj.trims;                    % 25 x 20 = [X0(12); U0(13)]
n = size(T, 2);
R = struct();
fn = {'u','Vair','alpha','theta','q_bar','lift_rpm','push_rpm', ...
      'Fz_air','Fz_lift','Fz_push','Fx_air','Fx_lift','Fx_push'};
for a = 1:numel(fn), R.(fn{a}) = zeros(1, n); end

for k = 1:n
    X0 = T(1:12, k);   U0 = T(13:25, k);
    % XU0 state order is [u v w p q r x y z phi theta psi]
    st = zeros(12,1);
    st(1:3)   = [0; 0; -80];
    st(4:6)   = X0(1:3);
    st(7:9)   = X0(10:12);
    st(10:12) = X0(4:6);

    eng   = U0(5:13);
    srf   = [U0(1)-U0(2); U0(1)+U0(2); U0(3); U0(3); U0(4)];
    [rho, sos] = guam.environment.atmosphere(-st(3));
    Rib   = RSLQR.rotm_i2b(st(7), st(8), st(9));
    v_air = guam.environment.airspeed_body(st(4:6), Rib);
    x_aero = [v_air; st(10:12)];

    F0 = run_LpC_aero(x_aero, zeros(9,1),            srf, rho, sos, guam.units, guam.vehicleConfig);
    F1 = run_LpC_aero(x_aero, [eng(1:8); 0],         srf, rho, sos, guam.units, guam.vehicleConfig);
    F2 = run_LpC_aero(x_aero, eng,                   srf, rho, sos, guam.units, guam.vehicleConfig);

    R.u(k)     = X0(1);
    R.Vair(k)  = norm(v_air);
    R.alpha(k) = atan2(v_air(3), v_air(1));
    R.theta(k) = X0(11);
    R.q_bar(k) = 0.5 * rho * norm(v_air)^2;
    R.lift_rpm(k) = mean(eng(1:8)) * 60/(2*pi);
    R.push_rpm(k) = eng(9) * 60/(2*pi);
    R.Fz_air(k)  = F0(3);            R.Fx_air(k)  = F0(1);
    R.Fz_lift(k) = F1(3) - F0(3);    R.Fx_lift(k) = F1(1) - F0(1);
    R.Fz_push(k) = F2(3) - F1(3);    R.Fx_push(k) = F2(1) - F1(1);
end

%% Vertical load share (body -Z is up; report each as a fraction of weight)
fprintf(' k |   u   | V_air | alpha | theta | q_bar |  lift rpm | push rpm |  wing  | rotors | pusher\n');
fprintf('   | [f/s] | [f/s] | [deg] | [deg] | [psf] |           |          |  %%W    |  %%W    |  %%W\n');
fprintf('%s\n', repmat('-',1,108));
for k = 1:n
    fprintf('%2d | %5.1f | %5.1f | %+5.2f | %+5.2f | %5.2f | %9.0f | %8.0f | %+6.1f | %+6.1f | %+6.1f\n', ...
        k, R.u(k), R.Vair(k), rad2deg(R.alpha(k)), rad2deg(R.theta(k)), R.q_bar(k), ...
        R.lift_rpm(k), R.push_rpm(k), ...
        -100*R.Fz_air(k)/W, -100*R.Fz_lift(k)/W, -100*R.Fz_push(k)/W);
end

%% Longitudinal force share
fprintf('\n k |   u   |  wing Fx | rotors Fx | pusher Fx |  total Fx  [lbf]\n');
fprintf('%s\n', repmat('-',1,66));
for k = 1:n
    fprintf('%2d | %5.1f | %8.0f | %9.0f | %9.0f | %8.0f\n', k, R.u(k), ...
        R.Fx_air(k), R.Fx_lift(k), R.Fx_push(k), ...
        R.Fx_air(k)+R.Fx_lift(k)+R.Fx_push(k));
end

%% Where does the wing take over?
share = -R.Fz_air / W;
kx = find(share >= 0.5, 1, 'first');
if ~isempty(kx)
    fprintf('\nwing carries 50%% of the weight at trim %d, u = %.1f ft/s\n', kx, R.u(kx));
end
fprintf('wing share: %.1f%% at hover -> %.1f%% at cruise\n', 100*share(1), 100*share(end));
save(fullfile(root,'logger','transition_physics.mat'), 'R', 'W', 'S');
