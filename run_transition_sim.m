% run_transition_sim - Flat-earth m-code port of the GUAM hover-to-cruise
% transition demo. Common to both scenarios:
%   0-20 s : vertical climb 0 -> 80 ft (initial climb rate 8 ft/s)
%   20-40 s: accelerate to 15 ft/s forward flight
% Cruise altitude depends on the scenario:
%   'althold' (default) : hold 80 ft through cruise
%   'climb'             : keep climbing to 100 ft (original)
%
% Pipeline: build a central Config -> construct LpC_GUAM(cfg) -> run the
% explicit closed-loop step loop here -> plot results vs. the reference.
% Run this script from anywhere; it adds the repo to the path.

here = fileparts(mfilename('fullpath'));
addpath(genpath(here));

%% Overriding simulation parameters (optional)
scenario = 'althold';   % TODO: 'althold' or 'climb' or 'lon_brt_verify' (추가예정)
params = struct();
params.steps = 4000;
params.target_vel = 15;

params.filter_mode = 'off';
params.filter_wh_anchor = [];

%% 1) Config assembly (single entry point)
cfg = Config(scenario, params);

%% 2) Dynamics simulation object
guam = LpC_GUAM(cfg);

%% 3) Simulation loop
rt = guam.refTraj;  N = size(rt.pos, 2);  dt = guam.simConfig.dt;
out.time    = (0:N-1) .* dt;
out.state   = zeros(12, N);
out.engine  = zeros(9, N);
out.surface = zeros(5, N);
out.alpha   = zeros(1, N);
out.beta    = zeros(1, N);
out.V       = zeros(1, N);
out.ref_pos = rt.pos;
out.ref_vel = rt.vel;

guam.reset();
for k = 1:N
    ref.pos = rt.pos(:, k);  ref.vel = rt.vel(:, k);
    ref.chi = rt.chi(k);     ref.chi_dot = rt.chidot(k);

    out.state(:, k) = guam.state;
    [out.alpha(k), out.beta(k), out.V(k)] = ...
        guam.aeroFrame.compute(guam.state(4), guam.state(5), guam.state(6));

    [engine, surface] = guam.step(ref);
    out.engine(:, k)  = engine;
    out.surface(:, k) = surface;
end

%% 4) Plots
t = out.time;

figure('Name', 'Position (NED)');
lbl = {'North [ft]', 'East [ft]', 'Down [ft]'};
for i = 1:3
    subplot(3, 1, i);
    plot(t, out.state(i, :), 'b', t, out.ref_pos(i, :), 'r--');
    ylabel(lbl{i}); grid on;
    if i == 1, legend('sim', 'ref'); title('Inertial position'); end
end
xlabel('Time [s]');

figure('Name', 'Body velocity');
lbl = {'u [ft/s]', 'v [ft/s]', 'w [ft/s]'};
for i = 1:3
    subplot(3, 1, i);
    plot(t, out.state(3 + i, :), 'b', t, out.ref_vel(i, :), 'r--');
    ylabel(lbl{i}); grid on;
    if i == 1, legend('sim', 'ref (heading frame)'); title('Velocity'); end
end
xlabel('Time [s]');

figure('Name', 'Attitude');
lbl = {'\phi [deg]', '\theta [deg]', '\psi [deg]'};
for i = 1:3
    subplot(3, 1, i);
    plot(t, rad2deg(out.state(6 + i, :)), 'b');
    ylabel(lbl{i}); grid on;
    if i == 1, title('Euler angles'); end
end
xlabel('Time [s]');

figure('Name', 'Effectors');
subplot(3, 1, 1);
plot(t, out.engine(1:8, :) .* (60 / (2 * pi)));
ylabel('Lift rotors [RPM]'); grid on; title('Effectors');
subplot(3, 1, 2);
plot(t, out.engine(9, :) .* (60 / (2 * pi)));
ylabel('Pusher [RPM]'); grid on;
subplot(3, 1, 3);
plot(t, rad2deg(out.surface));
ylabel('Surfaces [deg]'); grid on;
legend('LA', 'RA', 'LE', 'RE', 'RUD');
xlabel('Time [s]');
