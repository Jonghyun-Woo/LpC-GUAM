% run_transition_sim - Flat-earth m-code port of the GUAM hover-to-cruise
% transition demo. Common to both scenarios:
%   0-20 s : vertical climb 0 -> 80 ft (initial climb rate 8 ft/s)
%   20-40 s: accelerate to 15 ft/s forward flight
% Cruise altitude depends on the scenario:
%   'althold' (default) : hold 80 ft through cruise
%   'climb'             : keep climbing to 100 ft (original)
%
% Usage: run this script from anywhere; it adds Refactoring/ to the path,
% runs the closed-loop simulation, and plots results vs. the reference.

here = fileparts(mfilename('fullpath'));
addpath(genpath(pwd));

scenario = 'althold';   % 'althold' (cruise altitude hold, default) | 'climb' (original)
cfg.M          = 4000;  % sim steps -> T = M*dt = 40 s (0-20 s climb, 20-40 s cruise)
cfg.target_vel = 15;    % cruise forward speed [ft/s] (unused by 'althold' reference build)
guam = LpC_GUAM(scenario, cfg);

% Transition demo runs the nominal RSLQR allocation only (no CBF liveness
% intervention). The allocation frame code in RSLQR.pseudo_alloc still needs
% a valid WH anchor even with the filter off, so set both here.
guam.controller.liveness_lon.mode  = 'off';
guam.controller.liveness_wh_anchor = guam.controller.WH(3);

out  = guam.run();

%% Plots
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
