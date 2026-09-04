% run_transition_evidence_fig - time histories that show the vehicle actually
% transitions, not just accelerates.
%
% Four panels over the governed mission:
%   1  speed        : plan r, governed command v, flown u
%   2  propulsion   : pusher rpm against the lift-rotor rows (the handover)
%   3  load share   : vertical force carried by the rotors vs by the airframe,
%                     as a fraction of weight (the load transfer)
%   4  surfaces     : flap retracting and the elevator waking up
%
% Panel 3 is obtained by re-running the aero model at logged samples with the
% propellers zeroed, so the airframe contribution can be separated. The GUAM
% aero database blends prop-on/prop-off effects, so below ~80 ft/s that split
% attributes wing lift generated in the prop wash to the rotor term. Read the
% low-speed end as a bound, not an exact share.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt_sim = 0.01;
params = struct('filter_mode', 'off');  params.T_seg = 2.0;
GOV_HZ = 10;  N_dec = round(1/(GOV_HZ*dt_sim));

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
guam = LpC_GUAM(Config('trim_schedule', params));
gov  = TrimRefGovernor(traj, brtV, guam, dt_sim, ...
        struct('T', 6.0, 'delta', 0.3, 'eps', 0.02, 'M', 40));
guam.reset();  gov.reset();

%% Governed run, logging the effectors
SIM_MAX = 90;  Ng = round(SIM_MAX/dt_sim);
L = struct('t', zeros(1,Ng), 'st', zeros(12,Ng), 'v', zeros(1,Ng), ...
           'r', zeros(1,Ng), 'eng', zeros(9,Ng), 'srf', zeros(5,Ng));
v_final = traj.trim_lon(1,end);  t_done = NaN;  last = [];
tic;
for i = 1:Ng
    tt = (i-1)*dt_sim;
    if mod(i-1, N_dec) == 0, [ref, gi] = gov.step(tt);  last = gi;
    else,                    ref = gov.hold();  end
    L.t(i)   = tt;                L.st(:,i)  = guam.state;
    L.v(i)   = gov.v;             L.r(i)     = last.r;
    L.eng(:,i) = guam.engineDynamics.pos;
    L.srf(:,i) = guam.surfaceDynamics.pos;
    if isnan(t_done) && gov.v >= v_final - 1e-6, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done + 3, break; end
end
n = i;  fn = fieldnames(L);
for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:,1:n); end
L.t_done = t_done;
fprintf('governed run: final command at %.2f s (plan %.2f s), wall %.0f s\n', ...
        t_done, traj.t_node(end), toc);

%% Load split at decimated samples
step = 20;                                   % every 0.2 s
idx  = 1:step:n;
W    = guam.vehicleConfig.mass * 32.174;
Fz_air = zeros(1,numel(idx));  Fz_rot = zeros(1,numel(idx));
Fx_psh = zeros(1,numel(idx));  Fx_rot = zeros(1,numel(idx));
for a = 1:numel(idx)
    j = idx(a);  st = L.st(:,j);
    [rho, sos] = guam.environment.atmosphere(-st(3));
    Rib   = RSLQR.rotm_i2b(st(7), st(8), st(9));
    x_aero = [guam.environment.airspeed_body(st(4:6), Rib); st(10:12)];
    e = L.eng(:,j);  s = L.srf(:,j);
    F0 = run_LpC_aero(x_aero, zeros(9,1),      s, rho, sos, guam.units, guam.vehicleConfig);
    F1 = run_LpC_aero(x_aero, [e(1:8); 0],     s, rho, sos, guam.units, guam.vehicleConfig);
    F2 = run_LpC_aero(x_aero, e,               s, rho, sos, guam.units, guam.vehicleConfig);
    Fz_air(a) = F0(3);          Fz_rot(a) = F1(3) - F0(3);
    Fx_rot(a) = F1(1) - F0(1);  Fx_psh(a) = F2(1) - F1(1);
end
td = L.t(idx);
r2rpm = 60/(2*pi);

%% Figure
f = figure('Name','GUAM transition evidence','Position',[60 30 1150 950]);
ax = gobjects(1,4);

ax(1) = subplot(4,1,1); hold on; grid on;
plot(ax(1), L.t, L.r, '--', 'Color',[.15 .15 .15], 'LineWidth',1.2, 'DisplayName','plan r');
plot(ax(1), L.t, L.v, '-',  'Color',[0 .2 .85], 'LineWidth',1.8, 'DisplayName','governed command v');
plot(ax(1), L.t, L.st(4,:), '-', 'Color',[.85 .2 .2], 'LineWidth',1.5, 'DisplayName','flown u');
ylabel(ax(1),'speed [ft/s]'); legend(ax(1),'Location','southeast');
title(ax(1), sprintf('Hover-to-cruise transition under the governor (final command at %.1f s)', t_done));

ax(2) = subplot(4,1,2); hold on; grid on;
plot(ax(2), L.t, L.eng(9,:)*r2rpm, '-', 'Color',[.85 .35 0], 'LineWidth',2.0, ...
     'DisplayName','pusher prop (9)');
plot(ax(2), L.t, mean(L.eng(1:4,:),1)*r2rpm, '-', 'Color',[0 .45 .75], 'LineWidth',1.6, ...
     'DisplayName','lift rotors, front row (1-4)');
plot(ax(2), L.t, mean(L.eng(5:8,:),1)*r2rpm, '-', 'Color',[.3 .7 .35], 'LineWidth',1.6, ...
     'DisplayName','lift rotors, rear row (5-8)');
ylabel(ax(2),'rotor speed [rpm]'); legend(ax(2),'Location','east');
title(ax(2),'Propulsion handover: the pusher spins up as the lift rotors spin down');

ax(3) = subplot(4,1,3); hold on; grid on;
plot(ax(3), td, -100*Fz_rot/W, '-', 'Color',[0 .45 .75], 'LineWidth',1.8, ...
     'DisplayName','carried by lift rotors');
plot(ax(3), td, -100*Fz_air/W, '-', 'Color',[.55 .15 .6], 'LineWidth',1.8, ...
     'DisplayName','carried by wing + tail');
yline(ax(3), 100, 'k:', 'HandleVisibility','off');
ylabel(ax(3),'share of weight [%]'); legend(ax(3),'Location','east');
title(ax(3),'Load transfer: thrust-borne \rightarrow wing-borne');

ax(4) = subplot(4,1,4); hold on; grid on;
yyaxis(ax(4),'left');
plot(ax(4), L.t, rad2deg(mean(L.srf(1:2,:),1)), '-', 'LineWidth',1.8, ...
     'DisplayName','flap (mean of flaperons)');
ylabel(ax(4),'flap [deg]');
yyaxis(ax(4),'right');
plot(ax(4), L.t, rad2deg(L.srf(3,:)), '-', 'LineWidth',1.8, 'DisplayName','elevator');
ylabel(ax(4),'elevator [deg]');
xlabel(ax(4),'time [s]'); legend(ax(4),'Location','east');
title(ax(4),'Aerodynamic surfaces: flap retracts, elevator takes over pitch');

linkaxes(ax,'x'); xlim(ax(1), [0 L.t(end)]);

out = fullfile(root,'logger','transition_evidence');
savefig(f, [out '.fig']);
exportgraphics(f, [out '.png'], 'Resolution', 150);
save([out '.mat'], 'L', 'td', 'Fz_air', 'Fz_rot', 'Fx_rot', 'Fx_psh', 'W', 'traj', '-v7.3');
fprintf('saved %s(.fig/.png/.mat)\n', out);
