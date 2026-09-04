% diag_osc_source - is the slow altitude/speed oscillation the outer position
% loop, or a mode of the airframe?
%
% Test: settle at a fixed speed, displace the altitude, release, and watch the
% free response with the command held. Sweep the outer position-loop gain.
%
%   period moves with the gain -> the loop is making the oscillation
%   period stays put           -> it is an airframe mode the loop only rides
%
% The mission run cannot answer this: it lasts ~45 s and the oscillation has a
% ~25 s period, so barely one cycle is visible and the launch transient sits
% on top of it. Here the vehicle is settled first and then given a clean
% displacement, so several cycles can be counted.
%
% k_pos = 0.1 is the shipped value and is left as the default in RSLQR; this
% only sets it on the local copies.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
KP     = [0.05 0.10 0.20 0.40];   % 0.10 is the shipped value
V_TEST = 100;                     % hold this forward speed [ft/s]
DZ     = 20;                      % altitude displacement [ft]
T_FREE = 150;                     % watch this long [s]

params = struct('filter_mode','off');  params.T_seg = 2.0;
g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

% tuned longitudinal gains from the Q sweep
Q = [0.08 0.01 1000 2.0 0 0]';
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);
for jj = 1:M
    for ii = 1:N
        lon = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end

fprintf('\nhold u = %d ft/s, displace altitude %+d ft, release\n', V_TEST, DZ);
fprintf('q1 = 0.08, q4 = 2.0\n\n');
fprintf('%7s | %8s %8s | %10s %10s | %s\n', ...
        'k_pos','period','damping','peak 1','peak 2','peak times [s]');
fprintf('%s\n', repmat('-', 1, 74));

S = cell(1,numel(KP));
for a = 1:numel(KP)
    S{a} = free_response(params, Ki, Kx, KP(a), V_TEST, DZ, dt, T_FREE);
    r = S{a};
    if isnan(r.period)
        fprintf('%7.2f | %8s %8s | %10.1f %10s | %s\n', KP(a), ...
                'n/a','n/a', r.p1, 'n/a', 'fewer than 2 peaks');
    else
        fprintf('%7.2f | %8.1f %8.3f | %10.1f %10.1f | %s\n', KP(a), ...
                r.period, r.zeta, r.p1, r.p2, mat2str(round(r.tp(1:min(4,end)),1)));
    end
end

fprintf(['\nperiod is the mean gap between successive altitude peaks;\n' ...
         'damping is estimated from the ratio of consecutive peak heights.\n']);

%% figure
f = figure('Name','free response after an altitude displacement', ...
           'Position',[60 60 1100 720]);
ax(1) = subplot(2,1,1); hold on; grid on;
ax(2) = subplot(2,1,2); hold on; grid on;
co = lines(numel(KP));
for a = 1:numel(KP)
    nm = sprintf('k_{pos} = %.2f', KP(a));
    if abs(KP(a)-0.1) < 1e-9, nm = [nm '  (shipped)']; end %#ok<AGROW>
    plot(ax(1), S{a}.t, S{a}.alt - S{a}.alt0, '-', 'Color', co(a,:), ...
         'LineWidth', 1.5, 'DisplayName', nm);
    plot(ax(2), S{a}.t, S{a}.u - V_TEST, '-', 'Color', co(a,:), 'LineWidth', 1.5);
end
yline(ax(1), 0, 'k-', 'HandleVisibility','off');
yline(ax(2), 0, 'k-');
ylabel(ax(1), 'altitude deviation [ft]');  legend(ax(1),'Location','northeast');
ylabel(ax(2), 'speed deviation [ft/s]');   xlabel(ax(2),'time after release [s]');
title(ax(1), sprintf('free response, u held at %d ft/s, altitude displaced %+d ft', ...
                     V_TEST, DZ));
linkaxes(ax,'x');  xlim(ax(1), [0 T_FREE]);
exportgraphics(f, fullfile(root,'logger','osc_source.png'),'Resolution',140);
savefig(f, fullfile(root,'logger','osc_source.fig'));
fprintf('saved logger/osc_source(.fig/.png)\n');

% =========================================================================
function R = free_response(params, Ki, Kx, kp, v, dz, dt, T)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;  bc.k_pos = kp;

% Start AT the test trim, not at hover. reset() initialises from the first
% sample of refTraj, so hand it a constant-speed reference; otherwise the
% position reference runs away at v while the vehicle sits at 0 and the
% position feedback commands thousands of ft/s.
rt = guam.refTraj;  m = size(rt.pos, 2);
rt.time   = (0:m-1)*dt;
rt.vel    = repmat([v; 0; 0], 1, m);
rt.pos    = [(0:m-1)*dt*v; zeros(1,m); -80*ones(1,m)];
rt.chi    = zeros(1,m);
rt.chidot = zeros(1,m);
guam.refTraj = rt;
guam.reset();

mk = @(p,vv) struct('pos',[p;0;-80],'vel',[vv;0;0],'chi',0,'chi_dot',0);

% let any residual settle
p = 0;
for i = 1:round(25/dt), p = p + v*dt;  guam.step(mk(p,v)); end
R.alt0 = -guam.state(3);
if abs(R.alt0 - 80) > 20
    warning('settled %.1f ft off the reference altitude', R.alt0 - 80);
end

% displace the altitude and release
s = guam.state;  s(3) = s(3) - dz;  guam.state = s;

n = round(T/dt);  R.t = (1:n)*dt;  R.alt = zeros(1,n);  R.u = zeros(1,n);
for i = 1:n
    p = p + v*dt;  guam.step(mk(p,v));
    R.alt(i) = -guam.state(3);  R.u(i) = guam.state(4);
end

% peaks of the altitude deviation
d = R.alt - R.alt0;
[pk, ip] = findpeaks_local(d);
R.tp = R.t(ip);
if numel(pk) >= 2
    R.period = mean(diff(R.tp));
    ratio = abs(pk(2)/pk(1));
    lg = log(1/max(ratio, 1e-9));
    R.zeta = lg / sqrt(4*pi^2 + lg^2);
    R.p1 = pk(1);  R.p2 = pk(2);
else
    R.period = NaN;  R.zeta = NaN;
    if isempty(pk), R.p1 = max(d); else, R.p1 = pk(1); end
    R.p2 = NaN;
end
end

function [pk, ip] = findpeaks_local(y)
% local maxima that clear a small threshold, so noise is not counted
thr = 0.02 * max(abs(y));
ip = [];
for i = 2:numel(y)-1
    if y(i) > y(i-1) && y(i) >= y(i+1) && y(i) > thr
        ip(end+1) = i; %#ok<AGROW>
    end
end
pk = y(ip);
end
