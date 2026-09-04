% sweep_q4 - does proportional action on speed damp the slow altitude/speed
% oscillation?
%
% With q1 = 0.16 the loop tracks u well in the middle of the mission but
% altitude swings 80 -> 50 -> 130 -> 40 ft with a period around 25 s, and the
% pitch attitude swings with it. That period matches the classic slow
% speed/altitude exchange mode (pi*sqrt(2)*V/g ~ 21 s at 150 ft/s), which is
% normally damped by proportional feedback on speed - here that is Kx, and
% Kx = sqrt(q4 + 2*Ki), so q4 is the knob.
%
% Earlier I argued from a reduced model that q4 would hurt ramp tracking.
% That argument was not measured and the reduced models have disagreed with
% each other, so this just measures it.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  Q1 = 0.16;
Q4 = [0 0.1 0.5 2 8];
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));

g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

fprintf('\nq1 = %.2f fixed, sweeping q4\n\n', Q1);
fprintf('%6s | %6s %6s | %9s %9s %9s | %10s %10s\n', ...
        'q4','Ki','Kx','max|u-r|','e_pos','alt dev','alt swing','max|dth|');
fprintf('%s\n', repmat('-', 1, 76));

for a = 1:numel(Q4)
    Q = [Q1 0.01 1000 Q4(a) 0 0]';
    N = cfg.N_trim;  M = cfg.M_trim;
    Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);  bad = false;
    for jj = 1:M
        for ii = 1:N
            [lon, err] = bc0.ctrl_lon(bc0.trim_xu_eq(bc0.XU0(:,ii,jj)), ...
                                      Q, cfg.Rlon0, cfg.Wlon0);
            if err, bad = true; end
            Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
        end
    end
    kI = Ki(1,1,8,2);  kX = Kx(1,1,8,2);
    if bad
        fprintf('%6.2f | %6.3f %6.3f | %s\n', Q4(a), kI, kX, 'UNSTABLE at design');
        continue;
    end

    S = run_case(params, traj, Ki, Kx, dt);
    if S.blown
        fprintf('%6.2f | %6.3f %6.3f | %s\n', Q4(a), kI, kX, 'DIVERGED');
    else
        fprintf('%6.2f | %6.3f %6.3f | %9.2f %9.1f %9.1f | %10.1f %10.2f\n', ...
                Q4(a), kI, kX, S.eu, S.ep, S.ad, S.swing, S.eth);
    end
end

fprintf(['\nalt swing = peak-to-peak altitude after t = 15 s, i.e. the slow\n' ...
         'oscillation with the launch transient excluded.\n']);

% =========================================================================
function S = run_case(params, traj, Ki, Kx, dt)
guam = LpC_GUAM(Config('trim_schedule', params));
guam.controller.baseline_controller.LON.Ki = Ki;
guam.controller.baseline_controller.LON.Kx = Kx;
guam.reset();
rt = guam.refTraj;  n = size(rt.pos,2);
t = rt.time;  u = zeros(1,n); z = zeros(1,n); th = zeros(1,n);
thc = zeros(1,n); ep = zeros(1,n);  S.blown = false;
for i = 1:n
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, S.blown = true; return; end
    [X0, ~] = guam.controller.interp_xu0(rt.vel(1,i), rt.vel(3,i));
    R = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e = R * (rt.pos(:,i) - s(1:3));
    u(i)=s(4); z(i)=s(3); th(i)=s(8); thc(i)=X0(11); ep(i)=e(1);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
keep = t <= traj.t_node(end) + 3;
t=t(keep); u=u(keep); z=z(keep); th=th(keep); thc=thc(keep); ep=ep(keep);
alt = -z;
late = t >= 15;
S.eu    = max(abs(u - rt.vel(1,keep)));
S.ep    = max(abs(ep));
S.ad    = max(abs(alt - 80));
S.swing = max(alt(late)) - min(alt(late));
S.eth   = max(abs(rad2deg(th - thc)));
end
