% sweep_q1_q4 - q1 and q4 together.
%
% q1 (integral action on forward-speed error) took the speed error from 11.30
% to 5.79 ft/s; q4 (proportional action on the same error) then took it to
% 4.34. Neither was swept with the other moving, so this fills in the grid.
%
% Kx = sqrt(q4 + 2*sqrt(q1)) and Ki = sqrt(q1), so both weights push Kx up
% while only q1 moves Ki. The pair that matters is therefore not obvious from
% the formulas alone.
%
% q2 stays 0.01 (raising it doubles the altitude error), q3 stays 1000 (the
% attitude loop is the fast inner loop of the cascade), q5 = q6 = 0.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q1 = [0.08 0.16 0.24 0.40 0.64];
Q4 = [0 0.5 2.0 4.0];
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

EU = nan(numel(Q1), numel(Q4));  EP = EU;  AD = EU;  SW = EU;  VI = EU;
fprintf('\nq2 = 0.01, q3 = 1000, q5 = q6 = 0\n');
fprintf('%6s %6s | %6s %6s | %8s %8s %8s %9s %8s\n', ...
        'q1','q4','Ki','Kx','max|u-r|','e_pos','alt dev','alt swing','BRT>=0');
fprintf('%s\n', repmat('-', 1, 80));
for a = 1:numel(Q1)
  for b = 1:numel(Q4)
    Q = [Q1(a) 0.01 1000 Q4(b) 0 0]';
    [Ki, Kx, bad] = design(bc0, cfg, Q);
    kI = Ki(1,1,8,2);  kX = Kx(1,1,8,2);
    if bad
        fprintf('%6.2f %6.2f | %6.3f %6.3f | %s\n', Q1(a), Q4(b), kI, kX, 'UNSTABLE');
        continue;
    end
    S = run_case(params, traj, brtV, Ki, Kx);
    if S.blown
        fprintf('%6.2f %6.2f | %6.3f %6.3f | %s\n', Q1(a), Q4(b), kI, kX, 'DIVERGED');
        continue;
    end
    EU(a,b)=S.eu; EP(a,b)=S.ep; AD(a,b)=S.ad; SW(a,b)=S.swing; VI(a,b)=S.viol;
    fprintf('%6.2f %6.2f | %6.3f %6.3f | %8.2f %8.1f %8.1f %9.1f %8.1f\n', ...
            Q1(a), Q4(b), kI, kX, S.eu, S.ep, S.ad, S.swing, S.viol);
  end
  fprintf('%s\n', repmat('-', 1, 80));
end

%% grids
show = @(nm, X) print_grid(nm, X, Q1, Q4);
show('max |u - r|  [ft/s]', EU);
show('max position error  [ft]', EP);
show('altitude swing after 15 s  [ft]', SW);
show('BRT departure  [% of mission]', VI);

[~, i] = min(EU(:));  [ia, ib] = ind2sub(size(EU), i);
fprintf(['\nbest speed tracking : q1 = %.2f, q4 = %.2f  ->  %.2f ft/s\n' ...
         '   e_pos %.1f ft | alt dev %.1f ft | alt swing %.1f ft | BRT %.1f %%\n'], ...
        Q1(ia), Q4(ib), EU(ia,ib), EP(ia,ib), AD(ia,ib), SW(ia,ib), VI(ia,ib));
fprintf('shipped design (q1 = 0.01, q4 = 0): 11.30 ft/s | 38.3 ft | 58.2 ft | 25.7 %%\n');
save(fullfile(root,'logger','sweep_q1_q4.mat'),'Q1','Q4','EU','EP','AD','SW','VI');

% =========================================================================
function print_grid(nm, X, Q1, Q4)
fprintf('\n%s\n%8s', nm, 'q1 \ q4');
fprintf('%9.2f', Q4);  fprintf('\n');
for a = 1:numel(Q1)
    fprintf('%8.2f', Q1(a));
    for b = 1:numel(Q4)
        if isnan(X(a,b)), fprintf('%9s','--'); else, fprintf('%9.2f', X(a,b)); end
    end
    fprintf('\n');
end
end

function [Ki, Kx, bad] = design(bc, cfg, Q)
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);  bad = false;
for jj = 1:M
    for ii = 1:N
        [lon, err] = bc.ctrl_lon(bc.trim_xu_eq(bc.XU0(:,ii,jj)), Q, cfg.Rlon0, cfg.Wlon0);
        if err, bad = true; end
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end
end

function S = run_case(params, traj, brtV, Ki, Kx)
guam = LpC_GUAM(Config('trim_schedule', params));
guam.controller.baseline_controller.LON.Ki = Ki;
guam.controller.baseline_controller.LON.Kx = Kx;
guam.reset();
rt = guam.refTraj;  n = size(rt.pos,2);
t = rt.time;  u = zeros(1,n); z = zeros(1,n); ep = zeros(1,n);  nv = 0;
S.blown = false;
for i = 1:n
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, S.blown = true; return; end
    R = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e = R * (rt.pos(:,i) - s(1:3));
    u(i)=s(4); z(i)=s(3); ep(i)=e(1);
    nv = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
keep = t <= traj.t_node(end) + 3;
alt = -z(keep);  late = t(keep) >= 15;
S.eu    = max(abs(u(keep) - rt.vel(1,keep)));
S.ep    = max(abs(ep(keep)));
S.ad    = max(abs(alt - 80));
S.swing = max(alt(late)) - min(alt(late));
S.viol  = 100*nv/n;
end
