% sweep_q1_fine - find how far q1 can go, with q2 held at its original value.
%
% Separating the channels showed that q1 (accumulated forward-speed error) is
% the one that helps: raising it alone took BRT departure from 25.7 % to
% 11.5 % and improved speed tracking AND altitude at the same time. q2
% (vertical speed) did nothing for BRT and doubled the altitude deviation, so
% it stays at 0.01.
%
% q1 = 0.64 was rejected by the design's own pole check, so the ceiling is
% between 0.16 and 0.64. This walks that interval and reports which trim
% points fail first.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
Q1 = [0.05 0.06 0.07 0.08 0.09 0.10];
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

fprintf('\n%7s | %6s %6s %7s | %8s %8s %8s %8s | %s\n', ...
        'q1','Ki','Kx','1/Kx[s]','|u-r|','e_pos','BRT>=0','alt dev','design');
fprintf('%s\n', repmat('-', 1, 88));
fprintf('%7.3f | %6.3f %6.3f %7.2f | %8.2f %8.1f %8.1f %8.1f | %s\n', ...
        0.010, 0.099, 0.381, 2.63, 11.30, 38.3, 25.7, 58.2, 'ok (shipped)');

for a = 1:numel(Q1)
    Q = [Q1(a) 0.01 1000 0 0 0]';
    [Ki, Kx, badlist] = design(bc0, cfg, Q);
    if ~isempty(badlist)
        fprintf('%7.3f | %6s %6s %7s | %8s %8s %8s %8s | unstable at UH %s\n', ...
                Q1(a), '--','--','--','--','--','--','--', ...
                mat2str(unique(badlist)));
        continue;
    end
    [eu, ep, viol, ad, blown] = mission(params, brtV, dt, Ki, Kx);
    kI = Ki(1,1,8,2);  kX = Kx(1,1,8,2);
    if blown
        fprintf('%7.3f | %6.3f %6.3f %7.2f | %8s %8s %8s %8s | DIVERGED in sim\n', ...
                Q1(a), kI, kX, 1/kX, '--','--','--','--');
    else
        fprintf('%7.3f | %6.3f %6.3f %7.2f | %8.2f %8.1f %8.1f %8.1f | ok\n', ...
                Q1(a), kI, kX, 1/kX, eu, ep, viol, ad);
    end
end

% =========================================================================
function [Ki, Kx, badlist] = design(bc, cfg, Q)
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);  badlist = [];
for jj = 1:M
    for ii = 1:N
        tp = bc.trim_xu_eq(bc.XU0(:,ii,jj));
        [lon, err] = bc.ctrl_lon(tp, Q, cfg.Rlon0, cfg.Wlon0);
        if err, badlist(end+1) = ii; end %#ok<AGROW>
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end
end

function [eu, ep, viol, ad, blown] = mission(params, brtV, dt, Ki, Kx)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;
guam.reset();
rt = guam.refTraj;  N = size(rt.pos,2);
eu = 0; ep = 0; nv = 0; ad = 0; blown = false;
for i = 1:N
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, blown = true; break; end
    R = RSLQR.rotm_i2b(s(7), s(8), s(9));
    e = R * (rt.pos(:,i) - s(1:3));
    eu = max(eu, abs(s(4) - rt.vel(1,i)));
    ep = max(ep, abs(e(1)));
    ad = max(ad, abs(-s(3) - 80));
    nv = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
viol = 100*nv/N;
end
