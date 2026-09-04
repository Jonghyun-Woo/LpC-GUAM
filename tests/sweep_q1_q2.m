% sweep_q1_q2 - which weight buys the BRT improvement, and which one wrecks
% the altitude?
%
% The first sweep raised q1 and q2 together (the accumulated errors in
% forward and vertical speed) and got BRT departure 25.7 % -> 16.3 %, but
% altitude deviation went 58 ft -> 122 ft. Those are different channels, so
% they are separated here:
%
%   variant A : raise q1 only   (forward speed)   q2 stays 0.01
%   variant B : raise q2 only   (vertical speed)  q1 stays 0.01
%   variant C : raise both      (what the first sweep did)
%
% Response speed is read off Kx rather than from a step test - the step test
% in the previous sweep produced nonsense (thousands of percent overshoot),
% and 1/Kx is the time scale directly: an error of 1 ft/s commands Kx ft/s^2.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
VALS = [0.01 0.04 0.16];
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

g0  = LpC_GUAM(Config('trim_schedule', params));
bc0 = g0.controller.baseline_controller;
cfg = RSLQRConfig;

fprintf('\n%-8s %6s %6s | %6s %6s | %8s %8s %8s %8s\n', ...
        'variant','q1','q2','Ki','Kx','1/Kx[s]','|u-r|','BRT>=0','alt dev');
fprintf('%s\n', repmat('-', 1, 78));

for variant = 'ABC'
  for a = 1:numel(VALS)
    if variant == 'A', q1 = VALS(a);  q2 = 0.01;
    elseif variant == 'B', q1 = 0.01; q2 = VALS(a);
    else,                  q1 = VALS(a); q2 = VALS(a);  end
    if a > 1 || variant == 'A'      % skip duplicate baselines
        Q = [q1 q2 1000 0 0 0]';
        [Ki, Kx, bad] = design(bc0, cfg, Q);
        if bad
            fprintf('%-8s %6.3f %6.3f | %s\n', variant, q1, q2, 'UNSTABLE');
            continue;
        end
        [eu, viol, ad, blown] = mission(params, brtV, dt, Ki, Kx);
        kI = Ki(1,1,8,2);  kX = Kx(1,1,8,2);
        if blown
            fprintf('%-8s %6.3f %6.3f | %6.3f %6.3f | %8s DIVERGED\n', ...
                    variant, q1, q2, kI, kX, '--');
        else
            fprintf('%-8s %6.3f %6.3f | %6.3f %6.3f | %8.2f %8.2f %8.1f %8.1f\n', ...
                    variant, q1, q2, kI, kX, 1/kX, eu, viol, ad);
        end
    end
  end
  fprintf('%s\n', repmat('-', 1, 78));
end

fprintf(['\nbaseline q1 = q2 = 0.010 : |u-r| 11.30, BRT 25.7 %%, alt dev 58.2 ft\n' ...
         'A raises forward-speed tracking, B raises vertical-speed tracking.\n']);

% =========================================================================
function [Ki, Kx, bad] = design(bc, cfg, Q)
N = cfg.N_trim;  M = cfg.M_trim;
Ki = zeros(3,3,N,M);  Kx = zeros(3,4,N,M);  bad = false;
for jj = 1:M
    for ii = 1:N
        tp = bc.trim_xu_eq(bc.XU0(:,ii,jj));
        [lon, err] = bc.ctrl_lon(tp, Q, cfg.Rlon0, cfg.Wlon0);
        if err, bad = true;  return;  end
        Ki(:,:,ii,jj) = lon.Ki;  Kx(:,:,ii,jj) = lon.Kx;
    end
end
end

function [eu, viol, ad, blown] = mission(params, brtV, dt, Ki, Kx)
guam = LpC_GUAM(Config('trim_schedule', params));
bc = guam.controller.baseline_controller;
bc.LON.Ki = Ki;  bc.LON.Kx = Kx;
guam.reset();
rt = guam.refTraj;  N = size(rt.pos,2);
eu = 0; nv = 0; ad = 0; blown = false;
for i = 1:N
    s = guam.state;
    if any(~isfinite(s)) || abs(s(4)) > 1e4, blown = true; break; end
    eu = max(eu, abs(s(4) - rt.vel(1,i)));
    ad = max(ad, abs(-s(3) - 80));
    nv = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
    guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                     'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
end
viol = 100*nv/N;
end
