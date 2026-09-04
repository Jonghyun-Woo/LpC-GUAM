% tune_nominal_gains - why the nominal loop never tracks the plan, and which
% gain fixes it.
%
% Diagnosis. During a segment the plan is a ramp in speed (a = 8.4 ft/s per
% T_seg). Model the inner velocity loop as first order with time constant tau
% tracking the command  c = r + k_pos*e_pos:
%
%     u_dot  = (c - u)/tau ,      e_pos_dot = r - u
%
% In steady state e_pos_dot = 0 forces u = r, and u_dot = 0 then gives
%
%     e_pos = a*tau / k_pos ,      k_pos*e_pos = a*tau
%
% So the SPEED-COMMAND BIAS is a*tau and does not depend on k_pos at all:
% raising the position gain shrinks the position error and leaves the bias
% untouched. Only a faster inner loop (smaller tau) removes the bias, and the
% only knob for that here is a uniform scale on the scheduled Kx and Ki.
%
% This script runs the UNGOVERNED mission (plan straight into RSLQR) and
% measures both, so the prediction above can be checked rather than asserted.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);

a_plan = (traj.trim_lon(1,2) - traj.trim_lon(1,1)) / params.T_seg;
fprintf('plan acceleration a = %.2f ft/s^2\n', a_plan);

K_POS  = [0.1 0.3 1.0];
G_SCL  = [1.0 1.5 2.0 3.0];

fprintf(['\n k_pos | gain | max|u-r| | max|e_pos| | bias k*e | tau=bias/a | ' ...
         'V>=0 [%%] | alt dev | note\n']);
fprintf('%s\n', repmat('-', 1, 100));

for kp = K_POS
  for gs = G_SCL
    guam = LpC_GUAM(Config('trim_schedule', params));
    guam.controller.k_pos      = kp;
    guam.controller.gain_scale = gs;
    guam.reset();

    rt = guam.refTraj;  N = size(rt.pos,2);
    eu = 0; ep = 0; bias = 0; nv = 0; ad = 0; blown = false;
    for i = 1:N
        s = guam.state;
        if any(~isfinite(s)) || abs(s(4)) > 1e4, blown = true; break; end
        R = RSLQR.rotm_i2b(s(7), s(8), s(9));
        e = R * (rt.pos(:,i) - s(1:3));
        eu   = max(eu,   abs(s(4) - rt.vel(1,i)));
        ep   = max(ep,   abs(e(1)));
        bias = max(bias, abs(kp*e(1)));
        ad   = max(ad,   abs(-s(3) - 80));
        nv   = nv + (brtV.value(s([4 6 11 8]), rt.vel(1,i)) >= 0);
        guam.step(struct('pos',rt.pos(:,i),'vel',rt.vel(:,i), ...
                         'chi',rt.chi(i),'chi_dot',rt.chidot(i)));
    end
    if blown
        fprintf('%6.2f | %4.1f |    --    |     --     |    --    |     --     |    --    |    --   | DIVERGED\n', kp, gs);
    else
        fprintf('%6.2f | %4.1f | %8.2f | %10.1f | %8.2f | %10.2f | %8.1f | %7.1f |\n', ...
                kp, gs, eu, ep, bias, bias/a_plan, 100*nv/N, ad);
    end
  end
  fprintf('%s\n', repmat('-', 1, 100));
end

fprintf(['\nbaseline is k_pos = 0.10, gain = 1.0 (what every result so far used)\n' ...
         'read the tau column: if it is flat across k_pos, the bias really is\n' ...
         'set by the inner loop and the position gain cannot fix it.\n']);
