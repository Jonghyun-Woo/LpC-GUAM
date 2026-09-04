% check_alloc_theta - did raising the allocation weight on theta actually do
% anything, and if not, why not?
%
% Prints, at a low-speed trim and a cruise trim:
%   - the forward-acceleration row of the effectiveness matrix B, so it is
%     visible which effectors can produce u_dot at all
%   - the theta row of the allocation M for a range of weights, i.e. how much
%     theta the allocator asks for per unit of demanded acceleration
clear;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

params = struct('filter_mode','off');  params.T_seg = 2.0;
g0  = LpC_GUAM(Config('trim_schedule', params));
bc  = g0.controller.baseline_controller;
cfg = RSLQRConfig;

fprintf('\nstored W(12,12) before = %g\n', bc.LON.W(12,12,1,1));
bc.LON.W(12,12,:,:) = bc.LON.W(12,12,:,:) * 200;
fprintf('stored W(12,12) after  = %g   (assignment %s)\n\n', ...
        bc.LON.W(12,12,1,1), ...
        string(abs(bc.LON.W(12,12,1,1) - 20) < 1e-9));

NAMES = {'omp1','omp2','omp3','omp4','omp5','omp6','omp7','omp8', ...
         'pusher','elev','flap','THETA'};
IDX = [1 12];  LBL = {'hover-ish trim (i=1)','cruise trim (i=12)'};

for a = 1:2
    i = IDX(a);
    lon = bc.ctrl_lon(bc.trim_xu_eq(bc.XU0(:,i,1)), cfg.Qlon0, cfg.Rlon0, cfg.Wlon0);
    B = lon.B;   % 3 x 12, rows = [u_dot; w_dot; q_dot]
    fprintf('=== %s ===\n', LBL{a});
    fprintf('forward-acceleration authority  dU_dot/d(effector):\n');
    for k = 1:12
        fprintf('   %-8s %12.4f', NAMES{k}, B(1,k));
        if mod(k,3)==0, fprintf('\n'); end
    end
    fprintf('\ntheta asked for per unit demand, as the weight on theta changes:\n');
    fprintf('   %8s | %10s %10s %10s | %10s\n', ...
            'w_theta','d/dUdot','d/dWdot','d/dQdot','pusher d/dUdot');
    for wt = [0.1 1 20 1000 1e6]
        W = lon.W;  W(12,12) = wt;
        Mx = (W\lon.B')/(lon.B*(W\lon.B'));
        fprintf('   %8.3g | %10.4f %10.4f %10.4f | %10.4f\n', ...
                wt, Mx(12,1), Mx(12,2), Mx(12,3), Mx(9,1));
    end
    fprintf('\n');
end
