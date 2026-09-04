% sweep_linear_T_margin - horizon and constraint tightening together.
%
% The linear model's error grows with the horizon, so a shorter T should need
% a smaller margin. And the exit diagnostic showed the overshoot that breaks
% the certificate develops within ~2 s, so a short horizon may lose nothing.
% Both knobs pull the same way, hence the joint sweep.
%
% Reference (nonlinear rollout, T = 6, margin = 0):
%   51.70 s | 0.0 % outside | worst V -0.000 | alt dev 29.2 ft | wall 1285 s
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  N_dec = 10;  SIM_MAX = 120;
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
S = load(fullfile(root,'logger','closed_loop_model.mat'), 'M');

T_LIST  = [2.0 3.0 4.0 6.0];
MG_LIST = [0 0.05 0.10 0.15];

fprintf('\n   T | margin | flown  | V>=0 [%%] | worst V | alt dev | wall | verdict\n');
fprintf('%s\n', repmat('-', 1, 76));
best = struct('t', inf);
for Th = T_LIST
  for mg = MG_LIST
    guam = LpC_GUAM(Config('trim_schedule', params));
    gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
            struct('T',Th,'delta',0.3,'eps',0.02,'M',40, ...
                   'lin',S.M,'margin',mg));
    guam.reset();  gov.reset();

    v_final = traj.trim_lon(1,end);  t_done = NaN;
    Ng = round(SIM_MAX/dt);  Vw = -inf;  ad = 0;  nv = 0;  n = 0;
    tic;
    for i = 1:Ng
        tt = (i-1)*dt;
        if mod(i-1,N_dec)==0, ref = gov.step(tt); else, ref = gov.hold(); end
        Vn = brtV.value(guam.state([4 6 11 8]), gov.v);
        Vw = max(Vw, Vn);  nv = nv + (Vn >= 0);  n = n + 1;
        ad = max(ad, abs(-guam.state(3) - 80));
        if isnan(t_done) && gov.v >= v_final-1e-6, t_done = tt; end
        guam.step(ref);
        if ~isnan(t_done) && tt > t_done+2, break; end
    end
    w = toc;

    if isnan(t_done)
        fprintf('%4.1f | %6.2f |   --   | %8.1f | %+7.3f |    --   | %4.1f | TIMEOUT\n', ...
                Th, mg, 100*nv/n, Vw, w);
    else
        if Vw >= 0
            vd = 'UNSAFE';
        else
            vd = 'ok';
            if t_done < best.t
                best = struct('t',t_done,'T',Th,'mg',mg,'ad',ad,'V',Vw,'w',w);
            end
        end
        fprintf('%4.1f | %6.2f | %6.2f | %8.1f | %+7.3f | %7.1f | %4.1f | %s\n', ...
                Th, mg, t_done, 100*nv/n, Vw, ad, w, vd);
    end
  end
  fprintf('%s\n', repmat('-', 1, 76));
end

if isfinite(best.t)
    fprintf(['\nfastest SAFE setting: T = %.1f s, margin = %.2f\n' ...
             '   flown %.2f s | worst V %+.3f | alt dev %.1f ft | wall %.1f s\n'], ...
            best.T, best.mg, best.t, best.V, best.ad, best.w);
    fprintf('   vs nonlinear: %+.2f s slower, %+.1f ft altitude, %.0fx faster\n', ...
            best.t - 51.70, best.ad - 29.2, 1285/best.w);
else
    fprintf('\nno safe setting found in this grid\n');
end
