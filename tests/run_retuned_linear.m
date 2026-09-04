% run_retuned_linear - the linear-prediction governor on the RETUNED plant.
%
% gain_scale = 1.5 cuts the inner-loop time constant from 0.91 s to 0.66 s,
% which drops the ungoverned BRT violation from 25.7 % to 7.3 %. With less
% for the governor to clean up, both the required margin and the mission time
% should fall. Phi describes the CLOSED loop, so it has to be rebuilt for the
% new gains before any of this means anything.
%
% Phase 1  rebuild the closed-loop model at gain_scale = 1.5
% Phase 2  sweep horizon T against margin, find the safe settings
% Phase 3  take the best pair and push T_seg down for the fastest transition
%
% Reference, old gains: nonlinear 51.70 s / 0.0 %; linear+margin 62.30 s / 0.0 %
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  N_dec = 10;  SIM_MAX = 140;
GS = 1.5;                                  % retuned inner-loop gain scale

mk = @(p) local_make(root, p, GS);

%% ---------------- Phase 1: rebuild Phi at the new gains ----------------
params = struct('filter_mode','off');  params.T_seg = 2.0;
traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
guam = mk(params);

h0 = zeros(ClosedLoopModel.N_XI,1);
h0(1:3)=0.5; h0(4:6)=0.1; h0(7:9)=1e-3; h0(10:12)=1e-3;
h0(13:21)=0.5; h0(22:26)=1e-3; h0(27:32)=1e-2; h0(33:34)=1e-3;

fprintf('=== Phase 1: closed-loop model at gain_scale = %.1f ===\n', GS);
M = ClosedLoopModel(dt);
tic; M.build(guam, traj.trim_lon(1,:), h0, 15.0, false); tb = toc;
lam = cellfun(@(P) max(abs(eig(P))), M.Phi);
fprintf('built in %.0f s | max|eig| = %.6f | %s\n\n', tb, max(lam), ...
        ternary(max(lam) < 1, 'all Schur', 'NOT SCHUR - stop here'));
save(fullfile(root,'logger','closed_loop_model_gs15.mat'), 'M', 'GS', '-v7.3');

%% ---------------- Phase 2: T vs margin ----------------
fprintf('=== Phase 2: horizon vs margin (T_seg = 2.0) ===\n');
fprintf('   T | margin | flown  | V>=0 [%%] | worst V | alt dev | wall | verdict\n');
fprintf('%s\n', repmat('-',1,74));
best = struct('t', inf);
for Th = [4.0 6.0 8.0]
  for mg = [0 0.02 0.05 0.10 0.15]
    R = run_once(root, params, GS, brtV, traj, M, Th, mg, dt, N_dec, SIM_MAX);
    if isnan(R.t_done)
        fprintf('%4.1f | %6.2f |   --   | %8.1f | %+7.3f |    --   | %4.1f | TIMEOUT\n', ...
                Th, mg, R.viol, R.Vw, R.wall);
    else
        if R.Vw >= 0, vd = 'UNSAFE'; else, vd = 'ok'; end
        fprintf('%4.1f | %6.2f | %6.2f | %8.1f | %+7.3f | %7.1f | %4.1f | %s\n', ...
                Th, mg, R.t_done, R.viol, R.Vw, R.ad, R.wall, vd);
        if R.Vw < 0 && R.t_done < best.t
            best = struct('t',R.t_done,'T',Th,'mg',mg,'ad',R.ad,'V',R.Vw);
        end
    end
  end
  fprintf('%s\n', repmat('-',1,74));
end
if ~isfinite(best.t)
    fprintf('\nno safe setting found - stopping\n');  return;
end
fprintf('\nbest safe: T = %.1f, margin = %.2f -> %.2f s (old gains gave 62.30 s)\n\n', ...
        best.T, best.mg, best.t);

%% ---------------- Phase 3: push T_seg down ----------------
fprintf('=== Phase 3: T_seg sweep at T = %.1f, margin = %.2f ===\n', best.T, best.mg);
fprintf('T_seg | plan  | flown  | V>=0 [%%] | worst V | alt dev | verdict\n');
fprintf('%s\n', repmat('-',1,66));
fastest = struct('t', inf);
for ts = [2.0 1.5 1.0 0.5 0.25]
    p2 = params;  p2.T_seg = ts;
    tj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, p2));
    R = run_once(root, p2, GS, brtV, tj, M, best.T, best.mg, dt, N_dec, SIM_MAX);
    if isnan(R.t_done)
        fprintf('%5.2f | %5.1f |   --   | %8.1f | %+7.3f |    --   | TIMEOUT\n', ...
                ts, tj.t_node(end), R.viol, R.Vw);
    else
        if R.Vw >= 0, vd = 'UNSAFE'; else, vd = 'ok'; end
        fprintf('%5.2f | %5.1f | %6.2f | %8.1f | %+7.3f | %7.1f | %s\n', ...
                ts, tj.t_node(end), R.t_done, R.viol, R.Vw, R.ad, vd);
        if R.Vw < 0 && R.t_done < fastest.t
            fastest = struct('t',R.t_done,'ts',ts,'ad',R.ad);
        end
    end
end
if isfinite(fastest.t)
    fprintf(['\nfastest safe transition: %.2f s at T_seg = %.2f\n' ...
             '   old gains, nonlinear rollout: 51.70 s\n' ...
             '   old gains, linear + margin  : 62.30 s\n'], fastest.t, fastest.ts);
end

% =========================================================================
function g = local_make(root, params, GS) %#ok<INUSL>
g = LpC_GUAM(Config('trim_schedule', params));
g.controller.gain_scale = GS;
end

function R = run_once(root, params, GS, brtV, traj, M, Th, mg, dt, N_dec, SIM_MAX)
guam = local_make(root, params, GS);
gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
        struct('T',Th,'delta',0.3,'eps',0.02,'M',40,'lin',M,'margin',mg));
guam.reset();  gov.reset();
v_final = traj.trim_lon(1,end);  t_done = NaN;
Ng = round(SIM_MAX/dt);  Vw = -inf;  ad = 0;  nv = 0;  n = 0;
tic;
for i = 1:Ng
    tt = (i-1)*dt;
    if mod(i-1,N_dec)==0, ref = gov.step(tt); else, ref = gov.hold(); end
    Vn = brtV.value(guam.state([4 6 11 8]), gov.v);
    Vw = max(Vw,Vn);  nv = nv + (Vn>=0);  n = n+1;
    ad = max(ad, abs(-guam.state(3) - 80));
    if isnan(t_done) && gov.v >= v_final-1e-6, t_done = tt; end
    guam.step(ref);
    if ~isnan(t_done) && tt > t_done+2, break; end
end
R = struct('t_done',t_done,'Vw',Vw,'ad',ad,'viol',100*nv/n,'wall',toc);
end

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
