% sweep_tseg_rg - how fast can the transition actually be flown?
%
% T_seg only shapes the WISH r(t); the governor decides what is applied. So
% shrinking T_seg does not shorten the mission by itself - it just makes the
% plan more impatient and leaves more of the work to the governor. Pushed far
% enough (T_seg -> 0) the plan degenerates to "go to 160 ft/s right now", and
% what comes out is the fastest transition the vehicle can fly while keeping
% the liveness certificate.
%
% Without a governor the same sweep diverged below T_seg = 1.25 s
% (tests/sweep_tseg.m). The question here is whether the governor removes that
% cliff, and where the completion time bottoms out.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

T_SEGS  = [2.0 1.5 1.0 0.5];
dt_sim  = 0.01;
GOV_HZ  = 10;   N_dec = round(1/(GOV_HZ*dt_sim));
SIM_MAX = 120;

fprintf('\nT_seg | plan  | flown | delay  | not-followed | worst V | alt dev | verdict\n');
fprintf('  [s] |  [s]  |  [s]  |  [s]   |     [%%]      |         |  [ft]   |\n');
fprintf('%s\n', repmat('-', 1, 88));

R = struct('T_seg', {}, 't_done', {}, 'L', {});
for ii = 1:numel(T_SEGS)
    TS = T_SEGS(ii);
    params = struct('filter_mode', 'off');  params.T_seg = TS;
    traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
    brtV = BRTValue(fullfile(root, 'data'), traj.trim_lon, traj.wh_idx);
    cfg  = Config('trim_schedule', params);
    guam = LpC_GUAM(cfg);
    gov  = TrimRefGovernor(traj, brtV, guam, dt_sim, ...
             struct('T', 6.0, 'delta', 0.3, 'eps', 0.02, 'M', 10));
    guam.reset();  gov.reset();

    Ng = round(SIM_MAX/dt_sim);
    tv = zeros(1,Ng);  vv = tv;  kv = tv;  Vv = tv;  alt = tv;
    v_final = traj.trim_lon(1,end);
    t_done = NaN;  diverged = false;
    for i = 1:Ng
        tt = (i-1)*dt_sim;
        if mod(i-1, N_dec) == 0, [ref, gi] = gov.step(tt); kk = gi.kappa;
        else,                     ref = gov.hold();
        end
        tv(i) = tt;  vv(i) = gov.v;  kv(i) = kk;
        Vv(i) = brtV.value(guam.state([4 6 11 8]), gov.v);
        alt(i) = -guam.state(3);
        if isnan(t_done) && gov.v >= v_final - 1e-6, t_done = tt; end
        guam.step(ref);
        if abs(-guam.state(3)-80) > 500 || guam.state(4) > 300
            diverged = true;  break;
        end
        if ~isnan(t_done) && tt > t_done + 2, break; end
    end
    m = 1:i;
    if diverged,            verdict = 'DIVERGED';
    elseif isnan(t_done),   verdict = 'timeout';
    else,                   verdict = 'ok';
    end
    fprintf('%5.2f | %5.1f | %5.2f | %+6.2f |    %5.1f     | %+6.3f |  %5.1f  | %s\n', ...
        TS, traj.t_node(end), t_done, t_done - traj.t_node(end), ...
        100*nnz(kv(m) < 1-1e-9)/numel(m), max(Vv(m)), ...
        max(abs(alt(m)-80)), verdict);

    R(ii).T_seg = TS;  R(ii).t_done = t_done;
    R(ii).L = struct('t', tv(m), 'v', vv(m), 'V', Vv(m));
end

ok = ~isnan([R.t_done]);
if any(ok)
    [tb, jb] = min([R(ok).t_done]);
    Rok = R(ok);
    fprintf(['\nfastest completion: %.2f s at T_seg = %.2f s\n' ...
             'plan time at that T_seg was %.1f s, so the governor gave back %.1f s\n'], ...
            tb, Rok(jb).T_seg, 5 + 19*Rok(jb).T_seg, tb - (5 + 19*Rok(jb).T_seg));
end

%% Figure: the applied command for each schedule aggressiveness
f = figure('Name','T_seg sweep under the governor','Position',[60 60 1150 640]);
hold on; grid on;
col = turbo(numel(R));
for ii = 1:numel(R)
    plot(R(ii).L.t, R(ii).L.v, '-', 'Color', col(ii,:), 'LineWidth', 1.6, ...
         'DisplayName', sprintf('T_{seg} = %.2f s  ->  %.1f s', ...
                                R(ii).T_seg, R(ii).t_done));
end
xlabel('time [s]'); ylabel('applied command v [ft/s]');
title(['Applied command for increasingly impatient schedules ' ...
       '(the governor converges to the same envelope)']);
legend('Location','southeast');
