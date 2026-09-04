% sweep_tseg - how fast can the transition be scheduled and still stay
% BRT-certified under the governor?
%
% T_seg sets the commanded acceleration: each segment covers 8.44 ft/s of
% forward speed, so the schedule demands 8.44/T_seg ft/s^2. Below some
% T_seg the vehicle simply cannot follow, and no reference governor can
% fix that - a governor only reshapes the reference, it has no authority
% over the plant. This sweep locates that boundary empirically.
%
% Divergence guard: the (u,w,q,theta) BRT carries no altitude state, so a
% dive is invisible to it. The loop therefore watches altitude and speed
% directly and aborts the run when the plant is clearly gone.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

T_SEGS = [3.0 2.5 2.0 1.5 1.25 1.0];
dt_sim = 0.01;

spec = FilterConfig.channelSpec('lon');
gv = cell(1, 4);
for d = 1:4
    gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
end
brt = cell(1, 20);
for k = 1:20
    S = load(fullfile(root, 'data', ...
        sprintf('GUAM_LON_BRT_HJIR_UH%d_WH3.mat', k)), 'data');
    brt{k} = S.data;
end

fprintf('\nT_seg | accel   | nominal | reached | delay  | freeze | out-both | alt dev | verdict\n');
fprintf(' [s]  | [ft/s^2]|   [s]   |   [s]   |  [s]   |  [s]   |   [s]    |  [ft]   |\n');
fprintf('%s\n', repmat('-', 1, 92));

R = struct('T_seg', {}, 'reached', {}, 'ok', {});
for ii = 1:numel(T_SEGS)
    TS = T_SEGS(ii);
    params = struct('filter_mode', 'off');  params.T_seg = TS;
    traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt_sim, params));
    cfg  = Config('trim_schedule', params);
    guam = LpC_GUAM(cfg);
    gov  = TrimProgressGovernor(traj, brt, gv, dt_sim, 0.05);
    guam.reset();  gov.reset();

    T_max = 3 * traj.t_node(end) + 20;
    Ng    = round(T_max / dt_sim);
    t_fin = NaN;  n_out = 0;  alt_dev = 0;  diverged = false;  tt = 0;
    for i = 1:Ng
        tt = (i - 1) * dt_sim;
        [ref, gi] = gov.step(guam.state, tt);
        n_out = n_out + (gi.V_cur >= 0 && gi.V_next >= 0);
        alt_dev = max(alt_dev, abs(-guam.state(3) - 80));
        if isnan(t_fin) && gi.progress >= traj.n_trim, t_fin = tt; end
        guam.step(ref);
        if abs(-guam.state(3) - 80) > 500 || guam.state(4) > 300
            diverged = true;  break;
        end
        if ~isnan(t_fin) && tt > t_fin + 2, break; end
    end

    accel = (traj.UH(2) - traj.UH(1)) / TS;
    if diverged
        verdict = 'DIVERGED';
    elseif isnan(t_fin)
        verdict = 'timeout';
    else
        verdict = 'ok';
    end
    fprintf('%5.2f | %7.2f | %7.2f | %7.2f | %+6.2f | %6.2f | %8.2f | %7.0f | %s\n', ...
        TS, accel, traj.t_node(end), t_fin, t_fin - traj.t_node(end), ...
        gov.n_hold * dt_sim, n_out * dt_sim, alt_dev, verdict);

    R(ii).T_seg = TS;  R(ii).reached = t_fin;  R(ii).ok = strcmp(verdict, 'ok');
end

ok = [R.ok];
if any(ok)
    fastest = min([R(ok).T_seg]);
    best = R(ok);
    [~, jb] = min([best.reached]);
    fprintf(['\nfastest T_seg that still completes: %.2f s ' ...
             '(schedule accel %.2f ft/s^2)\n'], fastest, ...
            (8.439) / fastest);
    fprintf('shortest actual mission: %.2f s at T_seg = %.2f s\n', ...
            best(jb).reached, best(jb).T_seg);
    fprintf(['NOTE: a shorter T_seg does not always finish sooner - the ' ...
             'governor gives back\n      whatever the schedule over-demands.\n']);
end
