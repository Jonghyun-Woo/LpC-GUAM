% sweep_linear_margin - can a constraint tightening buy back the safety the
% linear predictor lost?
%
% With margin = 0 the linear-prediction governor ran 77x faster but left the
% tube for 7.0 % of the mission (worst V +0.126). Tightening the admissible
% set to V <= -margin should pay for the prediction error; the question is
% how much mission time that costs.
%
% Reference (nonlinear rollout): 51.70 s, 0.0 % outside, worst V -0.000,
% altitude deviation 29.2 ft, wall clock 1285 s.
clear; close all;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

dt = 0.01;  N_dec = 10;  SIM_MAX = 120;
params = struct('filter_mode','off');  params.T_seg = 2.0;

traj = TrimScheduleTrajectory.build(ReferenceTrajectory.trimOpts(dt, params));
brtV = BRTValue(fullfile(root,'data'), traj.trim_lon, traj.wh_idx);
S = load(fullfile(root,'logger','closed_loop_model.mat'), 'M');

MARGINS = [0 0.05 0.10 0.15 0.20 0.30];
fprintf('\n margin | flown  | V>=0 [%%] | worst V | alt dev | wall  | verdict\n');
fprintf('%s\n', repmat('-', 1, 72));

for mg = MARGINS
    guam = LpC_GUAM(Config('trim_schedule', params));
    gov  = TrimRefGovernor(traj, brtV, guam, dt, ...
            struct('T',6.0,'delta',0.3,'eps',0.02,'M',40, ...
                   'lin',S.M,'margin',mg));
    guam.reset();  gov.reset();

    v_final = traj.trim_lon(1,end);  t_done = NaN;
    Ng = round(SIM_MAX/dt);  Vw = -inf;  ad = 0;  nviol = 0;  n = 0;
    tic;
    for i = 1:Ng
        tt = (i-1)*dt;
        if mod(i-1,N_dec)==0, ref = gov.step(tt); else, ref = gov.hold(); end
        Vn = brtV.value(guam.state([4 6 11 8]), gov.v);
        Vw = max(Vw, Vn);  nviol = nviol + (Vn >= 0);  n = n + 1;
        ad = max(ad, abs(-guam.state(3) - 80));
        if isnan(t_done) && gov.v >= v_final-1e-6, t_done = tt; end
        guam.step(ref);
        if ~isnan(t_done) && tt > t_done+2, break; end
    end
    w = toc;

    if isnan(t_done)
        fprintf('%7.2f |   --   | %8.1f | %+7.3f | %7.1f | %5.1f | TIMEOUT\n', ...
                mg, 100*nviol/n, Vw, ad, w);
    else
        ok = 'ok';
        if Vw >= 0, ok = 'UNSAFE'; end
        fprintf('%7.2f | %6.2f | %8.1f | %+7.3f | %7.1f | %5.1f | %s\n', ...
                mg, t_done, 100*nviol/n, Vw, ad, w, ok);
    end
end

fprintf(['\nnonlinear reference: 51.70 s, 0.0 %%, -0.000, 29.2 ft, ' ...
         '1285 s wall\n']);
