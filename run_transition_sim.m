% run_transition_sim - run the closed-loop transition sim (filter ON/OFF) and
% save the results to a .mat. No plotting; visualize separately from the file:
%   plot_sim_diagnostics | visualize_tube_overlay_trace | visualize_tube_timeslices
clear all; close all;
here = fileparts(mfilename('fullpath'));
addpath(genpath(here));
results_mat = fullfile(here, 'reachable_data', 'transition_results.mat');

modes   = {'blend', 'off'};
loggers = struct();
configs = struct();
for i = 1:numel(modes)
    [loggers.(modes{i}), configs.(modes{i})] = run_once(struct('filter_mode', modes{i}));
end
cfg = configs.blend;

results = struct('dt', cfg.sim.dt, ...
    'meta', struct('scenario', cfg.sim.scenario, ...
                   'target_vel', cfg.controller.target_vel, ...
                   'steps', cfg.sim.steps));
for i = 1:numel(modes)
    m = modes{i};
    results.(m).trace = loggers.(m).exportTrace();
    results.(m).data  = loggers.(m).exportData();
end
if ~isfolder(fileparts(results_mat)), mkdir(fileparts(results_mat)); end
save(results_mat, 'results');
fprintf('Saved results to %s\n', results_mat);

% Console summary (blend vs off; both runs share the same length)
blend = results.blend.trace;  off = results.off.trace;
active = blend.active == 1;
ft2m = 0.3048;
fprintf('\n=== Liveness filter effect (blend vs off) ===\n');
fprintf('filter active steps   : %d / %d (%.1f%%)\n', nnz(active), numel(active), 100*mean(active));
fprintf('mean cmd change|active: %.4g\n', mean(blend.cmdChange(active), 'omitnan'));
fprintf('blend inside-tube     : %.1f%% (V<=0)\n', 100*mean(blend.V <= 0, 'omitnan'));
fprintf('max |traj deviation|  : u %.3g m/s | w %.3g m/s | q %.3g deg/s | theta %.3g deg\n', ...
    ft2m*max(abs(blend.uBody    - off.uBody)), ...
    ft2m*max(abs(blend.w        - off.w)), ...
    rad2deg(max(abs(blend.q     - off.q))), ...
    max(abs(blend.thetaDeg      - off.thetaDeg)));

%% -------------------------------------------------------------------------
function [logger, cfg] = run_once(overrides)
    cfg        = Config([], overrides);
    controller = Controller(cfg.controller, cfg.sim.dt);
    guam       = LpC_GUAM(cfg);
    logger     = SimLogger(cfg.logger, cfg.sim);

    ref_traj = cfg.controller.getReferenceTrajectory();
    N = size(ref_traj.pos, 2);

    [x0, engine0, surface0] = controller.initial_condition(ref_traj);
    guam.reset(x0, engine0, surface0);
    controller.reset();

    for k = 1:N
        ref.pos = ref_traj.pos(:, k);  ref.vel = ref_traj.vel(:, k);
        ref.chi = ref_traj.chi(k);     ref.chi_dot = ref_traj.chidot(k);
        logger.addState(guam, k, ref);
        [engine_cmd, surface_cmd] = controller.control(guam.state, ref);
        guam.step(engine_cmd, surface_cmd);
        logger.addInputs(controller, guam.engine, guam.surface);
    end
    logger.finalize();
end
