% run_transition_sim - Flat-earth m-code port of the GUAM hover-to-cruise
% transition demo; single entry point for basic runs and liveness-filter (BRT)
% verification. Scenarios: 'althold' (hold alt), 'climb' (climb 20 ft above alt),
% 'brt_verify' (level forward transition, BRT filter study).
% Pipeline: Config -> LpC_GUAM(cfg) -> closed-loop step loop (run_once) -> SimLogger.
% All simulation parameters are hardcoded in config/ (SimConfig, ControllerConfig,
% FilterConfig); this script only overrides the filter mode to build the ON/OFF pair.
clear all; close all;
here = fileparts(mfilename('fullpath'));
addpath(genpath(here));

%% 1) Run the filter-ON (blend) and filter-OFF baseline, keyed by mode
modes   = {'blend', 'off'};
loggers = struct();
configs = struct();
for i = 1:length(modes)
    mode = modes{i};
    [loggers.(mode), configs.(mode)] = run_once(struct('filter_mode', mode));
end

%% 2) Single-run plots from the filter-ON run
logger = loggers.blend;
cfg    = configs.blend;
logger.plot();

%% 3) Filter verification: ON vs OFF trajectory overlay on the LON BRT corridor
if isfield(loggers, 'off')
    trace_filtered = loggers.blend.exportTrace();
    trace_nominal  = loggers.off.exportTrace();

    opts = struct();
    opts.axis               = 'lon';
    opts.tube               = 'brt';
    opts.wh_idx             = 3;
    opts.uh_list            = 1:20;
    opts.coordMode          = 'absolute';
    opts.shiftByTrim        = true;
    opts.plot_q_theta_2d    = false;
    opts.compareOnlyUW      = false;
    opts.mainLabel          = 'Filtered trajectory';

    opts.extraTraces = struct();
    opts.extraTraces(1).trace     = trace_nominal;
    opts.extraTraces(1).label     = 'Nominal RSLQR trajectory';
    opts.extraTraces(1).lineStyle = '--';
    opts.extraTraces(1).lineWidth = 2.0;
    opts.extraTraces(1).color     = [0.85 0.10 0.10];

    % u-w 2D + lon airspeed-schedule 3D corridor (u-q-theta).
    figs = plot_tube_overlay_pair(trace_filtered, opts);

    % Lateral airspeed-schedule 3D corridors: u-v-p and u-phi-r.
    opts_lat = opts;
    opts_lat.axis = 'lat';
    opts_lat.plane_dims = [1 2];
    opts_lat.figTitle = 'LAT BRT corridor vs airspeed: u, v, p';
    fig_uvp = visualize_tube_overlay_trace(trace_filtered, opts_lat);

    opts_lat.plane_dims = [4 3];
    opts_lat.figTitle = 'LAT BRT corridor vs airspeed: u, phi, r';
    fig_uphir = visualize_tube_overlay_trace(trace_filtered, opts_lat);

    if cfg.logger.saveFigures
        logger.saveFigure(figs.uw,        'lon_overlay_uw');
        logger.saveFigure(figs.u_q_theta, 'lon_overlay_u_q_theta');
        logger.saveFigure(fig_uvp,        'lat_overlay_u_v_p');
        logger.saveFigure(fig_uphir,      'lat_overlay_u_phi_r');
    end
end

%% 4) Liveness-filter effect summary (blend vs off)
if isfield(loggers, 'off')
    trace_filtered = loggers.blend.exportTrace();
    trace_nominal  = loggers.off.exportTrace();

    n_common    = min(numel(trace_filtered.k), numel(trace_nominal.k));
    active_mask = trace_filtered.active(1:n_common) == 1;
    dev_u     = trace_filtered.uBody(1:n_common)    - trace_nominal.uBody(1:n_common);
    dev_w     = trace_filtered.w(1:n_common)        - trace_nominal.w(1:n_common);
    dev_q     = trace_filtered.q(1:n_common)        - trace_nominal.q(1:n_common);
    dev_theta = trace_filtered.thetaDeg(1:n_common) - trace_nominal.thetaDeg(1:n_common);

    fprintf('\n=== Liveness filter effect (blend vs off) ===\n');
    fprintf('filter active steps   : %d / %d (%.1f%%)\n', nnz(active_mask), n_common, 100*nnz(active_mask)/n_common);
    fprintf('mean cmd change|active: %.4g\n', mean(trace_filtered.cmdChange(active_mask), 'omitnan'));
    fprintf('blend inside-tube     : %.1f%% (V<=0)\n', 100*mean(trace_filtered.V(1:n_common) <= 0, 'omitnan'));
    ft2m = 0.3048;   % dev_u/dev_w are body velocities in ft/s; report in m/s
    fprintf('max |traj deviation|  : u %.3g m/s | w %.3g m/s | q %.3g deg/s | theta %.3g deg\n', ...
            ft2m*max(abs(dev_u)), ft2m*max(abs(dev_w)), rad2deg(max(abs(dev_q))), max(abs(dev_theta)));
end

% -------------------------------------------------------------------------
function [logger, cfg] = run_once(overrides)
    % One closed-loop transition run: the controller produces effector commands and
    % the plant (LpC_GUAM) integrates under them. Returns the finalized logger + config.
    % Scenario and all parameters come from config/ defaults; `overrides` carries
    % only the per-run filter mode.
    cfg        = Config([], overrides);
    controller = Controller(cfg.controller, cfg.sim.dt);
    guam       = LpC_GUAM(cfg);
    logger     = SimLogger(cfg.logger, cfg.sim);

    ref_traj       = cfg.controller.getReferenceTrajectory();
    num_ref_points = size(ref_traj.pos, 2);

    % Initialize the plant at the trim of the first reference point.
    [x0, engine0, surface0] = controller.initial_condition(ref_traj);
    guam.reset(x0, engine0, surface0);
    controller.reset();

    for k = 1:num_ref_points
        ref.pos = ref_traj.pos(:, k);  ref.vel = ref_traj.vel(:, k);
        ref.chi = ref_traj.chi(k);     ref.chi_dot = ref_traj.chidot(k);

        logger.addState(guam, k, ref);                              % pre-step state
        [engine_cmd, surface_cmd] = controller.control(guam.state, ref);
        guam.step(engine_cmd, surface_cmd);
        logger.addInputs(controller, guam.engine, guam.surface);    % post-step effectors + filter info
    end
    logger.finalize();
end
