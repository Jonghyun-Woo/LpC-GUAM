% run_transition_sim - Flat-earth m-code port of the GUAM hover-to-cruise
% transition demo, and the single entry point for both basic runs and
% liveness-filter (BRT) verification. Common climb/accel profile:
%   0-20 s : vertical climb 0 -> 80 ft (initial climb rate 8 ft/s)
%   20-40 s: accelerate to 15 ft/s forward flight
% Cruise altitude depends on the scenario:
%   'althold' (default) : hold 80 ft through cruise
%   'climb'             : keep climbing to 100 ft (original)
%   'lon_brt_verify'    : descending WH3 verification profile (filter study)
%
% Pipeline: assemble a central Config -> LpC_GUAM(cfg) -> run the explicit
% closed-loop step loop (local run_once) -> SimLogger buffers, plots, saves.
% When the liveness filter is enabled (params.filter_mode ~= 'off'), the
% script additionally runs a filter-OFF companion and overlays the filtered
% vs. nominal trajectory on the LON BRT corridor (replacing diag_e3_*).
% Run this script from anywhere; it adds the repo to the path.
clear all; close all;
here = fileparts(mfilename('fullpath'));
addpath(genpath(here));

%% Overriding simulation parameters (optional)
% Pre-set `scenario` / `params` before running to override these defaults.
if ~exist('scenario', 'var') || isempty(scenario)
    scenario = 'lon_brt_verify';   % 'althold' | 'climb' | 'lon_brt_verify'
end
if ~exist('params', 'var') || isempty(params)
    params = struct();
    params.steps = 10000;
    params.target_vel = 150;

    params.filter_mode = 'blend';   % 'off' | 'blend' | 'lr'
    params.filter_wh_anchor = [];

    % Logger options (optional):
    % params.saveFigures = true;  % save figures as PNG
    % params.saveDir     = '';    % '' -> logs/<timestamp>
end

%% 1) Run both filter modes in one loop: primary (requested) + OFF baseline.
% modes = unique({params.filter_mode, 'off'}, 'stable');
modes = {'blend', 'off'}; % least restrictive mode is added soon.
run_loggers = struct();
run_cfgs    = struct();
for i = 1:length(modes)
    m = modes{i};
    params.filter_mode = m;
    [run_loggers.(m), run_cfgs.(m)] = run_once(scenario, params);
end


%% 2) Single-run plots from the filter-ON (primary) run only
logger  = run_loggers.('blend');
cfg     = run_cfgs.('blend');
logger.plot();

%% 3) Filter verification: ON vs OFF trajectory overlay on the LON BRT corridor
if isfield(run_loggers, 'off')
    tr_on  = run_loggers.blend.exportTrace();      % filtered
    tr_off = run_loggers.off.exportTrace();        % nominal

    % Plot options
    opts = struct();
    opts.tube        = 'brt';
    opts.wh_idx      = 3;
    opts.uh_list     = 1:20;
    opts.coordMode   = 'absolute';
    opts.shiftByTrim = true;
    opts.plot_q_theta_2d = false;
    opts.compareOnlyUW   = false;
    opts.mainLabel   = 'Filtered trajectory';

    opts.extraTraces = struct();
    opts.extraTraces(1).trace     = tr_off;
    opts.extraTraces(1).label     = 'Nominal RSLQR trajectory';
    opts.extraTraces(1).lineStyle = '--';
    opts.extraTraces(1).lineWidth = 2.0;
    opts.extraTraces(1).color     = [0.85 0.10 0.10];

    figs = plot_lon_overlay_pair(tr_on, opts);

    if cfg.logger.saveFigures
        logger.saveFigure(figs.uw,        'lon_overlay_uw');
        logger.saveFigure(figs.u_q_theta, 'lon_overlay_u_q_theta');
    end
end

% -------------------------------------------------------------------------
function [L, cfg] = run_once(scenario, params)
% One closed-loop transition run. The step loop stays inline here
% SimLogger only buffers/plots. Returns the finalized logger
% and the config hub used.
cfg  = Config(scenario, params);
guam = LpC_GUAM(cfg);
L    = SimLogger(cfg.logger, cfg.sim);

rt = guam.refTraj;  N = size(rt.pos, 2);
guam.reset();
for k = 1:N
    % Update Reference (TODO: Move this into LpC_GUAM.step() to avoid exposing the refTraj and improve this logic)
    ref.pos = rt.pos(:, k);  ref.vel = rt.vel(:, k);
    ref.chi = rt.chi(k);     ref.chi_dot = rt.chidot(k);

    % Record the per-step state
    L.addState(guam, k, ref);
    % Physics step
    [engine, surface] = guam.step(ref);
    % Record the post-step control outputs including the safety filter states.
    L.addInputs(guam, engine, surface);
end
L.finalize();
end
