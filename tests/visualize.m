clear; close all;

addpath(genpath('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM\Refactoring'));
addpath(genpath('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM'));
cd('C:\Users\grape\OneDrive\CAU\AISL\LpC-GUAM\Refactoring');

% Filtered trajectory: ON
D_on = load('diag_plots/full_guam_caseA_currentU_linear/trace_q_0.1_ON_ref_150.0_dec_8.mat');
R_on = D_on.R;

% Nominal-controller trajectory: OFF
D_off = load('diag_plots/full_guam_caseA_currentU_linear/trace_q_0.1_OFF_ref_150.0_dec_8.mat');
R_off = D_off.R;

opts = struct();

opts.tube      = 'brt';
opts.wh_idx    = 3;
opts.uh_list   = 1:20;

opts.coordMode   = 'absolute';
opts.shiftByTrim = true;
opts.filePrefix = 'q_0_0_ON_vs_OFF';

opts.plot_q_theta_2d = false;

% Main trajectory label: ON / filtered
opts.mainLabel = 'Filtered trajectory';

% Additional trajectory: OFF / nominal
opts.extraTraces = struct();
opts.extraTraces(1).trace     = R_off.trace;
opts.extraTraces(1).label     = 'Nominal RSLQR trajectory';
opts.extraTraces(1).lineStyle = '--';
opts.extraTraces(1).lineWidth = 2.0;
opts.extraTraces(1).color     = [0.85 0.10 0.10];
opts.compareOnlyUW = false;

figs = plot_lon_overlay_pair(R_on, opts);