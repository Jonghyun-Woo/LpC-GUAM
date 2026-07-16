function figs = plot_lon_overlay_pair(R_or_trace, opts)
% plot_lon_overlay_pair
%
% Draw:
%   1) u-w 2D overlay
%   2) u-q-theta 3D overlay
% Optional:
%   3) q-theta 2D overlay

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    opts = fill_pair_defaults(opts);

    figs = struct();

    % ---------------------------------------------------------------------
    % 1) u-w overlay
    % ---------------------------------------------------------------------
    opts_uw = opts;
    opts_uw.keep_dims = [1 2];
    opts_uw.figTitle = sprintf('Full GUAM trajectory on LON %s corridor: u vs w', ...
        upper(opts.tube));
    opts_uw.filePrefix = [opts.filePrefix '_uw'];

    figs.uw = visualize_lon_tube_overlay_trace(R_or_trace, opts_uw);

    % ---------------------------------------------------------------------
    % 2) u-q-theta overlay
    % ---------------------------------------------------------------------
    opts_uqth = opts;
    opts_uqth.keep_dims = [1 3 4];
    opts_uqth.figTitle = sprintf('Full GUAM trajectory on LON %s corridor: u vs q vs theta', ...
        upper(opts.tube));
    opts_uqth.filePrefix = [opts.filePrefix '_u_q_theta'];

    % 첫 번째 u-w plot에만 nominal trajectory를 비교하고 싶으면,
    % 두 번째 plot에서는 extra trajectory를 제거.
    if opts.compareOnlyUW
        opts_uqth.extraTraces = [];
    end

    figs.u_q_theta = visualize_lon_tube_overlay_trace(R_or_trace, opts_uqth);

    % ---------------------------------------------------------------------
    % 3) Optional q-theta 2D overlay
    % ---------------------------------------------------------------------
    if opts.plot_q_theta_2d
        opts_qth = opts;
        opts_qth.keep_dims = [3 4];
        opts_qth.figTitle = sprintf('Full GUAM trajectory on LON %s corridor: q vs theta', ...
            upper(opts.tube));
        opts_qth.filePrefix = [opts.filePrefix '_q_theta'];

        if opts.compareOnlyUW
            opts_qth.extraTraces = [];
        end

        figs.q_theta = visualize_lon_tube_overlay_trace(R_or_trace, opts_qth);
    end
end

% -------------------------------------------------------------------------
function opts = fill_pair_defaults(opts)
    if ~isfield(opts, 'tube'),        opts.tube = 'brt'; end
    if ~isfield(opts, 'dataDir'),     opts.dataDir = 'tables'; end
    if ~isfield(opts, 'wh_idx'),      opts.wh_idx = 3; end
    if ~isfield(opts, 'uh_list'),     opts.uh_list = []; end
    if ~isfield(opts, 'coordMode'),   opts.coordMode = 'absolute'; end
    if ~isfield(opts, 'shiftByTrim'), opts.shiftByTrim = true; end
    if ~isfield(opts, 'filePrefix'),  opts.filePrefix = 'full_guam_overlay'; end
    if ~isfield(opts, 'plot_q_theta_2d'), opts.plot_q_theta_2d = false; end
    if ~isfield(opts, 'mainLabel'), opts.mainLabel = 'Full GUAM trajectory'; end
    if ~isfield(opts, 'extraTraces'), opts.extraTraces = []; end
    if ~isfield(opts, 'compareOnlyUW'), opts.compareOnlyUW = true; end
end