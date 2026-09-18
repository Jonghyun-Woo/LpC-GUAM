function figs = plot_tube_overlay_pair(R_or_trace, opts)
    % plot_tube_overlay_pair
    %
    % Draw, for opts.axis ('lon' or 'lat'):
    %   1) dim1-dim2 2D overlay               (lon: u-w,   lat: v-p)
    %   2) airspeed schedule 3D corridor      (lon: u-q-theta, lat: u-r-phi)
    % Optional:
    %   3) dim3-dim4 2D overlay               (lon: q-theta, lat: r-phi)

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    opts = fill_pair_defaults(opts);
    axU = upper(opts.axis);

    figs = struct();

    % ---------------------------------------------------------------------
    % 1) dim1-dim2 overlay
    % ---------------------------------------------------------------------
    opts_12 = opts;
    opts_12.keep_dims = [1 2];
    opts_12.plane_dims = [];
    opts_12.figTitle = sprintf('Full GUAM trajectory on %s %s corridor: dim1 vs dim2', ...
        axU, upper(opts.tube));
    opts_12.filePrefix = [opts.filePrefix '_d1d2'];

    figs.uw = visualize_tube_overlay_trace(R_or_trace, opts_12);

    % ---------------------------------------------------------------------
    % 2) airspeed schedule 3D corridor (dim3-dim4 stacked over airspeed u)
    % ---------------------------------------------------------------------
    opts_3d = opts;
    opts_3d.plane_dims = [3 4];
    opts_3d.figTitle = sprintf('%s %s corridor vs airspeed: u, dim3, dim4', ...
        axU, upper(opts.tube));
    opts_3d.filePrefix = [opts.filePrefix '_u_d3d4'];

    % compareOnlyUW: overlay the nominal trajectory on the dim1-dim2 plot only.
    if opts.compareOnlyUW
        opts_3d.extraTraces = [];
    end

    figs.u_q_theta = visualize_tube_overlay_trace(R_or_trace, opts_3d);

    % ---------------------------------------------------------------------
    % 3) Optional dim3-dim4 2D overlay
    % ---------------------------------------------------------------------
    if opts.plot_q_theta_2d
        opts_34 = opts;
        opts_34.keep_dims = [3 4];
        opts_34.plane_dims = [];
        opts_34.figTitle = sprintf('Full GUAM trajectory on %s %s corridor: dim3 vs dim4', ...
            axU, upper(opts.tube));
        opts_34.filePrefix = [opts.filePrefix '_d3d4'];

        if opts.compareOnlyUW
            opts_34.extraTraces = [];
        end

        figs.q_theta = visualize_tube_overlay_trace(R_or_trace, opts_34);
    end
end

% -------------------------------------------------------------------------
function opts = fill_pair_defaults(opts)
    if ~isfield(opts, 'axis'),              opts.axis = 'lon';                          end
    if ~isfield(opts, 'tube'),              opts.tube = 'brt';                          end
    if ~isfield(opts, 'wh_idx'),            opts.wh_idx = 3;                            end
    if ~isfield(opts, 'uh_list'),           opts.uh_list = [];                          end
    if ~isfield(opts, 'coordMode'),         opts.coordMode = 'absolute';                end
    if ~isfield(opts, 'shiftByTrim'),       opts.shiftByTrim = true;                    end
    if ~isfield(opts, 'filePrefix'),        opts.filePrefix = 'full_guam_overlay';      end
    if ~isfield(opts, 'plot_q_theta_2d'),   opts.plot_q_theta_2d = false;               end
    if ~isfield(opts, 'mainLabel'),         opts.mainLabel = 'Full GUAM trajectory';    end
    if ~isfield(opts, 'extraTraces'),       opts.extraTraces = [];                      end
    if ~isfield(opts, 'compareOnlyUW'),     opts.compareOnlyUW = true;                  end
end
