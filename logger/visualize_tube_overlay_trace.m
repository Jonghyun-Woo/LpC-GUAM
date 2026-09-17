function fig = visualize_tube_overlay_trace(R_or_trace, opts)
    % visualize_tube_overlay_trace
    % opts fields:
    %   opts.axis        : 'lon' or 'lat' (selects grid/trim/labels via axisSpec)
    %   opts.tube        : 'brt' or 'frt'
    %   opts.wh_idx      : WH index, default 3
    %   opts.uh_list     : UH index list
    %   opts.keep_dims   : [1 2], [3 4], or [1 3 4] (schedule slice stack)
    %   opts.coordMode   : 'absolute' or 'anchor'
    %   opts.shiftByTrim : true/false
    %   opts.figTitle    : custom title
    %
    % Dimension convention (per axis, from FilterConfig.axisSpec):
    %   lon: 1:u 2:w 3:q 4:theta      lat: 1:v 2:p 3:r 4:phi

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    if isfield(R_or_trace, 'trace')
        trace = R_or_trace.trace;
    else
        trace = R_or_trace;
    end

    opts = fill_overlay_defaults(opts);

    S = load('trim_table_Poly_ConcatVer4p0.mat');

    if isempty(opts.uh_list)
        opts.uh_list = 1:numel(1:20); % UH = 1~20
    end

    P = get_axis_plot_config(opts.axis);

    axes_base = cell(4, 1);
    for d = 1:4
        axes_base{d} = linspace(P.grid_min_plot(d), P.grid_max_plot(d), P.grid_num(d));
    end

    keep_dims = opts.keep_dims(:)';
    if ~(numel(keep_dims) == 2 || numel(keep_dims) == 3)
        error('opts.keep_dims must have length 2 or 3.');
    end

    if numel(keep_dims) == 3 && ~isequal(keep_dims, [1 3 4])
        error(['3D visualization supports only keep_dims = [1 3 4] ', ...
            '(dim1-dim3-dim4 schedule slice stack).']);
    end

    tube = upper(opts.tube);
    if ~any(strcmp(tube, {'BRT', 'FRT'}))
        error('opts.tube must be ''brt'' or ''frt''.');
    end
    axU = upper(opts.axis);
    tubePrefix = sprintf('reachable_data/guam_output/%s_%s/GUAM_%s_%s', axU, tube, axU, tube);
    if strcmp(tube, 'BRT'), tubeColor = 'b'; else, tubeColor = 'r'; end
    tubeLegend = tube;

    traj4 = build_trace_trajectory(trace, opts.axis, opts.coordMode, P);
    valid = all(isfinite(traj4), 2);

    traj4_valid = traj4(valid, :);
    traj = traj4_valid(:, keep_dims);
    k_valid = trace.k(valid);

    extraTrajs = build_extra_trajectories(opts.extraTraces, keep_dims, opts.axis, opts.coordMode, P);

    brtExitLocalIdx  = first_marker_index(trace.brtExitFlag, valid);
    gridExitLocalIdx = first_marker_index(trace.gridExitFlag, valid);

    fig = figure;
    clf;

    try
        fig.Theme = 'light';
    catch
    end


    fig.Position = [100, 100, 1200, 800];

    hold on;

    h_tube_first = [];
    h_target_first = [];

    for ii = 1:numel(opts.uh_list)
        uh_idx = opts.uh_list(ii);

        filename = fullfile(sprintf('%s_UH%d_WH%d.mat', tubePrefix, uh_idx, opts.wh_idx));

        if ~isfile(filename)
            warning('File not found: %s. Skipping.', filename);
            continue;
        end

        tubeFile = load(filename);
        if isfield(tubeFile, 'values')
            data = double(tubeFile.values);
        elseif isfield(tubeFile, 'data')
            data = double(tubeFile.data);
        else
            warning('File %s has neither "values" nor "data". Skipping.', filename);
            continue;
        end

        [axis_plot, target_lb_plot, target_ub_plot] = ...
            build_plot_axes_for_schedule( ...
            S, uh_idx, opts.wh_idx, axes_base, ...
            P.target_lb_plot, P.target_ub_plot, P, opts);

        if numel(keep_dims) == 2
            [X, Y, Z] = slice_value_function_2d(data, axes_base, axis_plot, keep_dims);

            [~, h_tube] = contour(X, Y, Z, [0 0], ...
                'Color', tubeColor, 'LineWidth', 0.9);

            if isempty(h_tube_first) && ~isempty(h_tube)
                h_tube_first = h_tube;
            end

            h_target = draw_target_box_2d( ...
                target_lb_plot(keep_dims), target_ub_plot(keep_dims), 'g');

            if isempty(h_target_first) && ~isempty(h_target)
                h_target_first = h_target;
            end

        else
            % 3D dim1-dim3-dim4 stack: per-UH dim3-dim4 V=0 contour at local
            % dim1=dim2=0, placed at that schedule's absolute trim dim1.
            [d1Plane, d3Plot, d4Plot, Vslice] = ...
                slice_dims34_at_center( ...
                data, axes_base, axis_plot);

            h_tube = draw_zero_contour_on_d1_plane( ...
                d1Plane, d3Plot, d4Plot, Vslice, ...
                tubeColor, 0.9);

            h_target = draw_target_rect_on_d1_plane( ...
                d1Plane, ...
                target_lb_plot([3 4]), ...
                target_ub_plot([3 4]), ...
                'g');

            if isempty(h_tube_first) && ~isempty(h_tube)
                h_tube_first = h_tube;
            end

            if isempty(h_target_first) && ~isempty(h_target)
                h_target_first = h_target;
            end
        end
    end

    % ---------------------------------------------------------------------
    % Overlay trajectory
    % ---------------------------------------------------------------------
    if numel(keep_dims) == 2
        h_traj = plot(traj(:,1), traj(:,2), 'k-', 'LineWidth', 2.0);

        h_extra = gobjects(0);
        extra_names = {};

        for ei = 1:numel(extraTrajs)
            h_extra(end+1) = plot(extraTrajs(ei).traj(:,1), ...
                extraTrajs(ei).traj(:,2), ...
                'LineStyle', extraTrajs(ei).lineStyle, ...
                'Color', extraTrajs(ei).color, ...
                'LineWidth', extraTrajs(ei).lineWidth);
            extra_names{end+1} = extraTrajs(ei).label;
        end

        h_start = plot(traj(1,1), traj(1,2), 'ko', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 7);

        h_end = plot(traj(end,1), traj(end,2), 'ks', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 7);

        h_brtExit = [];
        if ~isempty(brtExitLocalIdx)
            h_brtExit = plot(traj(brtExitLocalIdx,1), traj(brtExitLocalIdx,2), 'ro', ...
                'MarkerFaceColor', 'r', 'MarkerSize', 8);
        end

        h_gridExit = [];
        if ~isempty(gridExitLocalIdx)
            h_gridExit = plot(traj(gridExitLocalIdx,1), traj(gridExitLocalIdx,2), 'mo', ...
                'MarkerFaceColor', 'm', 'MarkerSize', 8);
        end

        xlabel(P.labels{keep_dims(1)});
        ylabel(P.labels{keep_dims(2)});
        axis square;

    else
        h_traj = plot3(traj(:,1), traj(:,2), traj(:,3), 'k-', 'LineWidth', 2.0);

        h_extra = gobjects(0);
        extra_names = {};

        for ei = 1:numel(extraTrajs)
            h_extra(end+1) = plot3(extraTrajs(ei).traj(:,1), ...
                extraTrajs(ei).traj(:,2), ...
                extraTrajs(ei).traj(:,3), ...
                'LineStyle', extraTrajs(ei).lineStyle, ...
                'Color', extraTrajs(ei).color, ...
                'LineWidth', extraTrajs(ei).lineWidth);
            extra_names{end+1} = extraTrajs(ei).label;
        end

        h_start = plot3(traj(1,1), traj(1,2), traj(1,3), 'ko', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 7);

        h_end = plot3(traj(end,1), traj(end,2), traj(end,3), 'ks', ...
            'MarkerFaceColor', 'y', 'MarkerSize', 7);

        h_brtExit = [];
        if ~isempty(brtExitLocalIdx)
            h_brtExit = plot3(traj(brtExitLocalIdx,1), ...
                traj(brtExitLocalIdx,2), ...
                traj(brtExitLocalIdx,3), ...
                'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        end

        h_gridExit = [];
        if ~isempty(gridExitLocalIdx)
            h_gridExit = plot3(traj(gridExitLocalIdx,1), ...
                traj(gridExitLocalIdx,2), ...
                traj(gridExitLocalIdx,3), ...
                'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
        end

        xlabel(P.labels{keep_dims(1)});
        ylabel(P.labels{keep_dims(2)});
        zlabel(P.labels{keep_dims(3)});

        view(30, 24);
        axis tight;
        pbaspect([1.6 1.0 1.0]);
    end

    if isempty(opts.figTitle)
        opts.figTitle = make_default_title(opts, P.short_labels, keep_dims);
    end

    title(opts.figTitle);
    grid on;

    legend_handles = gobjects(0);
    legend_names = {};

    if ~isempty(h_tube_first)
        legend_handles(end+1) = h_tube_first;
        legend_names{end+1} = tubeLegend;
    end

    if ~isempty(h_target_first)
        legend_handles(end+1) = h_target_first;
        legend_names{end+1} = 'Target';
    end

    legend_handles(end+1) = h_traj;
    legend_names{end+1} = opts.mainLabel;

    for ei = 1:numel(h_extra)
        if isgraphics(h_extra(ei))
            legend_handles(end+1) = h_extra(ei);
            legend_names{end+1} = extra_names{ei};
        end
    end

    legend_handles(end+1) = h_start;
    legend_names{end+1} = sprintf('Start, k=%d', k_valid(1));

    legend_handles(end+1) = h_end;
    legend_names{end+1} = sprintf('End, k=%d', k_valid(end));

    if ~isempty(h_brtExit)
        legend_handles(end+1) = h_brtExit;
        legend_names{end+1} = sprintf('BRT exit, k=%d', k_valid(brtExitLocalIdx));
    end

    if ~isempty(h_gridExit)
        legend_handles(end+1) = h_gridExit;
        legend_names{end+1} = sprintf('Grid exit, k=%d', k_valid(gridExitLocalIdx));
    end

    legend(legend_handles, legend_names, 'Location', 'best');
    set(gcf, 'Color', 'w');
end

% =========================================================================
% Helper functions
% =========================================================================

function opts = fill_overlay_defaults(opts)
    if ~isfield(opts, 'axis'),        opts.axis = 'lon';                        end
    if ~isfield(opts, 'tube'),        opts.tube = 'brt';                        end
    if ~isfield(opts, 'wh_idx'),      opts.wh_idx = 3;                          end
    if ~isfield(opts, 'uh_list'),     opts.uh_list = [];                        end
    if ~isfield(opts, 'keep_dims'),   opts.keep_dims = [1 2];                   end
    if ~isfield(opts, 'coordMode'),   opts.coordMode = 'absolute';              end
    if ~isfield(opts, 'shiftByTrim'), opts.shiftByTrim = true;                  end
    if ~isfield(opts, 'figTitle'),    opts.figTitle = '';                       end
    if ~isfield(opts, 'filePrefix'),  opts.filePrefix = 'overlay';              end
    if ~isfield(opts, 'mainLabel'),   opts.mainLabel = 'Full GUAM trajectory';  end
    if ~isfield(opts, 'extraTraces'), opts.extraTraces = [];                    end
end

function P = get_axis_plot_config(axis)
    % Grid/target/label config for the requested axis, driven by
    % FilterConfig.axisSpec so the plot grid always matches the BRT tables.
    spec = FilterConfig.axisSpec(axis);

    P = struct();
    P.ft2m    = 0.3048;
    P.rad2deg = 180 / pi;

    switch lower(axis)
        case 'lon'
            P.plot_scale  = [P.ft2m; P.ft2m; P.rad2deg; P.rad2deg];
            P.labels      = {'u (m/s)', 'w (m/s)', 'q (deg/s)', '\theta (deg)'};
            P.short_labels = {'u', 'w', 'q', 'theta'};
        case 'lat'
            P.plot_scale  = [P.ft2m; P.rad2deg; P.rad2deg; P.rad2deg];
            P.labels      = {'v (m/s)', 'p (deg/s)', 'r (deg/s)', '\phi (deg)'};
            P.short_labels = {'v', 'p', 'r', 'phi'};
        otherwise
            error('get_axis_plot_config: unknown axis "%s".', axis);
    end

    P.trim_rows = spec.trim_rows(:);
    P.grid_min  = spec.grid_min(:);
    P.grid_max  = spec.grid_max(:);
    P.grid_num  = spec.grid_num(:)';
    P.target_ub = spec.target_ub(:);
    P.target_lb = -P.target_ub;

    P.grid_min_plot  = P.grid_min  .* P.plot_scale;
    P.grid_max_plot  = P.grid_max  .* P.plot_scale;
    P.target_lb_plot = P.target_lb .* P.plot_scale;
    P.target_ub_plot = P.target_ub .* P.plot_scale;
end

function traj4 = build_trace_trajectory(trace, axis, coordMode, P)
    switch lower(coordMode)
        case 'absolute'
            switch lower(axis)
                case 'lon'
                    raw = [trace.uBody(:), trace.w(:), ...
                           trace.q(:), deg2rad(trace.thetaDeg(:))];
                case 'lat'
                    raw = [trace.v(:), trace.p(:), ...
                           trace.r(:), deg2rad(trace.phiDeg(:))];
                otherwise
                    error('build_trace_trajectory: unknown axis "%s".', axis);
            end
            traj4 = raw .* P.plot_scale(:)';

        case 'anchor'
            traj4 = trace.xF(:, 1:4) .* P.plot_scale(:)';

        otherwise
            error('coordMode must be ''absolute'' or ''anchor''.');
    end
end

function [axis_plot, target_lb_plot, target_ub_plot] = ...
    build_plot_axes_for_schedule(S, uh_idx, wh_idx, axes_base, ...
    target_lb_base, target_ub_base, P, opts)

    axis_plot = axes_base;
    target_lb_plot = target_lb_base;
    target_ub_plot = target_ub_base;

    if strcmpi(opts.coordMode, 'absolute') && opts.shiftByTrim
        X0 = S.XU0_interp(1:12, uh_idx, wh_idx);

        trim4 = X0(P.trim_rows) .* P.plot_scale;

        for d = 1:4
            axis_plot{d} = axes_base{d} + trim4(d);
        end

        target_lb_plot = target_lb_base + trim4;
        target_ub_plot = target_ub_base + trim4;
    end
end

function [X, Y, Z] = slice_value_function_2d(data, axes_base, axis_plot, keep_dims)
    d1 = keep_dims(1);
    d2 = keep_dims(2);

    x1 = axes_base{d1};
    x2 = axes_base{d2};

    [Q1, Q2] = ndgrid(x1, x2);

    Q = cell(4,1);
    for d = 1:4
        Q{d} = zeros(size(Q1));
    end

    Q{d1} = Q1;
    Q{d2} = Q2;

    Zraw = interpn(axes_base{1}, axes_base{2}, axes_base{3}, axes_base{4}, ...
        data, Q{1}, Q{2}, Q{3}, Q{4}, 'linear');

    % contour expects Z as length(y)-by-length(x)
    xPlot = axis_plot{d1};
    yPlot = axis_plot{d2};

    [X, Y] = meshgrid(xPlot, yPlot);
    Z = Zraw';
end

function [d1Plane, d3Plot, d4Plot, Vslice] = ...
    slice_dims34_at_center(data, axes_base, axis_plot)

    d3Base = axes_base{3};
    d4Base = axes_base{4};

    [D3, D4] = ndgrid(d3Base, d4Base);

    D1z = zeros(size(D3));
    D2z = zeros(size(D3));

    Vslice = interpn( ...
        axes_base{1}, ...
        axes_base{2}, ...
        axes_base{3}, ...
        axes_base{4}, ...
        data, ...
        D1z, D2z, D3, D4, ...
        'linear');

    d3Plot = axis_plot{3};
    d4Plot = axis_plot{4};

    d1Plane = interp1( ...
        axes_base{1}, ...
        axis_plot{1}, ...
        0, ...
        'linear', ...
        'extrap');
end

function hFirst = draw_zero_contour_on_d1_plane( ...
    d1Plane, d3Axis, d4Axis, Vslice, colorValue, lineWidth)
% Draw the dim3-dim4 slice V = 0 contour on the dim1 = d1Plane plane.

    hFirst = [];

    finiteValues = Vslice(isfinite(Vslice));
    if isempty(finiteValues)
        return;
    end

    vMin = min(finiteValues);
    vMax = max(finiteValues);

    % No V = 0 crossing inside the slice -> no contour to draw.
    if vMin > 0 || vMax < 0 || vMin == vMax
        return;
    end

    % contourc expects Z as length(y)-by-length(x).
    C = contourc(d3Axis, d4Axis, Vslice', [0 0]);

    col = 1;
    while col < size(C, 2)
        numPoints = C(2, col);

        if ~isfinite(numPoints) || numPoints < 2 || ...
                col + numPoints > size(C, 2)
            break;
        end

        points = C(:, col + 1 : col + numPoints);

        h = plot3( ...
            d1Plane * ones(1, numPoints), ...
            points(1, :), ...
            points(2, :), ...
            '-', ...
            'Color', colorValue, ...
            'LineWidth', lineWidth);

        if isempty(hFirst)
            hFirst = h;
        end

        col = col + numPoints + 1;
    end
end

function h = draw_target_box_2d(lb, ub, colorChar)
    x = [lb(1), ub(1), ub(1), lb(1), lb(1)];
    y = [lb(2), lb(2), ub(2), ub(2), lb(2)];

    h = plot(x, y, '-', 'Color', colorChar, 'LineWidth', 1.5);
end

function h = draw_target_rect_on_d1_plane( ...
    d1Plane, lb_34, ub_34, colorValue)

    d3 = [lb_34(1), ub_34(1), ub_34(1), lb_34(1), lb_34(1)];
    d4 = [lb_34(2), lb_34(2), ub_34(2), ub_34(2), lb_34(2)];

    d1 = d1Plane * ones(size(d3));

    h = plot3( ...
        d1, d3, d4, ...
        '-', ...
        'Color', colorValue, ...
        'LineWidth', 1.3);
end

function localIdx = first_marker_index(flag, valid)
    localIdx = [];

    if isempty(flag)
        return;
    end

    flag = flag(:);
    valid = valid(:);

    absIdx = find(flag == 1 & valid, 1, 'first');

    if isempty(absIdx)
        return;
    end

    validAbsIdx = find(valid);
    localIdx = find(validAbsIdx == absIdx, 1, 'first');
end

function titleStr = make_default_title(opts, short_labels, keep_dims)
    if numel(keep_dims) == 2
        titleStr = sprintf('%s %s overlay: %s vs %s, coord=%s', ...
            upper(opts.axis), ...
            upper(opts.tube), ...
            short_labels{keep_dims(1)}, ...
            short_labels{keep_dims(2)}, ...
            opts.coordMode);
    else
        titleStr = sprintf('%s %s overlay: %s, %s, %s, coord=%s', ...
            upper(opts.axis), ...
            upper(opts.tube), ...
            short_labels{keep_dims(1)}, ...
            short_labels{keep_dims(2)}, ...
            short_labels{keep_dims(3)}, ...
            opts.coordMode);
    end
end

function extraTrajs = build_extra_trajectories(extraSpecs, keep_dims, axis, coordMode, P)
    extraTrajs = struct([]);

    if isempty(extraSpecs)
        return;
    end

    for i = 1:numel(extraSpecs)
        spec = extraSpecs(i);

        if isfield(spec, 'trace')
            tr = spec.trace;
        elseif isfield(spec, 'R')
            tr = spec.R.trace;
        else
            warning('extraTraces(%d) has no trace or R field. Skipping.', i);
            continue;
        end

        traj4 = build_trace_trajectory(tr, axis, coordMode, P);
        valid = all(isfinite(traj4), 2);

        traj4 = traj4(valid, :);
        traj = traj4(:, keep_dims);

        if isempty(traj)
            continue;
        end

        extraTrajs(end+1).traj = traj; %#ok<AGROW>

        if isfield(spec, 'label')
            extraTrajs(end).label = spec.label;
        else
            extraTrajs(end).label = sprintf('extra trajectory %d', i);
        end

        if isfield(spec, 'lineStyle')
            extraTrajs(end).lineStyle = spec.lineStyle;
        else
            extraTrajs(end).lineStyle = '--';
        end

        if isfield(spec, 'lineWidth')
            extraTrajs(end).lineWidth = spec.lineWidth;
        else
            extraTrajs(end).lineWidth = 1.8;
        end

        if isfield(spec, 'color')
            extraTrajs(end).color = spec.color;
        else
            extraTrajs(end).color = [0.85 0.10 0.10];
        end
    end
end
