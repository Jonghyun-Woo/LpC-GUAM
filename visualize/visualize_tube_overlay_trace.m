function fig = visualize_tube_overlay_trace(trace, opts)
    % Overlay the BRT V=0 corridor (every UH tube at opts.wh_idx, each sliced at
    % its trim centre) with the flown trajectory and the nominal filter-off run.
    % opts.plane_dims picks the two in-plane dims; opts.stack=true places each
    % UH tube on its trim-airspeed plane for a 3D corridor-vs-airspeed view.
    %   LON: u-w (stack off), u-q-theta (stack on)
    %   LAT: u-v-p, u-r-phi (stack on)
    % trace: a trace struct, a results .mat path, or empty to load the default.
    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = fill_defaults(opts);

    nominal = [];
    if nargin < 1 || isempty(trace) || ischar(trace) || isstring(trace)
        if nargin < 1, trace = ''; end
        r       = load_sim_results(char(trace));
        trace   = r.blend.trace;
        nominal = r.off.trace;
    end

    P      = axis_plot_config(opts.axis);
    dim_x  = opts.plane_dims(1);
    dim_y  = opts.plane_dims(2);
    is3d   = opts.stack;
    ft2m   = 0.3048;
    S      = load('trim_table_Poly_ConcatVer4p0.mat');
    axU    = upper(opts.axis);
    tube   = upper(opts.tube);
    prefix = sprintf('reachability_data/guam_output/%s_%s/GUAM_%s_%s', axU, tube, axU, tube);
    tube_color = 'b';  if strcmp(tube, 'FRT'), tube_color = 'r'; end

    st      = trace_states_native(trace, opts.axis);
    valid   = all(isfinite(st), 2);
    plane   = st(valid, :) .* P.scale(:)';
    airspd  = trace.uBody(valid) * ft2m;
    k_valid = trace.k(valid);
    traj    = pack_coords(is3d, airspd, plane(:, dim_x), plane(:, dim_y));
    brt_exit  = first_valid_marker(trace.brtExitFlag,  valid);
    grid_exit = first_valid_marker(trace.gridExitFlag, valid);

    fig = figure; clf; set(gcf, 'Color', 'w', 'Position', [100 100 1200 800]);
    hold on; grid on; box on;

    h_tube = []; h_target = [];
    for uh_idx = opts.uh_list
        fname = sprintf('%s_UH%d_WH%d.mat', prefix, uh_idx, opts.wh_idx);
        if ~isfile(fname), continue; end
        tube_mat   = load(fname);
        value_grid = double(tube_mat.values);
        grid_axes  = tube_mat.grid_axes;

        X0 = S.XU0_interp(1:12, uh_idx, opts.wh_idx);
        u_plane = X0(1) * ft2m;
        [x_grid, y_grid, Vslice] = plane_slice_center(value_grid, grid_axes, dim_x, dim_y);
        x_grid = (x_grid + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
        y_grid = (y_grid + X0(P.trim_rows(dim_y))) * P.scale(dim_y);

        if is3d
            h = draw_contour_on_plane(u_plane, x_grid, y_grid, Vslice, tube_color);
        else
            [~, h] = contour(x_grid, y_grid, Vslice', [0 0], 'Color', tube_color, 'LineWidth', 0.9);
        end
        if isempty(h_tube) && ~isempty(h), h_tube = h; end

        xlb = (-P.target_ub(dim_x) + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
        xub = ( P.target_ub(dim_x) + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
        ylb = (-P.target_ub(dim_y) + X0(P.trim_rows(dim_y))) * P.scale(dim_y);
        yub = ( P.target_ub(dim_y) + X0(P.trim_rows(dim_y))) * P.scale(dim_y);
        box_pts = pack_coords(is3d, u_plane * ones(5, 1), ...
                              [xlb xub xub xlb xlb]', [ylb ylb yub yub ylb]');
        h_target = plot_pts(box_pts, 'g-', 'LineWidth', 1.5);
    end

    h_traj  = plot_pts(traj,        'k-', 'LineWidth', 2.0);
    h_start = plot_pts(traj(1, :),   'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    h_end   = plot_pts(traj(end, :), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);

    handles = h_traj;  names = {'Full GUAM trajectory'};
    if ~isempty(h_tube),   handles = [h_tube, handles];   names = [{tube}, names];     end
    if ~isempty(h_target), handles = [handles, h_target]; names = [names, {'Target'}]; end

    if ~isempty(nominal)
        ns    = trace_states_native(nominal, opts.axis);
        nv    = all(isfinite(ns), 2);
        nplane = ns(nv, :) .* P.scale(:)';
        nom   = pack_coords(is3d, nominal.uBody(nv) * ft2m, nplane(:, dim_x), nplane(:, dim_y));
        h_nom = plot_pts(nom, '--', 'Color', [0.85 0.10 0.10], 'LineWidth', 2.0);
    else
        h_nom = [];
    end

    handles = [handles, h_start, h_end];
    names   = [names, {sprintf('Start, k=%d', k_valid(1)), sprintf('End, k=%d', k_valid(end))}];
    if ~isempty(h_nom)
        handles = [handles, h_nom];  names = [names, {'Nominal RSLQR trajectory'}];
    end
    if ~isempty(brt_exit)
        h = plot_pts(traj(brt_exit, :), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
        handles = [handles, h];  names = [names, {sprintf('BRT exit, k=%d', k_valid(brt_exit))}];
    end
    if ~isempty(grid_exit)
        h = plot_pts(traj(grid_exit, :), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 8);
        handles = [handles, h];  names = [names, {sprintf('Grid exit, k=%d', k_valid(grid_exit))}];
    end

    if is3d
        xlabel('u (m/s)');
        ylabel(P.labels{dim_x});
        zlabel(P.labels{dim_y});
        view(-35, 22);
        pbaspect([5, 1, 1]);
    else
        xlabel(P.labels{dim_x}); ylabel(P.labels{dim_y});
    end
    legend(handles, names, 'Location', 'best');

    if opts.save_png
        outName = opts.png_file;
        if isempty(outName)
            if is3d
                dims = sprintf('u_%s_%s', P.short{dim_x}, P.short{dim_y});
            else
                dims = sprintf('%s_%s', P.short{dim_x}, P.short{dim_y});
            end
            outName = sprintf('reachability_data/tube_overlay_trace_%s_%s_%s.png', axU, tube, dims);
        end
        exportgraphics(fig, outName, 'Resolution', 150);
        fprintf('saved %s\n', outName);
    end
end

% -------------------------------------------------------------------------
function opts = fill_defaults(opts)
    opts.axis       = 'lon';
    opts.tube       = 'brt';
    opts.wh_idx     = 2;
    opts.uh_list    = 1:20;
    opts.plane_dims = [1, 2];
    opts.stack      = false;
    opts.save_png   = true;
    opts.png_file   = '';
end

function P = axis_plot_config(ax)
    spec = FilterConfig.axisSpec(ax);
    ft2m = 0.3048;  r2d = 180 / pi;
    switch lower(ax)
        case 'lon'
            P.scale  = [ft2m; ft2m; r2d; r2d];
            P.labels = {'u (m/s)', 'w (m/s)', 'q (deg/s)', '\theta (deg)'};
            P.short  = {'u', 'w', 'q', 'theta'};
        case 'lat'
            P.scale  = [ft2m; r2d; r2d; r2d];
            P.labels = {'v (m/s)', 'p (deg/s)', 'r (deg/s)', '\phi (deg)'};
            P.short  = {'v', 'p', 'r', 'phi'};
        otherwise
            error('visualize_tube_overlay_trace: unknown axis "%s".', ax);
    end
    P.trim_rows = spec.trim_rows(:);
    P.target_ub = spec.target_ub(:);
end

function state_native = trace_states_native(trace, ax)
    switch lower(ax)
        case 'lon'
            state_native = [trace.uBody(:), trace.w(:), trace.q(:), deg2rad(trace.thetaDeg(:))];
        case 'lat'
            state_native = [trace.v(:), trace.p(:), trace.r(:), deg2rad(trace.phiDeg(:))];
    end
end

function M = pack_coords(is3d, airspeed, a, b)
    % Trajectory/box coords: (airspeed, a, b) in 3D, (a, b) in 2D.
    if is3d, M = [airspeed(:), a(:), b(:)]; else, M = [a(:), b(:)]; end
end

function h = plot_pts(M, varargin)
    if size(M, 2) == 3
        h = plot3(M(:, 1), M(:, 2), M(:, 3), varargin{:});
    else
        h = plot(M(:, 1), M(:, 2), varargin{:});
    end
end

function [x_grid, y_grid, Vslice] = plane_slice_center(value_grid, grid_axes, dim_x, dim_y)
    % V over (dim_x, dim_y) with the complement dims held at the trim centre (0).
    x_grid = grid_axes{dim_x}(:)';
    y_grid = grid_axes{dim_y}(:)';
    [X, Y] = ndgrid(grid_axes{dim_x}, grid_axes{dim_y});
    Q = cell(4, 1);
    Q{dim_x} = X;  Q{dim_y} = Y;
    for d = setdiff(1:4, [dim_x dim_y])
        Q{d} = zeros(size(X));
    end
    Vslice = interpn(grid_axes{1}, grid_axes{2}, grid_axes{3}, grid_axes{4}, ...
                     value_grid, Q{1}, Q{2}, Q{3}, Q{4}, 'linear');
end

function h = draw_contour_on_plane(u_plane, x_grid, y_grid, Vslice, color)
    % Draw the V=0 contour of Vslice on the constant-airspeed plane u=u_plane.
    h = [];
    fin = Vslice(isfinite(Vslice));
    if isempty(fin) || min(fin) > 0 || max(fin) < 0, return; end
    C = contourc(x_grid, y_grid, Vslice', [0 0]);
    col = 1;
    while col < size(C, 2)
        n = C(2, col);
        if ~isfinite(n) || n < 2 || col + n > size(C, 2), break; end
        pts = C(:, col + 1 : col + n);
        hh = plot3(u_plane * ones(1, n), pts(1, :), pts(2, :), '-', 'Color', color, 'LineWidth', 0.9);
        if isempty(h), h = hh; end
        col = col + n + 1;
    end
end

function idx = first_valid_marker(flag, valid)
    idx = [];
    if isempty(flag), return; end
    abs_idx = find(flag(:) == 1 & valid(:), 1, 'first');
    if isempty(abs_idx), return; end
    idx = find(find(valid) == abs_idx, 1, 'first');
end
