function fig = visualize_tube_timeslices(trace, opts)
    % Sample the trajectory every opts.interval_s seconds and, at each sample,
    % draw the V=0 contour of the BRT the vehicle is in (nearest UH at
    % opts.wh_idx), sliced at that sample's logged off-plane deviation, coloured
    % by time. opts.stack=true places each contour on its sample's body-u plane
    % for a 3D time-vs-state view (u-q-theta, u-v-p, u-r-phi); otherwise the
    % contours overlay on the opts.plane_dims plane. The full trajectory and a
    % time-coloured scatter of the sampled states are overlaid. trace: a trace
    % struct, a results .mat path, or empty to load the default run.
    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = fill_defaults(opts);

    if nargin < 1 || isempty(trace) || ischar(trace) || isstring(trace)
        if nargin < 1, trace = ''; end
        r       = load_sim_results(char(trace));
        opts.dt = r.dt;
        trace   = r.blend.trace;
    end

    P      = axis_plot_config(opts.axis);
    dim_x  = opts.plane_dims(1);
    dim_y  = opts.plane_dims(2);
    is3d   = opts.stack;
    ft2m   = 0.3048;
    S      = load('trim_table_Poly_ConcatVer4p0.mat');
    axU    = upper(opts.axis);
    tube   = upper(opts.tube);
    prefix = sprintf('reachable_data/guam_output/%s_%s/GUAM_%s_%s', axU, tube, axU, tube);

    state_native = trace_states_native(trace, opts.axis);   % N x 4, native units
    valid        = all(isfinite(state_native), 2);
    state_native = state_native(valid, :);
    traj_plot    = state_native .* P.scale(:)';
    t_valid      = (trace.k(valid) - 1) * opts.dt;
    u_body       = trace.uBody(valid);                       % ft/s
    u_axis       = u_body * ft2m;                            % m/s, 3D stacking axis
    t_max        = max(t_valid);

    t_marks    = 0 : opts.interval_s : t_max;
    sample_idx = arrayfun(@(tm) find(t_valid >= tm, 1, 'first'), t_marks, 'UniformOutput', false);
    sample_idx = unique([sample_idx{:}]);

    fig = figure; clf; set(gcf, 'Color', 'w'); hold on; grid on; box on;
    cmap = turbo(256);
    colormap(cmap);
    clim([0 t_max]);

    h_tube = [];
    for i = sample_idx
        uh_idx = nearest_uh(S.UH, u_body(i), opts.uh_list);
        fname  = sprintf('%s_UH%d_WH%d.mat', prefix, uh_idx, opts.wh_idx);
        if ~isfile(fname), continue; end
        tube_mat   = load(fname);
        value_grid = double(tube_mat.values);
        grid_axes  = tube_mat.grid_axes;

        % BRT grid is in trim-deviation coords: slice at the logged deviation.
        X0        = S.XU0_interp(1:12, uh_idx, opts.wh_idx);
        state_dev = state_native(i, :)' - X0(P.trim_rows);

        [x_grid, y_grid, Vslice] = plane_slice(value_grid, grid_axes, dim_x, dim_y, state_dev);
        x_grid = (x_grid + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
        y_grid = (y_grid + X0(P.trim_rows(dim_y))) * P.scale(dim_y);

        col = cmap(color_index(t_valid(i), t_max), :);
        if is3d
            h = draw_contour_on_plane(u_axis(i), x_grid, y_grid, Vslice, col);
        else
            [~, h] = contour(x_grid, y_grid, Vslice', [0 0], 'Color', col, 'LineWidth', 1.4);
        end
        if isempty(h_tube) && ~isempty(h), h_tube = h; end
    end

    traj    = pack_coords(is3d, u_axis, traj_plot(:, dim_x), traj_plot(:, dim_y));
    h_traj  = plot_pts(traj,        'k-', 'LineWidth', 1.5);
    h_state = scatter_pts(is3d, traj(sample_idx, :), t_valid(sample_idx));
    h_start = plot_pts(traj(1, :),   'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
    h_end   = plot_pts(traj(end, :), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);

    if is3d
        xlabel('u (m/s)'); ylabel(P.labels{dim_x}); zlabel(P.labels{dim_y});
        view([-35, 22]); pbaspect([5, 1, 1]);
    else
        xlabel(P.labels{dim_x}); ylabel(P.labels{dim_y});
    end
    title(sprintf('%s %s tube time-slices every %g s: %s vs %s', ...
        axU, tube, opts.interval_s, P.short{dim_x}, P.short{dim_y}));
    cb = colorbar; cb.Label.String = 'time (s)';

    handles = h_traj;
    names = {'trajectory'};

    handles = [handles, h_state, h_start, h_end];
    names   = [names, {'sampled state', 'start', 'end'}];
    legend(handles, names, 'Location', 'southeast');

    exportgraphics(gcf, sprintf('reachable_data/tube_timeslices_%s_%s_%s.png', ...
        axU, P.short{dim_x}, P.short{dim_y}), 'Resolution', 150);
end

%% -------------------------------------------------------------------------
function opts = fill_defaults(opts)
    opts.axis       = 'lat';
    opts.plane_dims = [1, 2];
    opts.interval_s = 1.5;
    opts.dt         = 0.01;
    opts.wh_idx     = 2;
    opts.uh_list    = 1:20;
    opts.tube       = 'brt';
    opts.stack      = false;
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
            error('visualize_tube_timeslices: unknown axis "%s".', ax);
    end
    P.trim_rows = spec.trim_rows(:);
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
    if is3d, M = [airspeed(:), a(:), b(:)]; else, M = [a(:), b(:)]; end
end

function h = plot_pts(M, varargin)
    if size(M, 2) == 3
        h = plot3(M(:, 1), M(:, 2), M(:, 3), varargin{:});
    else
        h = plot(M(:, 1), M(:, 2), varargin{:});
    end
end

function h = scatter_pts(is3d, M, c)
    if is3d
        h = scatter3(M(:, 1), M(:, 2), M(:, 3), 40, c, 'filled', 'MarkerEdgeColor', 'k');
    else
        h = scatter(M(:, 1), M(:, 2), 40, c, 'filled', 'MarkerEdgeColor', 'k');
    end
end

function uh_idx = nearest_uh(UH, u_now, uh_list)
    [~, j] = min(abs(UH(uh_list) - u_now));
    uh_idx = uh_list(j);
end

function [x_grid, y_grid, Vslice] = plane_slice(value_grid, grid_axes, dim_x, dim_y, state_dev)
    % V over (dim_x, dim_y); complement dims held at state_dev (clamped to grid).
    x_grid = grid_axes{dim_x}(:)';
    y_grid = grid_axes{dim_y}(:)';
    [X, Y] = ndgrid(grid_axes{dim_x}, grid_axes{dim_y});
    Q = cell(4, 1);
    Q{dim_x} = X;  Q{dim_y} = Y;
    for d = setdiff(1:4, [dim_x dim_y])
        v    = min(max(state_dev(d), grid_axes{d}(1)), grid_axes{d}(end));
        Q{d} = v * ones(size(X));
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
        hh = plot3(u_plane * ones(1, n), pts(1, :), pts(2, :), '-', 'Color', color, 'LineWidth', 1.4);
        if isempty(h), h = hh; end
        col = col + n + 1;
    end
end

function idx = color_index(t, t_max)
    if t_max <= 0, idx = 1; return; end
    idx = max(1, min(256, round(1 + 255 * t / t_max)));
end
