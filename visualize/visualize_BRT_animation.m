function visualize_BRT_animation(trace, opts)
    % Animate how each plane's BRT V=0 contour deforms as the off-plane states
    % and the scheduled anchor trim evolve along a run, with the state
    % trajectory drawn as a growing trail + moving marker (as in
    % visualize_sim_animation). The drawn tube is the filter's selected curr/next
    % anchor (replicated from LivenessFilter), so the marker lies inside the
    % contour exactly when the logged V(t) <= 0. One GIF per view in opts.views:
    % longitudinal u-w and q-theta, lateral v-p and r-phi. trace: a trace
    % struct, a results .mat path, or empty to load the default run.
    if nargin < 2 || isempty(opts), opts = struct(); end
    opts = fill_defaults(opts);

    if nargin < 1 || isempty(trace) || ischar(trace) || isstring(trace)
        if nargin < 1, trace = ''; end
        r       = load_sim_results(char(trace));
        opts.dt = r.dt;
        trace   = r.blend.trace;
    end

    if ~isfolder(opts.out_dir), mkdir(opts.out_dir); end
    S   = load('trim_table_Poly_ConcatVer4p0.mat');
    sel = compute_selected(trace, S, opts);
    for v = opts.views
        animate_view(trace, opts, S, sel, v.axis, v.dims);
    end
end

%% -------------------------------------------------------------------------
function opts = fill_defaults(opts)
    opts.views      = struct('axis', {'lon', 'lon', 'lat', 'lat'}, ...
                             'dims', {[1 2], [3 4], [1 2], [3 4]});
    opts.dt         = 0.01;
    opts.wh_idx     = 2;
    opts.uh_list    = 1:20;
    opts.tube       = 'brt';
    opts.frame_step = 10;
    opts.fps        = 30;
    opts.out_dir    = 'visualize/figures';
end

function animate_view(trace, opts, S, sel, ax, dims)
    P      = axis_plot_config(ax);
    dim_x  = dims(1);  dim_y = dims(2);
    axU    = upper(ax);  tube = upper(opts.tube);
    prefix = tube_prefix(ax, opts.tube);

    state_native = trace_states_native(trace, ax);
    valid        = all(isfinite(state_native), 2);
    orig         = find(valid);                           % map filtered step -> trace index
    state_native = state_native(valid, :);
    traj         = state_native .* P.scale(:)';          % N x 4, plot units
    t_valid      = (trace.k(valid) - 1) * opts.dt;
    N            = size(state_native, 1);
    frames       = 1 : opts.frame_step : N;

    cache = containers.Map('KeyType', 'double', 'ValueType', 'any');
    get_table = @(uh) load_cached(cache, prefix, opts.wh_idx, uh);

    [xl, yl] = view_limits(S, opts, P, dim_x, dim_y, traj, sel(orig(frames)), get_table);

    fig = figure('Color', 'w', 'Position', [80 80 820 680], ...
                 'Visible', 'off', 'Renderer', 'painters', 'GraphicsSmoothing', 'off');
    axh = axes(fig); hold(axh, 'on'); grid(axh, 'on'); box(axh, 'on');
    xlim(axh, xl); ylim(axh, yl);
    xlabel(axh, P.labels{dim_x}); ylabel(axh, P.labels{dim_y});

    h_trail  = plot(axh, nan, nan, 'k-', 'LineWidth', 1.5);
    h_marker = scatter(axh, nan, nan, 70, [0.85 0.2 0.2], 'filled', 'MarkerEdgeColor', 'k');
    plot(axh, traj(1, dim_x), traj(1, dim_y), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    h_contour = gobjects(0);

    gif_path = fullfile(opts.out_dir, sprintf('brt_anim_%s_%s_%s.gif', axU, P.short{dim_x}, P.short{dim_y}));
    for i_frame = 1 : numel(frames)
        step   = frames(i_frame);
        uh_idx = sel(orig(step));
        tube_mat = get_table(uh_idx);
        delete(h_contour);  h_contour = gobjects(0);
        if ~isempty(tube_mat)
            X0        = S.XU0_interp(1:12, uh_idx, opts.wh_idx);
            state_dev = state_native(step, :)' - X0(P.trim_rows);
            [x_grid, y_grid, Vslice] = plane_slice(tube_mat.values, tube_mat.grid_axes, dim_x, dim_y, state_dev);
            x_grid = (x_grid + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
            y_grid = (y_grid + X0(P.trim_rows(dim_y))) * P.scale(dim_y);
            [~, h_contour] = contour(axh, x_grid, y_grid, Vslice', [0 0], 'Color', [0.1 0.5 0.9], 'LineWidth', 1.6);
        end
        set(h_trail,  'XData', traj(1:step, dim_x), 'YData', traj(1:step, dim_y));
        set(h_marker, 'XData', traj(step, dim_x),   'YData', traj(step, dim_y));
        title(axh, sprintf('%s %s  %s-%s   t = %.1f s   UH%d', ...
            axU, tube, P.short{dim_x}, P.short{dim_y}, t_valid(step), uh_idx));
        drawnow;
        write_gif(fig, gif_path, i_frame == 1, opts.fps);
    end
    fprintf('saved %s (%d frames)\n', gif_path, numel(frames));
    close(fig);
end

function [xl, yl] = view_limits(S, opts, P, dim_x, dim_y, traj, uh_used, get_table)
    xmin = min(traj(:, dim_x));  xmax = max(traj(:, dim_x));
    ymin = min(traj(:, dim_y));  ymax = max(traj(:, dim_y));
    for uh_idx = unique(uh_used(:)')
        tube_mat = get_table(uh_idx);
        if isempty(tube_mat), continue; end
        X0 = S.XU0_interp(1:12, uh_idx, opts.wh_idx);
        gx = (tube_mat.grid_axes{dim_x} + X0(P.trim_rows(dim_x))) * P.scale(dim_x);
        gy = (tube_mat.grid_axes{dim_y} + X0(P.trim_rows(dim_y))) * P.scale(dim_y);
        xmin = min(xmin, min(gx));  xmax = max(xmax, max(gx));
        ymin = min(ymin, min(gy));  ymax = max(ymax, max(gy));
    end
    xl = pad_range(xmin, xmax);  yl = pad_range(ymin, ymax);
end

function lim = pad_range(lo, hi)
    m = 0.05 * max(hi - lo, eps);
    lim = [lo - m, hi + m];
end

function tube_mat = load_cached(cache, prefix, wh_idx, uh_idx)
    if isKey(cache, uh_idx), tube_mat = cache(uh_idx); return; end
    fname = sprintf('%s_UH%d_WH%d.mat', prefix, uh_idx, wh_idx);
    if isfile(fname)
        tube_mat = load(fname);
        tube_mat.values = double(tube_mat.values);
    else
        tube_mat = [];
    end
    cache(uh_idx) = tube_mat;
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
            error('visualize_BRT_animation: unknown axis "%s".', ax);
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

function sel = compute_selected(trace, S, opts)
    % Replicate the filter's curr/next anchor scheduling so the drawn tube is the
    % one V(t) was logged against: nearest UH node is curr, one node up is next,
    % and next is taken once it is reachable (in-grid, V<0) on every checked axis.
    UH    = S.UH;
    ulist = opts.uh_list;
    chk_axes = {'lon', 'lat'};
    chk = struct('spec', {}, 'state', {}, 'prefix', {}, 'cache', {});
    for i = 1:numel(chk_axes)
        chk(i).spec   = FilterConfig.axisSpec(chk_axes{i});
        chk(i).state  = trace_states_native(trace, chk_axes{i});
        chk(i).prefix = tube_prefix(chk_axes{i}, opts.tube);
        chk(i).cache  = containers.Map('KeyType', 'double', 'ValueType', 'any');
    end

    u_body = trace.uBody(:);
    N      = numel(u_body);
    sel    = zeros(N, 1);
    for s = 1:N
        uhA      = max(UH(1), min(UH(end), u_body(s)));
        [~, j]   = min(abs(UH(ulist) - uhA));
        cid      = ulist(j);
        nid      = min(ulist(end), cid + 1);
        ready    = true;
        for i = 1:numel(chk)
            tube_mat = load_cached(chk(i).cache, chk(i).prefix, opts.wh_idx, nid);
            if isempty(tube_mat), ready = false; break; end
            ga   = tube_mat.grid_axes;
            X0   = S.XU0_interp(1:12, nid, opts.wh_idx);
            x    = chk(i).state(s, :)' - X0(chk(i).spec.trim_rows(:));
            gmin = cellfun(@(c) c(1),   ga(:));
            gmax = cellfun(@(c) c(end), ga(:));
            xc   = min(max(x, gmin), gmax);
            V    = interpn(ga{1}, ga{2}, ga{3}, ga{4}, tube_mat.values, xc(1), xc(2), xc(3), xc(4), 'linear');
            if ~(all(x >= gmin) && all(x <= gmax) && isfinite(V) && V < 0), ready = false; break; end
        end
        if (cid == nid) || ready, sel(s) = nid; else, sel(s) = cid; end
    end
end

function prefix = tube_prefix(ax, tube)
    axU = upper(ax);  tubeU = upper(tube);
    prefix = sprintf('reachability_data/guam_output/%s_%s/GUAM_%s_%s', axU, tubeU, axU, tubeU);
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

function write_gif(fig, gif_path, is_first, fps)
    [gif_frame, colour_map] = rgb2ind(frame2im(getframe(fig)), 256);
    if is_first
        imwrite(gif_frame, colour_map, gif_path, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
    else
        imwrite(gif_frame, colour_map, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end
