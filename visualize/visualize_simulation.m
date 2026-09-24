function figs = visualize_simulation(source)
    % Basic + filter diagnostic figures for a saved run. source: a data struct
    % (SimLogger.exportData), a results .mat path, or empty for the default path.
    % Every figure is also saved as a PNG under visualize/figures/.
    close all;
    if nargin < 1, source = ''; end
    if isstruct(source) && isfield(source, 'state')
        data = source;
        data_off = [];
    else
        r = load_sim_results(char(source));
        data = r.blend.data;
        if isfield(r, 'off'), data_off = r.off.data; else, data_off = []; end
    end
    figs = plot_basic(data, data_off);
    figs = plot_filter(data, figs);
    figs = plot_brt_envelope(data, figs);
    save_figures(figs);
end

function save_figures(figs)
    outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
    if ~isfolder(outdir), mkdir(outdir); end
    names = fieldnames(figs);
    for i = 1:numel(names)
        f = figs.(names{i});
        if ~isa(f, 'matlab.ui.Figure'), continue; end
        exportgraphics(f, fullfile(outdir, [names{i} '.png']), 'Resolution', 200);
    end
end

function figs = plot_basic(data, data_off)
    if nargin < 2, data_off = []; end
    time = data.t;
    state = data.state;
    ft2m = 0.3048;
    figs = struct();
    line_width = 1.2;
    fig_size = [640, 720];
    nom_style = {'Color', [0.35 0.35 0.35], 'LineWidth', 1.5};

    figs.position = figure('Name', 'Position (NED)'); figs.position.Position(3:4) = fig_size;
    tl = tiledlayout(figs.position, 3, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    labels = {'North [m]', 'East [m]', 'Down [m]'};
    for i = 1:3
        nexttile(tl); hold on; grid on;
        shade_active_band(time, data.active);
        plot(time, ft2m * data.ref_pos(:, i), 'r--', 'LineWidth', line_width, 'DisplayName', 'Reference Trajectory');
        plot(time, ft2m * state(:, i), 'b', 'LineWidth', line_width, 'DisplayName', 'Safety Filter');
        if ~isempty(data_off)
            plot(time, ft2m * data_off.state(:, i), nom_style{:}, 'DisplayName', 'Nominal');
        end
        ylabel(labels{i});
        if i == 1, legend('Location', 'southeast'); title('Inertial position'); end
    end
    xlabel('Time [s]');

    figs.velocity = figure('Name', 'Body velocity'); figs.velocity.Position(3:4) = fig_size;
    tl = tiledlayout(figs.velocity, 3, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    ax_vel = gobjects(1, 3);
    labels = {'u [m/s]', 'v [m/s]', 'w [m/s]'};
    for i = 1:3
        ax_vel(i) = nexttile(tl); hold on; grid on;
        shade_active_band(time, data.active);
        plot(time, ft2m * data.ref_vel(:, i), 'r--', 'LineWidth', line_width, 'DisplayName', 'Reference Trajectory');
        plot(time, ft2m * state(:, 3 + i), 'b', 'LineWidth', line_width, 'DisplayName', 'Safety Filter');
        if ~isempty(data_off)
            plot(time, ft2m * data_off.state(:, 3 + i), nom_style{:}, 'DisplayName', 'Nominal');
        end
        ylabel(labels{i});
        if i == 1, legend('Location', 'northwest'); title('Body-fixed frame Velocity'); end
    end
    xlabel('Time [s]');
    figs.ax.velocity = ax_vel;

    figs.attitude = figure('Name', 'Attitude'); figs.attitude.Position(3:4) = fig_size;
    tl = tiledlayout(figs.attitude, 3, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    ax_att = gobjects(1, 3);
    labels = {'\phi [deg]', '\theta [deg]', '\psi [deg]'};
    for i = 1:3
        ax_att(i) = nexttile(tl); hold on; grid on;
        shade_active_band(time, data.active);
        plot(time, rad2deg(state(:, 6 + i)), 'b', 'LineWidth', line_width, 'DisplayName', 'Safety Filter');
        if ~isempty(data_off)
            plot(time, rad2deg(data_off.state(:, 6 + i)), nom_style{:}, 'DisplayName', 'Nominal');
        end
        ylabel(labels{i});
        if i == 1, legend('Location', 'northwest'); title('Euler angles'); end
    end
    xlabel('Time [s]');
    figs.ax.attitude = ax_att;

    figs.actuators = figure('Name', 'Actuators'); figs.actuators.Position(3:4) = [1000, 630];
    tl = tiledlayout(figs.actuators, 2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    rad2rpm = 60 / (2 * pi);
    rpm = data.engine .* rad2rpm;
    nexttile(tl); hold on; grid on;
    yl = pad_lim(rpm(:));
    h_act = shade_active(time, data.active, yl);
    rotors = gobjects(1, 9);
    for i = 1:9
        rotors(i) = plot(time, rpm(:, i), 'LineWidth', line_width);
    end
    ylim(yl);
    ylabel('Rotors [RPM]');
    title('Actuators');
    names = {'\Pi_1', '\Pi_2', '\Pi_3', '\Pi_4', '\Pi_5', '\Pi_6', '\Pi_7', '\Pi_8', '\Pi_p'};
    if isempty(h_act)
        legend(rotors, names, 'Location', 'northeastoutside');
    else
        legend([rotors, h_act], [names, {'Filter active'}], 'Location', 'northeastoutside');
    end

    surf_deg = rad2deg(data.surface);
    nexttile(tl); hold on; grid on;
    yl = pad_lim(surf_deg(:));
    shade_active(time, data.active, yl);
    plot(time, surf_deg, 'LineWidth', line_width);
    ylim(yl);
    ylabel('Surfaces [deg]');  legend('LA', 'RA', 'LE', 'RE', 'RUD', 'Location', 'northeastoutside');
    xlabel('Time [s]');
end

function figs = plot_filter(data, figs)
    k = data.k;
    brtV = data.brtV;
    if all(isnan(brtV)), return; end
    line_width = 1.2;
    ylim_gap = 0.02;
    fig_size = [1000, 630]; 

    t = data.dt .* k;
    figs.brt_value = figure('Name', 'BRT value', 'Color', 'w'); figs.brt_value.Position(3:4) = fig_size;
    hold on; grid on;
    yline(0, '--', 'V=0');
    yl = ylim;
    h_viol = patch(nan, nan, [1 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h_actv = patch(nan, nan, [0 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    shade_intervals(t, data.active == 1, [0 0.6 0], 0.15, yl);
    shade_intervals(t, brtV > 0, [1 0 0], 0.5, yl);
    h_line = plot(t, brtV, 'LineWidth', line_width);

    legend([h_line, h_viol, h_actv], {'V(t)', 'V>0 (violation)', 'filter active'}, 'Location', 'northwest');
    ylim([min(brtV) - ylim_gap, max(brtV) + ylim_gap]);
    xlabel('Time [s]');
    ylabel('V');
    title('BRT value V(t)');
end

function figs = plot_brt_envelope(data, figs)
    % Overlay the per-state BRT safe interval (lb/ub where V=0 along that axis,
    % holding the other perturbation states fixed) onto the velocity/attitude
    % figures. Reconstructs the anchor scheduling and value function from a
    % default Config; skips silently if the BRT tables are unavailable. The
    % anchor is the filter's selected curr/next frame, so the bounds are
    % consistent with the logged V(t) in figs.brt_value.
    ctrl = try_build_controller();
    if isempty(ctrl), return; end

    ft2m = 0.3048;  r2d = 180 / pi;
    %  axis  slot  figure       subplot  scale
    targets = {
        'lon', 1, 'velocity', 1, ft2m
        'lat', 1, 'velocity', 2, ft2m
        'lon', 2, 'velocity', 3, ft2m
        'lat', 4, 'attitude', 1, r2d
        'lon', 4, 'attitude', 2, r2d };

    N      = size(data.state, 1);
    stride = max(1, round(N / 250));
    idx    = 1:stride:N;
    states = data.state(idx, :);
    t      = data.t(idx);
    sel    = selected_uh(ctrl, states);

    for row = 1:size(targets, 1)
        ax = targets{row, 1};  slot = targets{row, 2};
        fig = targets{row, 3};  sub = targets{row, 4};  sc = targets{row, 5};
        if ~isfield(ctrl.safety_filter.value_function, ax), continue; end

        [lo, hi] = brt_bounds(ctrl, states, sel, ax, slot);
        if all(isnan(lo)) && all(isnan(hi)), continue; end

        ax = figs.ax.(fig)(sub); hold(ax, 'on');
        plot(ax, t, sc * hi, 'm--', 'LineWidth', 1.0, 'DisplayName', 'BRT bound');
        plot(ax, t, sc * lo, 'm--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
end

function ctrl = try_build_controller()
    ctrl = [];
    try
        hub  = Config('', struct());
        ctrl = Controller(hub.controller, hub.sim.dt);
        if isempty(ctrl.safety_filter) || ~isstruct(ctrl.safety_filter.value_function)
            ctrl = [];
        end
    catch
        ctrl = [];
    end
end

function sel = selected_uh(ctrl, states)
    % Replicate the filter's joint curr/next anchor selection (LivenessFilter):
    % use_next = same_anchor || transition_ready, where transition_ready needs
    % the next anchor valid with V<0 across every active axis. Returns the
    % selected UH index per timestep so the bounds match the logged V(t).
    UH    = ctrl.baseline_controller.UH;
    wh    = ctrl.safety_filter_wh_anchor;
    axes  = ctrl.safety_filter.axes;
    n     = size(states, 1);
    sel   = zeros(n, 1);
    for m = 1:n
        s   = states(m, :)';
        uhA = max(UH(1), min(UH(end), s(4)));
        [~, cid] = min(abs(UH(:) - uhA));
        nid = min(numel(UH), cid + 1);

        transition_ready = true;
        for a = 1:numel(axes)
            vf   = ctrl.safety_filter.value_function.(axes{a});
            spec = FilterConfig.axisSpec(axes{a});
            uh   = UH(nid);
            [X0a, ~] = ctrl.interp_xu0(uh, wh);
            xa   = s(spec.state_rows) - X0a(spec.trim_rows);
            [V, ~, ok] = vf.query(xa, uh, wh);
            inside = all(xa >= vf.grid_min - 1e-12) && all(xa <= vf.grid_max + 1e-12);
            transition_ready = transition_ready && (ok && inside && isfinite(V) && V < 0);
        end
        if (cid == nid) || transition_ready, sel(m) = nid; else, sel(m) = cid; end
    end
end

function [lo, hi] = brt_bounds(ctrl, states, sel, ax, slot)
    NS   = 121;
    UH   = ctrl.baseline_controller.UH;
    wh   = ctrl.safety_filter_wh_anchor;
    vf   = ctrl.safety_filter.value_function.(ax);
    spec = FilterConfig.axisSpec(ax);

    n  = size(states, 1);
    lo = nan(n, 1);  hi = nan(n, 1);
    if ~vf.available, return; end

    gmin = vf.grid_min(slot);  gmax = vf.grid_max(slot);
    xd   = linspace(gmin, gmax, NS)';

    for m = 1:n
        s  = states(m, :)';
        uh = UH(sel(m));

        [X0a, ~] = ctrl.interp_xu0(uh, wh);
        xa = s(spec.state_rows) - X0a(spec.trim_rows);

        Vprof = nan(NS, 1);  ok_all = true;
        for p = 1:NS
            xq = xa;  xq(slot) = xd(p);
            [V, ~, ok] = vf.query(xq, uh, wh);
            if ~ok, ok_all = false; break; end
            Vprof(p) = V;
        end
        if ~ok_all, continue; end

        [lp, hp] = safe_interval(xd, Vprof, gmin, gmax);
        trimval  = X0a(spec.trim_rows(slot));
        lo(m) = trimval + lp;
        hi(m) = trimval + hp;
    end
end

function [lo, hi] = safe_interval(xd, V, gmin, gmax)
    % Tube boundary [lo, hi] along the axis: the outermost V=0 crossings of the
    % slice, independent of where the current state sits. NaN only if V>0 over
    % the whole axis (no tube at this slice).
    lo = nan;  hi = nan;
    zc = [];
    for p = 1:numel(xd) - 1
        if V(p) * V(p + 1) < 0
            zc(end + 1) = xd(p) + (0 - V(p)) * (xd(p + 1) - xd(p)) / (V(p + 1) - V(p)); %#ok<AGROW>
        end
    end
    if isempty(zc)
        if all(V <= 0), lo = gmin;  hi = gmax; end
        return;
    end
    lo = min(zc);  hi = max(zc);
end

function h = shade_active(t, active, yl)
    % Green background over the intervals where the liveness filter was active
    % (same convention as figs.brt_value). Drawn before the traces so it sits
    % behind them; returns a legend patch handle, or [] if never active.
    h = [];
    if isempty(active) || all(isnan(active)) || ~any(active == 1), return; end
    shade_intervals(t, active == 1, [0 0.6 0], 0.1, yl);
    h = patch(nan, nan, [0 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

function yl = pad_lim(vals)
    lo = min(vals);  hi = max(vals);
    if ~isfinite(lo) || ~isfinite(hi) || lo == hi, yl = [lo - 1, hi + 1]; return; end
    pad = 0.03 * (hi - lo);
    yl = [lo - pad, hi + pad];
end

function shade_active_band(t, active)
    % Full-height green bands over liveness-filter-active intervals. Uses
    % xregion so the bands track the axis limits (unlike the fixed-height
    % effector shading) and stay out of the legend.
    if isempty(active) || all(isnan(active)) || ~any(active == 1), return; end
    [x0, x1] = active_intervals(t, active == 1);
    if isempty(x0), return; end
    hr = xregion(x0, x1, 'FaceColor', [0 0.6 0], 'FaceAlpha', 0.15);
    set(hr, 'HandleVisibility', 'off');
end

function shade_intervals(t, mask, color, alpha, yl)
    [x0, x1] = active_intervals(t, mask);
    for i = 1:numel(x0)
        patch([x0(i) x1(i) x1(i) x0(i)], [yl(1) yl(1) yl(2) yl(2)], color, ...
            'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function [x0, x1] = active_intervals(t, mask)
    % Contiguous true-runs of mask as time intervals [x0, x1], each extended to
    % the midpoint of the neighbouring samples so adjacent bands meet cleanly.
    mask = logical(mask(:)');
    x0 = [];  x1 = [];
    if ~any(mask), return; end
    edges = diff([false, mask, false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    x0 = zeros(numel(starts), 1);  x1 = zeros(numel(starts), 1);
    for i = 1:numel(starts)
        a = starts(i);  b = stops(i);
        xa = t(a);  xb = t(b);
        if a > 1, xa = (t(a - 1) + t(a)) / 2; end
        if b < numel(t), xb = (t(b) + t(b + 1)) / 2; end
        x0(i) = xa;  x1(i) = xb;
    end
end
