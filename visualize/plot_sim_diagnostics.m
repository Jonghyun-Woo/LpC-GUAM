function figs = plot_sim_diagnostics(source)
    % Basic + filter diagnostic figures for a saved run. source: a data struct
    % (SimLogger.exportData), a results .mat path, or empty for the default path.
    if nargin < 1, source = ''; end
    if isstruct(source) && isfield(source, 'state')
        data = source;
    else
        r = load_sim_results(char(source));
        data = r.blend.data;
    end
    figs = plot_basic(data);
    figs = plot_filter(data, figs);
end

function figs = plot_basic(data)
    time = data.t;
    state = data.state;
    ft2m = 0.3048;
    figs = struct();
    line_width = 1.2;
    fig_size = [1000, 630];

    figs.position = figure('Name', 'Position (NED)', 'Position', [100 100 fig_size]);
    labels = {'North [m]', 'East [m]', 'Down [m]'};
    for i = 1:3
        subplot(3, 1, i);
        plot(time, ft2m * state(:, i), 'b', time, ft2m * data.ref_pos(:, i), 'r--', 'LineWidth', line_width);
        ylabel(labels{i}); grid on;
        if i == 1, legend('Response', 'Reference Trajectory', 'Location', 'southeast'); title('Inertial position'); end
    end
    xlabel('Time [s]');

    figs.velocity = figure('Name', 'Body velocity', 'Position', [100 100 fig_size]);
    labels = {'u [m/s]', 'v [m/s]', 'w [m/s]'};
    for i = 1:3
        subplot(3, 1, i);
        plot(time, ft2m * state(:, 3 + i), 'b', time, ft2m * data.ref_vel(:, i), 'r--', 'LineWidth', line_width);
        ylabel(labels{i}); grid on;
        if i == 1, legend('Response', 'Reference Trajectory', 'Location', 'southeast'); title('Body-fixed frame Velocity'); end
    end
    xlabel('Time [s]');

    figs.attitude = figure('Name', 'Attitude', 'Position', [100 100 fig_size]);
    labels = {'\phi [deg]', '\theta [deg]', '\psi [deg]'};
    for i = 1:3
        subplot(3, 1, i);
        plot(time, rad2deg(state(:, 6 + i)), 'b', 'LineWidth', line_width);
        ylabel(labels{i}); grid on;
        if i == 1, title('Euler angles'); end
    end
    xlabel('Time [s]');

    figs.effectors = figure('Name', 'Effectors', 'Position', [100 100 fig_size]);
    subplot(2, 1, 1);
    rad2rpm = 60 / (2 * pi);
    hold on; grid on;
    for i = 1:9
        plot(time, data.engine(:, i) .* rad2rpm, 'LineWidth', line_width);
    end
    ylabel('Rotors [RPM]');
    title('Actuator response');
    legend('\Pi_1', '\Pi_2', '\Pi_3', '\Pi_4', '\Pi_5', '\Pi_6', '\Pi_7', '\Pi_8', '\Pi_p', 'Location', 'northeastoutside');
    
    subplot(2, 1, 2);
    hold on; grid on;
    plot(time, rad2deg(data.surface), 'LineWidth', line_width);
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
    figs.brt_value = figure('Name', 'BRT value', 'Color', 'w', 'Position', [100 100 fig_size]);
    hold on; grid on;
    yline(0, '--', 'V=0');
    yl = ylim;
    h_viol = patch(nan, nan, [1 0 0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    h_actv = patch(nan, nan, [0 0.6 0], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    shade_intervals(t, data.active == 1, [0 0.6 0], 0.2, yl);
    shade_intervals(t, brtV > 0, [1 0 0], 0.5, yl);
    h_line = plot(t, brtV, 'LineWidth', line_width);

    legend([h_line, h_viol, h_actv], {'V(t)', 'V>0 (violation)', 'filter active'}, 'Location', 'northwest');
    ylim([min(brtV) - ylim_gap, max(brtV) + ylim_gap]);
    xlabel('Time [s]');
    ylabel('V');
    title('BRT value V(t)');
end

function shade_intervals(t, mask, color, alpha, yl)
    mask = logical(mask(:)');
    if ~any(mask), return; end
    edges = diff([false, mask, false]);
    starts = find(edges == 1);
    stops = find(edges == -1) - 1;
    for i = 1:numel(starts)
        a = starts(i);  b = stops(i);
        x0 = t(a);  x1 = t(b);
        if a > 1, x0 = (t(a - 1) + t(a)) / 2; end
        if b < numel(t), x1 = (t(b) + t(b + 1)) / 2; end
        patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], color, ...
            'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
