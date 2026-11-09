close all; clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);

idx_list = 1:28;
trimData = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), 'XU0_interp', 'WH');
wh_idx = find(trimData.WH == 0, 1);
trim = trimData.XU0_interp(:, idx_list, wh_idx);

ft2m = 0.3048;
r2d = 180 / pi;
u_trim = trim(1, :) * ft2m;
theta_trim = trim(11, :) * r2d;
surface_trim = trim(13:16, :)' * r2d;
lift_rpm_trim = trim(17:24, :)';
pusher_rpm_trim = trim(25, :);
line_width = 1.5;
fig_size = [800, 480];
tick_font_size = 13;
title_font_size = 15;
label_font_size = 13;
legend_font_size = 12;

thetaFig = figure();
thetaFig.Theme = 'light';
thetaFig.Position(3:4) = fig_size;
hold on; grid on;
set(gca, 'FontSize', tick_font_size);
plot(u_trim, theta_trim, 'o-', 'LineWidth', line_width);
xlabel('$u$ [m/s]', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel('$\theta^*$ [deg]', 'Interpreter', 'latex', 'FontSize', label_font_size);
title('Pitch Angle $\theta^*$', 'Interpreter', 'latex', 'FontSize', title_font_size);
phase_patches = shade_phases(gca, u_trim);
legend(phase_patches, {'Hovering', 'Transition', 'Cruise'}, 'Location', 'northwest', 'FontSize', legend_font_size);

ctrlsrfFig = figure();
ctrlsrfFig.Theme = 'light';
ctrlsrfFig.Position(3:4) = fig_size;
hold on; grid on;
set(gca, 'FontSize', tick_font_size);
plot(u_trim, surface_trim, 'o-', 'LineWidth', line_width);
xlabel('$u$ [m/s]', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel('$\delta^*$ [deg]', 'Interpreter', 'latex', 'FontSize', label_font_size);
legend({'flap', 'aileron', 'elevator', 'rudder'}, 'Location', 'southeast', 'FontSize', legend_font_size);
title('Control Surface', 'Interpreter', 'latex', 'FontSize', title_font_size);
shade_phases(gca, u_trim);

rotorFig = figure();
rotorFig.Theme = 'light';
rotorFig.Position(3:4) = fig_size;
hold on; grid on;
set(gca, 'FontSize', tick_font_size);
plot(u_trim, lift_rpm_trim, 'o-', 'LineWidth', line_width);
xlabel('$u$ [m/s]', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel('$\Pi_L^*$ [RPM]', 'Interpreter', 'latex', 'FontSize', label_font_size);
legend(compose('lift%d', 1:8), 'Location', 'northeast', 'FontSize', legend_font_size);
title('Lift Propeller', 'Interpreter', 'latex', 'FontSize', title_font_size);
shade_phases(gca, u_trim);

pusherFig = figure();
pusherFig.Theme = 'light';
pusherFig.Position(3:4) = fig_size;
hold on; grid on;
set(gca, 'FontSize', tick_font_size);
plot(u_trim, pusher_rpm_trim, 'o-', 'LineWidth', line_width);
xlabel('$u$ [m/s]', 'Interpreter', 'latex', 'FontSize', label_font_size);
ylabel('$\Pi_P^*$ [RPM]', 'Interpreter', 'latex', 'FontSize', label_font_size);
title('Pusher Propeller', 'Interpreter', 'latex', 'FontSize', title_font_size);
shade_phases(gca, u_trim);

figs = {thetaFig, ctrlsrfFig, rotorFig, pusherFig};
names = {'trim_theta', 'trim_control_surface', 'trim_lift_propeller', 'trim_pusher_propeller'};
for k = 1:numel(figs)
    exportgraphics(figs{k}, fullfile(here, 'figures', [names{k}, '.png']), 'Resolution', 300);
end

function patches = shade_phases(ax, u_trim)
    phase_last_idx = [1, 20, numel(u_trim)];
    phase_colors = [0.20 0.45 0.85; 0.95 0.70 0.15; 0.25 0.65 0.35];
    mid_u = (u_trim(1:end-1) + u_trim(2:end)) / 2;
    edges_u = [u_trim(1) - (mid_u(1) - u_trim(1)), mid_u(phase_last_idx(1:2)), u_trim(end) + (u_trim(end) - mid_u(end))];
    xlim(ax, [edges_u(1), edges_u(end)]);
    y_limits = ylim(ax);
    ylim(ax, y_limits);
    patches = gobjects(1, 3);
    for k = 1:3
        x_left = edges_u(k);
        x_right = edges_u(k + 1);
        patches(k) = patch(ax, [x_left x_right x_right x_left], y_limits([1 1 2 2]), phase_colors(k, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(patches(k), 'bottom');
    end
end
