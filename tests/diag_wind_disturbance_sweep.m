clear; close all; clc;

cfg = Config();
GUAM = LpC_GUAM(cfg);

load('controller/trim_table_Poly_ConcatVer4p0.mat');

max_dF = nan(3, 20);
min_dF = nan(3, 20);
max_dM = nan(3, 20);
min_dM = nan(3, 20);
% (alpha, beta) [deg] at which each extreme occurs: (component, idx, [alpha beta])
max_dF_ab = nan(3, 20, 2);
min_dF_ab = nan(3, 20, 2);
max_dM_ab = nan(3, 20, 2);
min_dM_ab = nan(3, 20, 2);
% body-frame wind unit direction [x_fwd; y_right; z_down] at each extreme (|wind| = 10 m/s):
% (component, idx, [wx wy wz])
max_dF_wind = nan(3, 20, 3);
min_dF_wind = nan(3, 20, 3);
max_dM_wind = nan(3, 20, 3);
min_dM_wind = nan(3, 20, 3);
WH = 3;

for idx = 1 : 20
    engine0 = XU0_interp(end-8:end, idx, WH);
    surface0 = XU0_interp(end-12:end-9, idx, WH);
    flp = surface0(1);
    ail = surface0(2);
    ele = surface0(3);
    rud = surface0(4);
    surface0 = [flp - ail;...
                flp + ail;...
                ele;...
                ele;...
                rud];
    Vb0     = XU0_interp(1:3, idx, WH);
    omega0  = XU0_interp(4:6, idx, WH);
    euler0  = XU0_interp(10:12, idx, WH);
    pos0    = [0; 0; -100];
    X0      = [pos0; Vb0; euler0; omega0];
    
    GUAM.reset(X0, engine0, surface0);
    dx0 = GUAM.state_derivative(X0, engine0, surface0);
    
    % Wind disturbance sweep at trim: 10 m/s body-frame wind, all incidence directions
    Vw      = 10 * GUAM.units.m;                    % 10 m/s expressed in ft/s (sim base unit)
    R_i2b   = RSLQR.rotm_i2b(X0(7), X0(8), X0(9));  % trim attitude (NED->body)
    R_b2i   = R_i2b';
    mass    = GUAM.vehicleConfig.mass;
    Inertia = GUAM.vehicleConfig.I;
    [rho, a_snd] = GUAM.environment.atmosphere(-X0(3));  % same atmosphere state_derivative uses
    
    alpha_deg = -180:5:180;    % wind AoA sweep
    beta_deg  = -180:5:180;    % wind AoS sweep
    N_alpha = numel(alpha_deg);  N_beta = numel(beta_deg);
    dVel_b = zeros(3, N_alpha, N_beta);   % disturbance acceleration           in body frame
    dOmg_b = zeros(3, N_alpha, N_beta);   % disturbance angular acceleration   in body frame
    valid = true(N_alpha, N_beta);      % aero-database validity of each sweep point

    for i_alpha = 1:N_alpha
        for i_beta = 1:N_beta
            alpha_rad = deg2rad(alpha_deg(i_alpha));
            beta_rad  = deg2rad(beta_deg(i_beta));
            % 10 m/s wind, direction (alpha, beta) in body axes (aero convention)
            w_body = Vw * [cos(alpha_rad)*cos(beta_rad); sin(beta_rad); sin(alpha_rad)*cos(beta_rad)];
            GUAM.environment.wind_ned = R_b2i * w_body;   % inject as body-frame wind

            ddx = GUAM.state_derivative(X0, engine0, surface0) - dx0;
            dVel_b(:, i_alpha, i_beta) = ddx(4:6);     % dF <- acc
            dOmg_b(:, i_alpha, i_beta) = ddx(10:12);   % dM <- omega_dot

            % Re-run the aero model to read back its validity flags.
            v_air = X0(4:6) - w_body;
            [~, ~, Validity] = run_LpC_aero([v_air; X0(10:12)], engine0, surface0, ...
                                            rho, a_snd, GUAM.units, GUAM.vehicleConfig);
            valid(i_alpha, i_beta) = ~any(Validity);
        end
    end
    GUAM.environment.wind_ned = zeros(3, 1);   % restore calm air
    
    % Drop the invalid aerodynamic points.
    dVel_b(:, ~valid) = NaN;
    dOmg_b(:, ~valid) = NaN;
    fprintf('aero-database validity: %d / %d sweep points excluded (%.1f%%)\n', ...
            nnz(~valid), numel(valid), 100 * nnz(~valid) / numel(valid));
    for comp = 1:3
        slice_F = reshape(dVel_b(comp, :, :), N_alpha, N_beta);
        slice_M = reshape(dOmg_b(comp, :, :), N_alpha, N_beta);

        % dF: on ties, take the smaller-sideslip direction (beta near 0 / 180 deg)
        [max_dF(comp, idx), a_at, b_at] = extreme_ab(slice_F, true,  'min', alpha_deg, beta_deg);
        max_dF_ab(comp, idx, :)   = [a_at, b_at];
        max_dF_wind(comp, idx, :) = wind_body_unit(a_at, b_at);
        [min_dF(comp, idx), a_at, b_at] = extreme_ab(slice_F, false, 'min', alpha_deg, beta_deg);
        min_dF_ab(comp, idx, :)   = [a_at, b_at];
        min_dF_wind(comp, idx, :) = wind_body_unit(a_at, b_at);

        % dM: on ties, take the larger-sideslip direction (beta near +/-90 deg)
        [max_dM(comp, idx), a_at, b_at] = extreme_ab(slice_M, true,  'max', alpha_deg, beta_deg);
        max_dM_ab(comp, idx, :)   = [a_at, b_at];
        max_dM_wind(comp, idx, :) = wind_body_unit(a_at, b_at);
        [min_dM(comp, idx), a_at, b_at] = extreme_ab(slice_M, false, 'max', alpha_deg, beta_deg);
        min_dM_ab(comp, idx, :)   = [a_at, b_at];
        min_dM_wind(comp, idx, :) = wind_body_unit(a_at, b_at);
    end
end

save('guam_disturbance_lb_ub_10mps.mat',...
        'max_dF',...
        'min_dF',...
        'max_dM',...
        'min_dM',...
        'max_dF_ab',...
        'min_dF_ab',...
        'max_dM_ab',...
        'min_dM_ab',...
        'max_dF_wind',...
        'min_dF_wind',...
        'max_dM_wind',...
        'min_dM_wind');

% Disturbance magnitude profiles over the (alpha, beta) grid
[BB, AA] = meshgrid(beta_deg, alpha_deg);
Acc_mag = squeeze(vecnorm(dVel_b, 2, 1));
Ang_mag = squeeze(vecnorm(dOmg_b, 2, 1));

font_size = 20;
figure();
hold on; grid on; 
scatter3(squeeze(dVel_b(1, :, :)), squeeze(dVel_b(2, :, :)), squeeze(dVel_b(3, :, :)), '.');
view(3);
xlabel('$\dot{u}$ [ft/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size); 
ylabel('$\dot{v}$ [ft/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size); 
zlabel('$\dot{w}$ [ft/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size);

figure();
hold on; grid on; 
scatter3(squeeze(dOmg_b(1, :, :)), squeeze(dOmg_b(2, :, :)), squeeze(dOmg_b(3, :, :)), '.');
view(3);
xlabel('$\dot{p}$ [rad/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('$\dot{q}$ [rad/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size);
zlabel('$\dot{r}$ [rad/$s^2$]', 'Interpreter', 'latex', 'FontSize', font_size);

font_size = 14;
figure('Name', 'Wind disturbance profile (10 m/s body-frame incidence sweep)');
subplot(1, 2, 1);
contourf(BB, AA, Acc_mag, 20, 'LineColor', 'none'); colorbar;
xlabel('wind AoS $\beta$ [deg]', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('wind AoA $\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', font_size);
title('$|\Delta a_b|$  Acceleration Difference', 'Interpreter', 'latex', 'FontSize', font_size); axis tight;

subplot(1, 2, 2);
contourf(BB, AA, Ang_mag, 20, 'LineColor', 'none'); colorbar;
xlabel('wind AoS $\beta$ [deg]', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('wind AoA $\alpha$ [deg]', 'Interpreter', 'latex', 'FontSize', font_size);
title('$|\Delta \dot{\omega}_b|$  Angular Acceleration Difference', 'Interpreter', 'latex', 'FontSize', font_size); axis tight;

font_size = 12;
figure('Name', 'Max / Min Wind Disturbance Magnitudes over Trim Sweep');
indices = 1:20;

subplot(2, 2, 1);
hold on; grid on;
plot(indices, max_dF(1, :), '-o', 'DisplayName', 'Max $\dot{u}$'); 
plot(indices, max_dF(2, :), '-o', 'DisplayName', 'Max $\dot{v}$'); 
plot(indices, max_dF(3, :), '-o', 'DisplayName', 'Max $\dot{w}$'); 
legend('Location', 'northwest', 'Interpreter', 'latex');
xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Max Disturbance Acceleration', 'Interpreter', 'latex', 'FontSize', font_size);

subplot(2, 2, 2);
hold on; grid on;
plot(indices, min_dF(1, :), '-o', 'DisplayName', 'Min $\dot{u}$'); 
plot(indices, min_dF(2, :), '-o', 'DisplayName', 'Min $\dot{v}$'); 
plot(indices, min_dF(3, :), '-o', 'DisplayName', 'Min $\dot{w}$'); 
legend('Location', 'best', 'Interpreter', 'latex');
xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Min Disturbance Acceleration', 'Interpreter', 'latex', 'FontSize', font_size);

subplot(2, 2, 3);
hold on; grid on;
plot(indices, max_dM(1, :), '-o', 'DisplayName', 'Max $\dot{p}$'); 
plot(indices, max_dM(2, :), '-o', 'DisplayName', 'Max $\dot{q}$'); 
plot(indices, max_dM(3, :), '-o', 'DisplayName', 'Max $\dot{r}$'); 
legend('Location', 'northeast', 'Interpreter', 'latex');
xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Max Disturbance Angular Acceleration', 'Interpreter', 'latex', 'FontSize', font_size);

subplot(2, 2, 4);
hold on; grid on;
plot(indices, min_dM(1, :), '-o', 'DisplayName', 'Min $\dot{p}$');
plot(indices, min_dM(2, :), '-o', 'DisplayName', 'Min $\dot{q}$');
plot(indices, min_dM(3, :), '-o', 'DisplayName', 'Min $\dot{r}$');
legend('Location', 'southeast', 'Interpreter', 'latex');
xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
ylabel('Min Disturbance Angular Acceleration', 'Interpreter', 'latex', 'FontSize', font_size);

%%
% Per-component figures: extreme (max/min) value and the body-frame wind direction
% (unit vector, |wind| = 10 m/s) that produced it.
comp_name = {'\dot{u}', '\dot{v}', '\dot{w}', '\dot{p}', '\dot{q}', '\dot{r}'};
comp_unit = {'ft/$s^2$', 'ft/$s^2$', 'ft/$s^2$', 'rad/$s^2$', 'rad/$s^2$', 'rad/$s^2$'};
wind_axis = {'$W_x$ (fwd)', '$W_y$ (right)', '$W_z$ (down)'};
Val_max  = [max_dF;             max_dM];             % 6 x 20
Val_min  = [min_dF;             min_dM];             % 6 x 20
Wind_max = cat(1, max_dF_wind, max_dM_wind);         % 6 x 20 x [wx wy wz]
Wind_min = cat(1, min_dF_wind, min_dM_wind);         % 6 x 20 x [wx wy wz]
font_size = 14;
for k = 1:6
    name = comp_name{k};
    figure('Name', sprintf('Extreme %s and body-frame wind direction that produced it', name));

    subplot(4, 1, 1);
    hold on; grid on;
    plot(indices, Val_max(k, :), '-o', 'DisplayName', ['Max $' name '$']);
    plot(indices, Val_min(k, :), '-o', 'DisplayName', ['Min $' name '$']);
    legend('Location', 'best', 'Interpreter', 'latex');
    % xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
    ylabel(['$' name '$ [' comp_unit{k} ']'], 'Interpreter', 'latex', 'FontSize', font_size);

    % Body-axis components of the incoming-wind unit vector at the max / min point
    for ax = 1:3
        subplot(4, 1, ax + 1);
        hold on; grid on;
        plot(indices, 10.*Wind_max(k, :, ax), '-o', 'DisplayName', ['at Max $' name '$']);
        plot(indices, 10.*Wind_min(k, :, ax), '-o', 'DisplayName', ['at Min $' name '$']);
        ylim([-10.1, 10.1]);
        legend('Location', 'best', 'Interpreter', 'latex');
        if ax == 3
            xlabel('Index of Trim Conditions', 'Interpreter', 'latex', 'FontSize', font_size);
        end
        ylabel(['wind ' wind_axis{ax}], 'Interpreter', 'latex', 'FontSize', font_size);
    end
end

function [val, alpha_at, beta_at] = extreme_ab(slice, is_max, beta_pref, alpha_deg, beta_deg)
% Extreme value of a (N_alpha x N_beta) slice and the (alpha, beta) [deg] at which
% it occurs. On ties, beta_pref selects by sideslip magnitude |sin(beta)|:
%   'min' -> beta closest to 0 / 180 deg (least sideways wind)
%   'max' -> beta closest to +/-90 deg  (most sideways wind)
    if all(isnan(slice(:)))
        val = NaN; alpha_at = NaN; beta_at = NaN; return;
    end
    if is_max
        val = max(slice, [], 'all');
    else
        val = min(slice, [], 'all');
    end
    [rows, cols] = find(slice == val);       % rows = alpha index, cols = beta index
    sideslip = abs(sind(beta_deg(cols)));    % |sin(beta)|: 0 at 0/180 deg, 1 at +/-90 deg
    if strcmp(beta_pref, 'min')
        [~, sel] = min(sideslip);
    else
        [~, sel] = max(sideslip);
    end
    alpha_at = alpha_deg(rows(sel));
    beta_at  = beta_deg(cols(sel));
end

function w = wind_body_unit(alpha_deg, beta_deg)
% Body-frame unit direction of the incoming wind for aero incidence (alpha, beta) [deg],
% matching the sweep convention w_body = Vw * [cos a cos b; sin b; sin a cos b].
    a = deg2rad(alpha_deg);
    b = deg2rad(beta_deg);
    w = [cos(a) * cos(b); sin(b); sin(a) * cos(b)];
end
