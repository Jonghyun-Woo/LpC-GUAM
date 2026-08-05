clear; close all; clc;

cfg = Config();
GUAM = LpC_GUAM(cfg);

load('controller/trim_table_Poly_ConcatVer4p0.mat');

max_dF = nan(3, 20);
min_dF = nan(3, 20);
max_dM = nan(3, 20);
min_dM = nan(3, 20);
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
    Na = numel(alpha_deg);  Nb = numel(beta_deg);
    dVel_b = zeros(3, Na, Nb);   % disturbance acceleration           in body frame
    dOmg_b = zeros(3, Na, Nb);   % disturbance angular acceleration   in body frame
    valid = true(Na, Nb);      % aero-database validity of each sweep point
    
    for ia = 1:Na
        for ib = 1:Nb
            a = deg2rad(alpha_deg(ia));
            b = deg2rad(beta_deg(ib));
            % 10 m/s wind, direction (a,b) in body axes (aero convention)
            w_body = Vw * [cos(a)*cos(b); sin(b); sin(a)*cos(b)];
            GUAM.environment.wind_ned = R_b2i * w_body;   % inject as body-frame wind
    
            ddx = GUAM.state_derivative(X0, engine0, surface0) - dx0;
            dVel_b(:, ia, ib) = ddx(4:6);     % dF <- acc
            dOmg_b(:, ia, ib) = ddx(10:12);   % dM <- omega_dot
    
            % Re-run the aero model to read back its validity flags.
            v_air = X0(4:6) - w_body;
            [~, ~, Validity] = run_LpC_aero([v_air; X0(10:12)], engine0, surface0, ...
                                            rho, a_snd, GUAM.units, GUAM.vehicleConfig);
            valid(ia, ib) = ~any(Validity);
        end
    end
    GUAM.environment.wind_ned = zeros(3, 1);   % restore calm air
    
    % Drop the invalid aerodynamic points.
    dVel_b(:, ~valid) = NaN;
    dOmg_b(:, ~valid) = NaN;
    fprintf('aero-database validity: %d / %d sweep points excluded (%.1f%%)\n', ...
            nnz(~valid), numel(valid), 100 * nnz(~valid) / numel(valid));
    max_dF(:, idx) = max(dVel_b, [], [2, 3]);
    min_dF(:, idx) = min(dVel_b, [], [2, 3]);
    max_dM(:, idx) = max(dOmg_b, [], [2, 3]);
    min_dM(:, idx) = min(dOmg_b, [], [2, 3]);
end

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
