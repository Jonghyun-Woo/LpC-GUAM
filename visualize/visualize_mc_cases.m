% Detailed view of ONE Monte-Carlo rollout from tests/test_mc_verify.m.
% Four figures: axis states, effector commands, the NED path (shifted so the
% t = T position sits at [0 0 -100]), and the airframe attitude at 0/0.5/1/1.5 s.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
CASE       = 'rphi';   % 'uw' | 'qtheta' (lon) | 'vp' | 'rphi' (lat)
UH_IDX     = 8;
SAMPLE     = 'auto';     % 'auto' picks the last sample to reach the target, or an index
MESH_TIMES = [0 0.5 1.0 1.5];
MESH_VIEW  = 'tail';     % lat only: 'tail' looks forward, 'front' looks aft
AXIS_EQUAL = false;      % path plot: true distorts nothing but flattens climb
% -------------------------------------------------------------------------

dataDir = fullfile(root, 'reachable_data', 'mc_verify_timestack');
figDir  = fullfile(dataDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

mc = load(fullfile(dataDir, sprintf('mc_%s_UH%d.mat', CASE, UH_IDX)));
ft2m = 0.3048;  r2d = 180 / pi;  radps2rpm = 60 / (2 * pi);

if strcmp(mc.axis, 'lon')
    plant_rows = [4; 6; 11; 8];   trim_rows = [1; 3; 5; 11];
    state_scale = [ft2m; ft2m; r2d; r2d];
    state_label = {'\Deltau [m/s]', '\Deltaw [m/s]', '\Deltaq [deg/s]', '\Delta\theta [deg]'};
else
    plant_rows = [5; 10; 12; 7];  trim_rows = [2; 4; 6; 10];
    state_scale = [ft2m; r2d; r2d; r2d];
    state_label = {'\Deltav [m/s]', '\Deltap [deg/s]', '\Deltar [deg/s]', '\Delta\phi [deg]'};
end

t = mc.t;
if strcmp(SAMPLE, 'auto')
    t_arrive = mc.t_reach;
    t_arrive(~mc.reach) = Inf;      % a sample that never arrives is the latest
    [~, SAMPLE] = max(t_arrive);
end

dev = zeros(numel(t), 4);
for d = 1:4
    dev(:, d) = (squeeze(mc.X(SAMPLE, :, plant_rows(d))) - mc.trim_state(trim_rows(d))) * state_scale(d);
end
euler   = squeeze(mc.X(SAMPLE, :, 7:9));
t_reach = mc.t_reach(SAMPLE);
tag     = sprintf('%s_UH%d_s%d', mc.case_name, mc.uh_idx, SAMPLE);
banner  = sprintf('%s  UH%d (u_{trim} = %.1f m/s)  sample %d', ...
                  mc.case_name, mc.uh_idx, mc.uh_vel * ft2m, SAMPLE);

cLine  = [0.00 0.45 0.74];
cBound = [0.85 0.10 0.10];

% --- states ---------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [40 40 1100 720]);
tl  = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for d = 1:4
    nexttile(tl); hold on; grid on; box on;
    plot(t, dev(:, d), '-', 'Color', cLine, 'LineWidth', 1.8);
    bound = mc.target_ub(d) * state_scale(d);
    yline( bound, '--', 'Color', cBound, 'LineWidth', 1.2);
    yline(-bound, '--', 'Color', cBound, 'LineWidth', 1.2);
    mark_reach(t_reach);
    xlim([t(1) t(end)]);  xlabel('t [s]');  ylabel(state_label{d});
end
title(tl, [banner '  --  states (dashed = target set)']);
exportgraphics(fig, fullfile(figDir, sprintf('case_%s_states.png', tag)), 'Resolution', 150);

% --- inputs ---------------------------------------------------------------
groups = effector_groups(mc.axis, radps2rpm, r2d);
fig = figure('Color', 'w', 'Position', [60 60 460 * numel(groups) 420]);
tl  = tiledlayout(fig, 1, numel(groups), 'Padding', 'compact', 'TileSpacing', 'compact');
for g = 1:numel(groups)
    nexttile(tl); hold on; grid on; box on;
    channels = groups(g).channels;
    plot(t, double(squeeze(mc.U(SAMPLE, :, channels))) * groups(g).scale, '-', 'LineWidth', 1.2);
    for c = channels
        yline(mc.input_lb(c) * groups(g).scale, ':', 'Color', cBound);
        yline(mc.input_ub(c) * groups(g).scale, ':', 'Color', cBound);
    end
    if numel(channels) <= 3
        legend(groups(g).names, 'Location', 'best', 'Box', 'off');
    end
    mark_reach(t_reach);
    xlim([t(1) t(end)]);  xlabel('t [s]');  ylabel(groups(g).label);
    title(groups(g).title);
end
title(tl, [banner '  --  effector commands (dotted = input box)']);
exportgraphics(fig, fullfile(figDir, sprintf('case_%s_inputs.png', tag)), 'Resolution', 150);

% --- NED path -------------------------------------------------------------
position = squeeze(mc.X(SAMPLE, :, 1:3));
position = position - position(end, :) + [0 0 -100];
mark_idx = round(MESH_TIMES / mc.dt) + 1;

fig = figure('Color', 'w', 'Position', [80 80 860 720]);
ax  = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
plot3(ax, position(:, 1), position(:, 2), position(:, 3), '-', 'Color', cLine, 'LineWidth', 1.8);
plot3(ax, position(mark_idx, 1), position(mark_idx, 2), position(mark_idx, 3), 'ko', ...
      'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 1.4);
for k = 1:numel(mark_idx)
    text(ax, position(mark_idx(k), 1), position(mark_idx(k), 2), position(mark_idx(k), 3), ...
         sprintf('  %.1f s', MESH_TIMES(k)), 'FontWeight', 'bold');
end
if ~isnan(t_reach)
    i_reach = round(t_reach / mc.dt) + 1;
    plot3(ax, position(i_reach, 1), position(i_reach, 2), position(i_reach, 3), 'p', ...
          'Color', [0.10 0.60 0.25], 'MarkerFaceColor', [0.10 0.60 0.25], 'MarkerSize', 13);
end
if AXIS_EQUAL, axis(ax, 'equal'); end
set(ax, 'ZDir', 'reverse');  view(ax, 135, 22);
xlabel(ax, 'north [ft]'); ylabel(ax, 'east [ft]'); zlabel(ax, 'down [ft]');
title(ax, sprintf('%s  --  NED path (t = %.1f s pinned to [0 0 -100])', banner, t(end)));
exportgraphics(fig, fullfile(figDir, sprintf('case_%s_path.png', tag)), 'Resolution', 150);

% --- attitude meshes ------------------------------------------------------
parts = load_lpc_geometry();
parts = parts(~contains({parts.name}, 'disk', 'IgnoreCase', true));
colours = part_colours(parts);
span = max(vecnorm(vertcat(parts.vertices), 2, 2));

fig = figure('Color', 'w', 'Position', [100 100 1100 900]);
tl  = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for k = 1:numel(mark_idx)
    step = mark_idx(k);
    R_body_to_ned = RSLQR.rotm_i2b(euler(step, 1), euler(step, 2), euler(step, 3)).';

    ax = nexttile(tl); hold(ax, 'on');
    for m = 1:numel(parts)
        patch(ax, 'Vertices', parts(m).vertices * R_body_to_ned.', 'Faces', parts(m).faces, ...
              'FaceColor', colours(m, :), 'EdgeColor', 'none', ...
              'FaceLighting', 'gouraud', 'AmbientStrength', 0.45);
    end
    axis(ax, 'equal');
    xlim(ax, [-span span]); ylim(ax, [-span span]); zlim(ax, [-span span]);
    set(ax, 'ZDir', 'reverse');
    apply_camera(ax, mc.axis, MESH_VIEW);
    camlight(ax, 'headlight'); camlight(ax, 'left');
    title(ax, sprintf('t = %.1f s   \\phi = %+.1f\\circ  \\theta = %+.1f\\circ  \\psi = %+.1f\\circ', ...
                      MESH_TIMES(k), rad2deg(euler(step, 1)), rad2deg(euler(step, 2)), ...
                      rad2deg(euler(step, 3))));
end
title(tl, [banner '  --  attitude, ' camera_name(mc.axis, MESH_VIEW)]);
exportgraphics(fig, fullfile(figDir, sprintf('case_%s_mesh.png', tag)), 'Resolution', 150);

fprintf('saved case_%s_{states,inputs,path,mesh}.png to %s\n', tag, figDir);


function groups = effector_groups(axis_name, radps2rpm, r2d)
    lift = struct('title', 'lift rotors 1-8', 'channels', 1:8, ...
                  'scale', radps2rpm, 'label', '\Delta\Omega [RPM]', 'names', {{}});
    if strcmp(axis_name, 'lon')
        groups = [lift, ...
                  struct('title', 'pusher', 'channels', 9, 'scale', radps2rpm, ...
                         'label', '\Delta\Omega [RPM]', 'names', {{'pusher'}}), ...
                  struct('title', 'control surfaces', 'channels', [10 11], 'scale', r2d, ...
                         'label', '\Delta\delta [deg]', 'names', {{'elevator', 'flap'}})];
    else
        groups = [lift, ...
                  struct('title', 'control surfaces', 'channels', [9 10], 'scale', r2d, ...
                         'label', '\Delta\delta [deg]', 'names', {{'aileron', 'rudder'}})];
    end
end

function apply_camera(ax, axis_name, mesh_view)
    if strcmp(axis_name, 'lon')
        view(ax, [0 -1 0]);          % side view: theta in the screen plane
    elseif strcmp(mesh_view, 'front')
        view(ax, [-1 0 0]);          % ahead of the aircraft, looking aft
    else
        view(ax, [1 0 0]);           % behind the aircraft, looking forward
    end
end

function name = camera_name(axis_name, mesh_view)
    if strcmp(axis_name, 'lon')
        name = 'side view';
    elseif strcmp(mesh_view, 'front')
        name = 'front view';
    else
        name = 'tail view';
    end
end

function mark_reach(t_reach)
    if ~isnan(t_reach)
        xline(t_reach, '-', 'Color', [0.10 0.60 0.25], 'LineWidth', 1.4);
    end
end

function colours = part_colours(parts)
    colours = repmat([0.72 0.74 0.78], numel(parts), 1);
    for k = 1:numel(parts)
        name = parts(k).name;
        if     contains(name, 'blades'),   colours(k, :) = [0.25 0.27 0.31];
        elseif contains(name, 'Wing'),     colours(k, :) = [0.42 0.58 0.78];
        elseif contains(name, 'Tail'),     colours(k, :) = [0.52 0.66 0.84];
        elseif contains(name, 'Fuselage'), colours(k, :) = [0.86 0.87 0.89];
        end
    end
end
