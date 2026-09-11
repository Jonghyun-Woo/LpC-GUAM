% Animates the LpC attitude from a Monte-Carlo rollout written by
% tests/test_mc_verify.m. The airframe is rotated in place by the sampled
% (phi, theta, psi); an animated GIF is written next to the source data.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
MC_FILE    = fullfile(root, 'reachable_data', 'mc_verify_timestack', 'mc_qtheta_UH8.mat');
SAMPLE     = 'auto';     % 'auto' picks the largest attitude excursion, or an index
FRAME_STEP = 3;          % rollout steps per animation frame
FPS        = 20;
GIF_PATH   = fullfile(root, 'reachable_data', 'mc_verify_timestack', 'figures', ...
                      'attitude_qtheta_UH8.gif');
% -------------------------------------------------------------------------

parts = load_lpc_geometry();
mc    = load(MC_FILE);

euler = permute(mc.X(:, :, 7:9), [2 3 1]);          % step x [phi theta psi] x sample
trim_euler = mc.trim_state(10:12).';
if strcmp(SAMPLE, 'auto')
    excursion = squeeze(max(vecnorm(euler - trim_euler, 2, 2), [], 1));
    [~, SAMPLE] = max(excursion);
end
euler = euler(:, :, SAMPLE);

if ~exist(fileparts(GIF_PATH), 'dir'), mkdir(fileparts(GIF_PATH)); end

% Actuator disks are analysis surfaces, not airframe; they would hide the blades.
keep  = ~contains({parts.name}, 'disk', 'IgnoreCase', true);
parts = parts(keep);
face_colour = part_colours(parts);

span = max(vecnorm(vertcat(parts.vertices), 2, 2));
fig  = figure('Color', 'w', 'Position', [80 80 900 760]);
ax   = axes(fig); hold(ax, 'on');
handles = gobjects(numel(parts), 1);
for k = 1:numel(parts)
    handles(k) = patch(ax, 'Vertices', parts(k).vertices, 'Faces', parts(k).faces, ...
                       'FaceColor', face_colour(k, :), 'EdgeColor', 'none', ...
                       'FaceLighting', 'gouraud', 'AmbientStrength', 0.45);
end
axis(ax, 'equal');
% box(ax, 'on');
% grid(ax, 'on');
xlim(ax, [-span span]); ylim(ax, [-span span]); zlim(ax, [-span span]);
set(ax, 'ZDir', 'reverse');
xlabel(ax, 'north [ft]'); ylabel(ax, 'east [ft]'); zlabel(ax, 'down [ft]');
view(ax, 135, 20); camlight(ax, 'headlight'); camlight(ax, 'left');

frames = 1:FRAME_STEP:size(euler, 1);
for i_frame = 1:numel(frames)
    step = frames(i_frame);
    phi = euler(step, 1);  theta = euler(step, 2);  psi = euler(step, 3);
    R_body_to_ned = RSLQR.rotm_i2b(phi, theta, psi).';

    for k = 1:numel(parts)
        handles(k).Vertices = parts(k).vertices * R_body_to_ned.';
    end
    title(ax, sprintf(['%s  UH%d  sample %d\n' ...
                       't = %.2f s   \\phi = %+.1f\\circ   \\theta = %+.1f\\circ   \\psi = %+.1f\\circ'], ...
                      mc.case_name, mc.uh_idx, SAMPLE, mc.t(step), ...
                      rad2deg(phi), rad2deg(theta), rad2deg(psi)));
    drawnow;

    [gif_frame, colour_map] = rgb2ind(frame2im(getframe(fig)), 256);
    if i_frame == 1
        imwrite(gif_frame, colour_map, GIF_PATH, 'gif', 'LoopCount', Inf, 'DelayTime', 1/FPS);
    else
        imwrite(gif_frame, colour_map, GIF_PATH, 'gif', 'WriteMode', 'append', 'DelayTime', 1/FPS);
    end
end
fprintf('saved %s (%d frames)\n', GIF_PATH, numel(frames));


function colours = part_colours(parts)
    colours = repmat([0.72 0.74 0.78], numel(parts), 1);
    for k = 1:numel(parts)
        name = parts(k).name;
        if     contains(name, 'blades'),  colours(k, :) = [0.25 0.27 0.31];
        elseif contains(name, 'Wing'),    colours(k, :) = [0.42 0.58 0.78];
        elseif contains(name, 'Tail'),    colours(k, :) = [0.52 0.66 0.84];
        elseif contains(name, 'Fuselage'),colours(k, :) = [0.86 0.87 0.89];
        end
    end
end
