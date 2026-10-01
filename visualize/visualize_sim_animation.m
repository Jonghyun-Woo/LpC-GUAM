% Animates a saved closed-loop run (run_transition_sim results). Two separate
% GIFs are written so each frame stays cheap to render: one for the inertial
% position (full trajectory drawn once, current step marked with a moving
% scatter) and one for the airframe attitude (mesh only, no axes decoration).

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
RESULTS_FILE = fullfile(root, 'reachability_data', 'transition_results.mat');
MODE         = 'blend';   % 'blend' or 'off'
FRAME_STEP   = 10;        % rollout steps per animation frame
FPS          = 30;
MESH_KEEP    = 0.08;      % fraction of airframe faces to keep (reducepatch)
OUT_DIR      = fullfile(root, 'visualize', 'figures');
% -------------------------------------------------------------------------

results = load_sim_results(RESULTS_FILE);
data    = results.(MODE).data;
if ~isfolder(OUT_DIR), mkdir(OUT_DIR); end

animate_position(data, FRAME_STEP, FPS, fullfile(OUT_DIR, 'sim_position.gif'));
animate_attitude(data, FRAME_STEP, FPS, MESH_KEEP, fullfile(OUT_DIR, 'sim_attitude.gif'));


function animate_position(data, frame_step, fps, gif_path)
    ft2m = 0.3048;
    pos  = ft2m * data.state(:, 1:3);   % [north east down] in m
    t    = data.t;

    fig = figure('Color', 'w', 'Position', [80 80 900 760]);
    ax  = axes(fig); hold(ax, 'on'); grid(ax, 'on');
    trj = plot3(ax, pos(:, 1), pos(:, 2), pos(:, 3), 'r-', 'LineWidth', 1.2);
    marker = scatter3(ax, pos(1, 1), pos(1, 2), pos(1, 3), 70, [0.85 0.2 0.2], 'o', 'filled');
    set(ax, 'ZDir', 'reverse');
    xlabel(ax, 'north [m]'); ylabel(ax, 'east [m]'); zlabel(ax, 'down [m]');
    view(ax, -25, 7); axis(ax, 'equal'); 
    xlim([0, 1600]); ylim([-400, 400]); zlim([-400, 200]);

    frames = 1:frame_step:size(pos, 1);
    for i_frame = 1:numel(frames)
        step = frames(i_frame);
        set(trj, 'XData', pos(1:step, 1), 'YData', pos(1:step, 2), 'ZData', pos(1:step, 3));
        set(marker, 'XData', pos(step, 1), 'YData', pos(step, 2), 'ZData', pos(step, 3));
        drawnow;
        write_gif(fig, gif_path, i_frame == 1, fps);
    end
    fprintf('saved %s (%d frames)\n', gif_path, numel(frames));
end


function animate_attitude(data, frame_step, fps, mesh_keep, gif_path)
    parts = load_lpc_geometry();
    keep  = ~contains({parts.name}, 'disk', 'IgnoreCase', true);
    parts = parts(keep);
    if mesh_keep < 1
        for k = 1:numel(parts)
            nfv = reducepatch(parts(k).faces, parts(k).vertices, mesh_keep);
            parts(k).faces = nfv.faces;  parts(k).vertices = nfv.vertices;
        end
    end
    face_colour = part_colours(parts);
    euler = data.state(:, 7:9);   % [phi theta psi]

    span = max(vecnorm(vertcat(parts.vertices), 2, 2));
    fig  = figure('Color', 'w', 'Position', [80 80 900 760]);
    ax   = axes(fig); hold(ax, 'on');
    handles = gobjects(numel(parts), 1);
    for k = 1:numel(parts)
        handles(k) = patch(ax, 'Vertices', parts(k).vertices, 'Faces', parts(k).faces, ...
                           'FaceColor', face_colour(k, :), 'EdgeColor', 'none', ...
                           'FaceLighting', 'gouraud', 'AmbientStrength', 0.45);
    end
    axis(ax, 'off');
    xlim(ax, [-span span]); ylim(ax, [-span span]); zlim(ax, [-span span]);
    set(ax, 'ZDir', 'reverse');
    view(ax, 0, 0); camlight(ax, 'headlight'); camlight(ax, 'left');

    frames = 1:frame_step:size(euler, 1);
    for i_frame = 1:numel(frames)
        step = frames(i_frame);
        R_body_to_ned = RSLQR.rotm_i2b(euler(step, 1), euler(step, 2), euler(step, 3)).';
        for k = 1:numel(parts)
            handles(k).Vertices = parts(k).vertices * R_body_to_ned.';
        end
        drawnow;
        write_gif(fig, gif_path, i_frame == 1, fps);
    end
    fprintf('saved %s (%d frames)\n', gif_path, numel(frames));
end


function write_gif(fig, gif_path, is_first, fps)
    [gif_frame, colour_map] = rgb2ind(frame2im(getframe(fig)), 256);
    if is_first
        imwrite(gif_frame, colour_map, gif_path, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
    else
        imwrite(gif_frame, colour_map, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end


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
