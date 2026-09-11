function parts = load_lpc_geometry(stp_path, samples_per_segment)
% Parses the CADDEE LPC STEP file into per-component patch surfaces expressed
% in the GUAM body frame (x forward, y right, z down), centred on the model
% bounding box. Result is cached next to this file and reused unless the STEP
% file is newer.
%
% parts : struct array with fields name, vertices (nv x 3), faces (nf x 4)

    if nargin < 1 || isempty(stp_path)
        stp_path = fullfile('C:', 'WooJongHyun', 'Git', 'CADDEE_alpha', ...
                            'examples', 'test_geometries', 'LPC_final_custom_blades.stp');
    end
    if nargin < 2 || isempty(samples_per_segment), samples_per_segment = 4; end

    geometry_version = 2;   % bump whenever the frame or tessellation changes
    cache_path = fullfile(fileparts(mfilename('fullpath')), 'lpc_geometry_cache.mat');
    stp_info   = dir(stp_path);
    if isempty(stp_info)
        error('load_lpc_geometry:missingStep', 'STEP file not found: %s', stp_path);
    end
    if isfile(cache_path)
        cached = load(cache_path);
        if isfield(cached, 'geometry_version') && ...
           cached.geometry_version == geometry_version && ...
           cached.stp_datenum == stp_info.datenum && ...
           cached.samples_per_segment == samples_per_segment
            parts = cached.parts;
            return;
        end
    end

    lines = read_lines(stp_path);
    [point_xyz, row_of_id] = parse_points(lines);
    surfaces   = parse_surfaces(lines);
    comp_of_id = parse_components(lines);

    grids = cell(numel(surfaces), 1);
    names = cell(numel(surfaces), 1);
    for k = 1:numel(surfaces)
        control_net = reshape(point_xyz(row_of_id(surfaces(k).point_ids(:)), :), ...
                              [size(surfaces(k).point_ids), 3]);
        grids{k} = eval_bezier_patch(control_net, surfaces(k).degree_u, ...
                                     surfaces(k).degree_v, samples_per_segment);
        names{k} = comp_of_id(surfaces(k).id);
    end

    parts = group_into_parts(grids, names);
    parts = to_body_frame(parts);

    stp_datenum = stp_info.datenum; %#ok<NASGU>
    save(cache_path, 'parts', 'stp_datenum', 'samples_per_segment', 'geometry_version');
end


function lines = read_lines(path)
    fid = fopen(path, 'r');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    scanned = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = scanned{1};
end

function [point_xyz, row_of_id] = parse_points(lines)
    text_ = strjoin(lines(contains(lines, '=CARTESIAN_POINT(')), '');
    text_ = regexprep(text_, '=CARTESIAN_POINT\(''[^'']*'',\(', ' ');
    text_ = strrep(strrep(strrep(text_, '#', ' '), ',', ' '), '));', ' ');
    id_xyz = reshape(sscanf(text_, '%f'), 4, []).';

    point_xyz = id_xyz(:, 2:4);
    row_of_id = zeros(max(id_xyz(:, 1)), 1);
    row_of_id(id_xyz(:, 1)) = 1:size(id_xyz, 1);
end

function surfaces = parse_surfaces(lines)
    tokens = regexp(lines(contains(lines, '=B_SPLINE_SURFACE_WITH_KNOTS(')), ...
                    '^#(\d+)=B_SPLINE_SURFACE_WITH_KNOTS\(''(.*?)'',(\d+),(\d+),\(\((.*?)\)\),\.', ...
                    'tokens', 'once');
    surfaces = struct('id', {}, 'degree_u', {}, 'degree_v', {}, 'point_ids', {});
    for k = 1:numel(tokens)
        net_rows = strsplit(tokens{k}{5}, '),(');
        ids = cellfun(@(r) id_list(r).', net_rows, 'UniformOutput', false);
        surfaces(k) = struct('id', str2double(tokens{k}{1}), ...
                             'degree_u', str2double(tokens{k}{3}), ...
                             'degree_v', str2double(tokens{k}{4}), ...
                             'point_ids', vertcat(ids{:}));
    end
end

function comp_of_id = parse_components(lines)
    tokens = regexp(lines(contains(lines, '=GEOMETRIC_SET(')), ...
                    '^#(\d+)=GEOMETRIC_SET\(''(.*?)'',\((.*?)\)\)', 'tokens', 'once');
    comp_of_id = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for k = 1:numel(tokens)
        % Set names read "<code>, <component>, <index>".
        fields = strsplit(tokens{k}{2}, ', ');
        name   = fields{min(2, numel(fields))};
        ids    = id_list(tokens{k}{3});
        for j = 1:numel(ids)
            comp_of_id(ids(j)) = name;
        end
    end
end

function ids = id_list(str)
    ids = sscanf(strrep(strrep(str, '#', ' '), ',', ' '), '%f');
end

function grid_xyz = eval_bezier_patch(control_net, degree_u, degree_v, ns)
    % The STEP uses PIECEWISE_BEZIER_KNOTS throughout, so the control net is
    % already in Bezier form and each segment evaluates from Bernstein bases.
    [n_u, n_v, ~] = size(control_net);
    n_seg_u = (n_u - 1) / degree_u;
    n_seg_v = (n_v - 1) / degree_v;
    basis_u = bernstein(degree_u, linspace(0, 1, ns));
    basis_v = bernstein(degree_v, linspace(0, 1, ns));

    grid_xyz = zeros(n_seg_u * ns, n_seg_v * ns, 3);
    for i = 1:n_seg_u
        for j = 1:n_seg_v
            block = control_net(1 + degree_u*(i-1) : 1 + degree_u*i, ...
                                1 + degree_v*(j-1) : 1 + degree_v*j, :);
            for c = 1:3
                grid_xyz((i-1)*ns + (1:ns), (j-1)*ns + (1:ns), c) = ...
                    basis_u * block(:, :, c) * basis_v.';
            end
        end
    end
end

function basis = bernstein(degree, t)
    t = t(:);
    basis = zeros(numel(t), degree + 1);
    for i = 0:degree
        basis(:, i+1) = nchoosek(degree, i) * t.^i .* (1 - t).^(degree - i);
    end
end

function parts = group_into_parts(grids, names)
    [unique_names, ~, group_of] = unique(names, 'stable');
    parts = struct('name', {}, 'vertices', {}, 'faces', {});
    for g = 1:numel(unique_names)
        vertices = [];
        faces    = [];
        for k = find(group_of(:).' == g)
            grid_xyz = grids{k};
            [n_row, n_col, ~] = size(grid_xyz);
            faces    = [faces; quad_faces(n_row, n_col) + size(vertices, 1)]; %#ok<AGROW>
            vertices = [vertices; reshape(grid_xyz, [], 3)];                  %#ok<AGROW>
        end
        parts(g) = struct('name', unique_names{g}, 'vertices', vertices, 'faces', faces);
    end
end

function faces = quad_faces(n_row, n_col)
    [row, col] = ndgrid(1:n_row-1, 1:n_col-1);
    node = @(i, j) (j - 1) * n_row + i;
    faces = [node(row(:), col(:)), node(row(:)+1, col(:)), ...
             node(row(:)+1, col(:)+1), node(row(:), col(:)+1)];
end

function parts = to_body_frame(parts)
    % CAD frame is x aft, y right, z up and already shares the GUAM body origin
    % (pusher hub sits at x_cad 31.94 <-> Prop_location(1,9) = -31.94), so only
    % the x/z signs flip. Vertices are then referred to the CG, which is what
    % the airframe rotates about.
    cg = VehicleConfig.cm_b(:).';
    for k = 1:numel(parts)
        v = parts(k).vertices;
        parts(k).vertices = [-v(:, 1), v(:, 2), -v(:, 3)] - cg;
    end
end
