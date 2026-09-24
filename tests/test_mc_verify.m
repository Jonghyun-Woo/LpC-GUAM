% Monte-Carlo reach test of the time-varying BRT stacks on the nonlinear plant.
% Cases: uw / qtheta (lon), vp / rphi (lat). Two states are sampled inside the
% full-horizon tube, the other two held at FIX; the rollout uses the BRT-optimal
% bang-bang input from the time-varying gradV. Success = target box within T.
% Per-case state histories are written for tests/viz_mc_verify.m.

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(root));

% ------------------------------- settings --------------------------------
DATASET  = 'guam_timestack';
UH_LIST  = [1 4 8 12 16 20];
WH_IDX   = 2;
N_SAMPLE = 100;
DT       = 0.01;
SEED     = 20260906;

FIX.lon = [0; 0; 0; 0];   % [u, w, q, theta] held values (deviation)
FIX.lat = [0; 0; 0; 0];   % [v, p, r, phi]
% -------------------------------------------------------------------------

runDir = fullfile(root, 'reachability_data', DATASET, 'BRT');
outDir = fullfile(root, 'reachability_data', 'mc_verify_timestack');
if ~exist(outDir, 'dir'), mkdir(outDir); end

yml       = read_yml(fullfile(runDir, 'guam_analysis_config.yml'));
T_HORIZON = yml.hj_analysis_config.time;
N_STEP    = round(T_HORIZON / DT);
axisCfg   = axis_config(yml);
axisCfg.lon.held_dev = FIX.lon;
axisCfg.lat.held_dev = FIX.lat;

cases = struct('name',      {'uw', 'qtheta', 'vp',  'rphi'}, ...
               'axis',      {'lon', 'lon',   'lat', 'lat'}, ...
               'free_dims', {[1 2], [3 4],   [1 2], [3 4]});

trimTable = load(fullfile(root, 'controller', 'trim_table_Poly_ConcatVer4p0.mat'), ...
                 'XU0_interp', 'UH', 'Bp_lon_interp', 'Bp_lat_interp');

fprintf('dataset  : %s\n', DATASET);
fprintf('horizon  : %.2f s (%d steps at dt = %.3f s)\n', T_HORIZON, N_STEP, DT);
fprintf('samples  : %d per (case, UH),  UH = %s\n\n', N_SAMPLE, mat2str(UH_LIST));

n_uh    = numel(UH_LIST);
uh_rows = cell(n_uh, 1);

% parfor over UH, not over samples: the lat stack is ~665 MB, so broadcasting
% it per case would cost more than the rollouts. Workers load their own.
parfor i_uh = 1:n_uh
    uh_idx = UH_LIST(i_uh);
    XU0    = trimTable.XU0_interp(:, uh_idx, WH_IDX);
    GUAM   = LpC_GUAM(Config('althold', struct('dt', DT)));
    case_rows = struct('case', {}, 'uh', {}, 'reach', {}, 't_reach', {}, 'left_grid', {});

    for i_axis = 1:2
        axis_name = 'lon';
        if i_axis == 2, axis_name = 'lat'; end
        gridCfg = axisCfg.(axis_name);

        V_stack = load_stack(fullfile(runDir, sprintf('BRT_%s_MAT', upper(axis_name)), ...
                             sprintf('GUAM_%s_BRT_UH%d_WH%d_stack.mat', ...
                                     upper(axis_name), uh_idx, WH_IDX)));
        gridCfg  = add_stack_indexing(gridCfg, size(V_stack));
        dt_slice = T_HORIZON / (gridCfg.n_slice - 1);

        sched = struct('U0', XU0(13:25), 'X0', XU0(1:12), ...
                       'axis', axis_name, 'plant_rows', gridCfg.plant_rows);
        if i_axis == 1
            sched.Bp = trimTable.Bp_lon_interp(:, :, uh_idx, WH_IDX);
        else
            sched.Bp = trimTable.Bp_lat_interp(:, :, uh_idx, WH_IDX);
        end
        u_trim         = sched.U0(gridCfg.u0_idx);
        sched.input_lb = max(gridCfg.phys_lb, u_trim + gridCfg.delta_lb) - u_trim;
        sched.input_ub = min(gridCfg.phys_ub, u_trim + gridCfg.delta_ub) - u_trim;

        for i_case = 1:numel(cases)
            if ~strcmp(cases(i_case).axis, axis_name), continue; end
            free_dims = cases(i_case).free_dims;

            rng(SEED + 1000 * uh_idx + i_case);
            x0_dev   = sample_inside_tube(V_stack, gridCfg, free_dims, ...
                                          gridCfg.held_dev, N_SAMPLE);
            n_sample = size(x0_dev, 1);

            X_hist = nan(n_sample, N_STEP + 1, 12);
            U_hist = nan(n_sample, N_STEP + 1, numel(sched.input_lb), 'single');
            V_hist = nan(n_sample, N_STEP + 1, 'single');
            reach  = false(n_sample, 1);
            t_reach   = nan(n_sample, 1);
            left_grid = false(n_sample, 1);

            for i_sample = 1:n_sample
                [X_hist(i_sample, :, :), U_hist(i_sample, :, :), V_hist(i_sample, :), ...
                 reach(i_sample), t_reach(i_sample), left_grid(i_sample)] = ...
                    rollout(GUAM, V_stack, gridCfg, sched, x0_dev(i_sample, :), ...
                            DT, N_STEP, dt_slice);
            end

            result = struct( ...
                'case_name', cases(i_case).name, 'axis', axis_name, ...
                'uh_idx', uh_idx, 'wh_idx', WH_IDX, 'uh_vel', trimTable.UH(uh_idx), ...
                'dt', DT, 'T', T_HORIZON, 't', (0:N_STEP) * DT, ...
                'free_dims', free_dims, 'x0_dev', x0_dev, ...
                'X', X_hist, 'U', U_hist, 'V', V_hist, ...
                'reach', reach, 't_reach', t_reach, 'left_grid', left_grid, ...
                'brt_contour', {tube_contour(V_stack, gridCfg, free_dims)}, ...
                'trim_state', sched.X0, 'trim_input', sched.U0, ...
                'input_lb', sched.input_lb, 'input_ub', sched.input_ub, ...
                'grid_min', gridCfg.grid_min, 'grid_max', gridCfg.grid_max, ...
                'target_ub', gridCfg.target_ub);
            save_result(fullfile(outDir, sprintf('mc_%s_UH%d.mat', ...
                                 cases(i_case).name, uh_idx)), result);

            case_rows(end + 1) = struct('case', cases(i_case).name, 'uh', uh_idx, ...
                                        'reach', reach, 't_reach', t_reach, ...
                                        'left_grid', left_grid);
        end
        V_stack = [];
    end
    uh_rows{i_uh} = case_rows;
end

summary = [uh_rows{:}];
fprintf('%-7s %-4s %6s %8s %10s %10s\n', ...
        'case', 'UH', 'n', 'reach%', 'mean t_r', 'leftgrid%');
for i_case = 1:numel(cases)
    rows = summary(strcmp({summary.case}, cases(i_case).name));
    [~, order] = sort([rows.uh]);
    for i_row = order
        row = rows(i_row);
        fprintf('%-7s %-4d %6d %8.1f %10.3f %10.1f\n', row.case, row.uh, ...
                numel(row.reach), 100 * mean(row.reach), ...
                mean(row.t_reach(row.reach)), 100 * mean(row.left_grid));
    end
end

save(fullfile(outDir, 'mc_summary.mat'), 'summary', 'DATASET', 'UH_LIST', ...
     'WH_IDX', 'N_SAMPLE', 'DT', 'T_HORIZON', 'FIX', 'SEED');
fprintf('\nsaved %s\n', fullfile(outDir, 'mc_summary.mat'));

%% Local functions
function [X_hist, U_hist, V_hist, reach, t_reach, left_grid] = ...
         rollout(GUAM, V_stack, gridCfg, sched, x0_dev, dt, n_step, dt_slice)
    state = [0; 0; -100; sched.X0(1:3); sched.X0(10:12); sched.X0(4:6)];
    trim_ref = state(sched.plant_rows);
    state(sched.plant_rows) = trim_ref + x0_dev(:);

    X_hist = nan(1, n_step + 1, 12);
    U_hist = nan(1, n_step + 1, numel(sched.input_lb), 'single');
    V_hist = nan(1, n_step + 1, 'single');
    reach  = false;  t_reach = NaN;  left_grid = false;

    for step = 0:n_step
        dev = state(sched.plant_rows) - trim_ref;
        X_hist(1, step + 1, :) = state;

        % Stack is stored in backward time (last slice = target set), so
        % elapsed time indexes forward from slice 1.
        slice = min(gridCfg.n_slice, max(1, 1 + round(step * dt / dt_slice)));
        [V, gradV] = value_and_gradient(V_stack, gridCfg, slice, dev);
        V_hist(1, step + 1) = V;

        if ~reach && all(abs(dev) <= gridCfg.target_ub)
            reach = true;  t_reach = step * dt;
        end
        if any(dev < gridCfg.grid_min | dev > gridCfg.grid_max), left_grid = true; end

        switch_fn = sched.Bp' * gradV;
        u = (switch_fn >= 0) .* sched.input_lb + (switch_fn < 0) .* sched.input_ub;
        U_hist(1, step + 1, :) = u;
        [engine_cmd, surface_cmd] = actuator_cmd(sched, u);

        if step == n_step, break; end
        k1 = GUAM.state_derivative(state,             engine_cmd, surface_cmd);
        k2 = GUAM.state_derivative(state + 0.5*dt*k1, engine_cmd, surface_cmd);
        k3 = GUAM.state_derivative(state + 0.5*dt*k2, engine_cmd, surface_cmd);
        k4 = GUAM.state_derivative(state +     dt*k3, engine_cmd, surface_cmd);
        state = state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

function [engine_cmd, surface_cmd] = actuator_cmd(sched, u)
    U0 = sched.U0;
    if strcmp(sched.axis, 'lon')  % u = [lift1..8, pusher, elevator, flap]
        engine_cmd = U0(5:13) + [u(1:8); u(9)];
        flap = U0(1) + u(11);  ail = U0(2);
        ele  = U0(3) + u(10);  rud = U0(4);
    else                          % u = [lift1..8, aileron, rudder]
        engine_cmd = U0(5:13) + [u(1:8); 0];
        flap = U0(1);          ail = U0(2) + u(9);
        ele  = U0(3);          rud = U0(4) + u(10);
    end
    surface_cmd = [flap - ail; flap + ail; ele; ele; rud];
end

function [V, gradV] = value_and_gradient(V_stack, gridCfg, slice, dev)
    % Multilinear interpolation straight on the 5-D stack (no per-slice copy).
    % Nine query points give V and its central-difference gradient in one pass.
    query = repmat(dev(:)', 9, 1);
    for d = 1:4
        query(2*d,     d) = dev(d) + gridCfg.grid_step(d);
        query(2*d + 1, d) = dev(d) - gridCfg.grid_step(d);
    end
    f_idx = (query - gridCfg.grid_min') ./ gridCfg.grid_step' + 1;
    f_idx = min(max(f_idx, 1), gridCfg.grid_num');
    i_lo  = min(floor(f_idx), gridCfg.grid_num' - 1);
    frac  = f_idx - i_lo;

    vals = zeros(9, 1);
    for i_corner = 1:16
        corner = gridCfg.corners(i_corner, :);
        weight = prod(corner .* frac + (1 - corner) .* (1 - frac), 2);
        vals   = vals + weight .* ...
                 double(V_stack(slice + (i_lo + corner - 1) * gridCfg.stride));
    end
    V     = vals(1);
    gradV = (vals(2:2:8) - vals(3:2:9)) ./ (2 * gridCfg.grid_step);
end

function x0_dev = sample_inside_tube(V_stack, gridCfg, free_dims, held_dev, n_want)
    x0_dev = zeros(0, 4);
    for attempt = 1:40 %#ok<NASGU>
        n_cand = 20 * n_want;
        cand   = repmat(held_dev(:)', n_cand, 1);
        for d = free_dims
            cand(:, d) = gridCfg.grid_min(d) + ...
                         (gridCfg.grid_max(d) - gridCfg.grid_min(d)) * rand(n_cand, 1);
        end
        inside = false(n_cand, 1);
        for i_cand = 1:n_cand
            inside(i_cand) = value_and_gradient(V_stack, gridCfg, 1, cand(i_cand, :)) < 0;
        end
        inside = inside & ~all(abs(cand) <= gridCfg.target_ub', 2);
        x0_dev = [x0_dev; cand(inside, :)]; %#ok<AGROW>
        if size(x0_dev, 1) >= n_want, break; end
    end
    x0_dev = x0_dev(1:min(n_want, size(x0_dev, 1)), :);
end

function axisCfg = axis_config(yml)
    phys = yml.dynamics;

    lon = yml.longitudinal;
    axisCfg.lon = grid_from_yml(lon, {'u', 'w', 'q', 'theta'});
    axisCfg.lon.u0_idx     = [5:12, 13, 3, 1];   % lift1..8, pusher, elevator, flap
    axisCfg.lon.plant_rows = [4; 6; 11; 8];      % plant rows for u, w, q, theta
    axisCfg.lon.delta_lb = [repmat(lon.input_min_Pi, 8, 1); lon.input_min_Pi_p; lon.input_min_delta_e; lon.input_min_delta_f];
    axisCfg.lon.delta_ub = [repmat(lon.input_max_Pi, 8, 1); lon.input_max_Pi_p; lon.input_max_delta_e; lon.input_max_delta_f];
    axisCfg.lon.phys_lb  = [repmat(phys.input_min_Pi, 8, 1); phys.input_min_Pi_p; phys.input_min_delta_e; phys.input_min_delta_f];
    axisCfg.lon.phys_ub  = [repmat(phys.input_max_Pi, 8, 1); phys.input_max_Pi_p; phys.input_max_delta_e; phys.input_max_delta_f];

    lat = yml.lateral;
    axisCfg.lat = grid_from_yml(lat, {'v', 'p', 'r', 'phi'});
    axisCfg.lat.u0_idx     = [5:12, 2, 4];       % lift1..8, aileron, rudder
    axisCfg.lat.plant_rows = [5; 10; 12; 7];     % plant rows for v, p, r, phi
    axisCfg.lat.delta_lb = [repmat(lat.input_min_Pi, 8, 1); lat.input_min_delta_a; lat.input_min_delta_r];
    axisCfg.lat.delta_ub = [repmat(lat.input_max_Pi, 8, 1); lat.input_max_delta_a; lat.input_max_delta_r];
    axisCfg.lat.phys_lb  = [repmat(phys.input_min_Pi, 8, 1); phys.input_min_delta_a; phys.input_min_delta_r];
    axisCfg.lat.phys_ub  = [repmat(phys.input_max_Pi, 8, 1); phys.input_max_delta_a; phys.input_max_delta_r];
end

function gridCfg = grid_from_yml(section, names)
    gridCfg.grid_min  = cellfun(@(s) section.(['grid_min_'    s]), names)';
    gridCfg.grid_max  = cellfun(@(s) section.(['grid_max_'    s]), names)';
    gridCfg.grid_num  = cellfun(@(s) section.(['grid_number_' s]), names)';
    gridCfg.target_ub = cellfun(@(s) section.(['target_max_'  s]), names)';
    gridCfg.grid_step = (gridCfg.grid_max - gridCfg.grid_min) ./ (gridCfg.grid_num - 1);
    gridCfg.grid_vec  = arrayfun(@(a, b, n) linspace(a, b, n), gridCfg.grid_min, ...
                                 gridCfg.grid_max, gridCfg.grid_num, 'UniformOutput', false)';
    gridCfg.held_dev  = zeros(4, 1);
end

function segs = tube_contour(V_stack, gridCfg, free_dims)
    % V = 0 boundary of the full-horizon tube (slice 1) on the free-dim plane,
    % held dims at their nearest node to held_dev. Used to frame the viz axes.
    idx = repmat({':'}, 1, 4);
    for d = setdiff(1:4, free_dims)
        [~, i_node] = min(abs(gridCfg.grid_vec{d} - gridCfg.held_dev(d)));
        idx{d} = i_node;
    end
    Z = squeeze(double(V_stack(1, idx{:})));
    segs = brt_zero_contour(gridCfg.grid_vec{free_dims(1)}, ...
                            gridCfg.grid_vec{free_dims(2)}, Z.');
end

function gridCfg = add_stack_indexing(gridCfg, stack_size)
    assert(isequal(stack_size(2:5), gridCfg.grid_num(:)'), ...
           'test_mc_verify:gridMismatch', ...
           'stack size %s does not match the yml grid %s.', ...
           mat2str(stack_size(2:5)), mat2str(gridCfg.grid_num(:)'));
    n = stack_size;
    gridCfg.n_slice = n(1);
    gridCfg.stride  = [n(1); n(1)*n(2); n(1)*n(2)*n(3); n(1)*n(2)*n(3)*n(4)];
    gridCfg.corners = double(dec2bin(0:15) - '0');
end

function V_stack = load_stack(fname)
    % Rebuild the 5-D [n_slice, Nu, Nw, Nq, Ntheta] array the rest of the code
    % indexes into, from the .mat time-stack's 1xK cell of 4-D slices.
    s       = load(fname, 'Vslices');
    V_stack = permute(cat(5, s.Vslices{:}), [5 1 2 3 4]);
end

function yml = read_yml(path)
    lines = regexp(fileread(path), '\r?\n', 'split');
    yml = struct();  section = '';
    for i_line = 1:numel(lines)
        line_ = regexprep(lines{i_line}, '#.*$', '');
        if isempty(strtrim(line_)), continue; end
        tok = regexp(line_, '^(\s*)([^:\s][^:]*):\s*(.*)$', 'tokens', 'once');
        if isempty(tok), continue; end
        indent = numel(tok{1});  key = strtrim(tok{2});  val = strtrim(tok{3});
        if indent == 0
            if isempty(val), section = key;  yml.(section) = struct(); end
            continue;
        end
        if isempty(section), continue; end
        num = str2double(val);
        if isnan(num), yml.(section).(key) = val; else, yml.(section).(key) = num; end
    end
end

function save_result(fname, result)
    save(fname, '-struct', 'result');
end
