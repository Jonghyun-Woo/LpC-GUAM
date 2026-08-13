classdef ValueFunction < handle
    % Loads <prefix>_UH{i}_WH{j}.mat BRT value functions (prefix per axis)
    % and provides scheduled, interpolated value V and gradient gradV at a
    % query perturbation state.
    %
    % Units: ft/s, ft/s, rad/s, rad. NO m/deg scaling (that scaling in helperOC
    % visualize_*_tube.m is visualization-only).
    %
    % Data convention (helperOC createGrid): data(i,j,k,l) <-> (x1_i,x2_j,x3_k,x4_l).
    % Value functions are time-converged (DtV = 0); the stored 'data' is 4D.
    %
    % Scheduling: BRT files exist per trim grid point (UH_idx, WH_idx). The
    % constructor parses which (UH_idx, WH_idx) are present, and query()
    % brackets (uh, wh) among the available breakpoints. Outside the available
    % coverage, ok=false so the caller passes through.
    %
    % Reference: docs/refs/liveness-filter.md, helperOC GUAM_HJIR.m / Config.m.

    properties
        axis_spec       % FilterConfig.axisSpec(axis) struct
        available       % true if at least one BRT table was loaded
        grid_vectors    % 1x4 cell of grid vectors [ft/s, ft/s, rad/s, rad]
        grid_min        % 4x1
        grid_max        % 4x1

        uh_breakpoint   % full UH breakpoint vector [ft/s] (trim table S.UH)
        wh_breakpoint   % full WH breakpoint vector [ft/s] (trim table S.WH)

        avail_uh_idx    % sorted available UH indices (into uh_breakpoint)
        avail_wh_idx    % sorted available WH indices (into wh_breakpoint)
        avail_uh_vel    % uh_breakpoint(avail_uh_idx)  (sorted velocities)
        avail_wh_vel    % wh_breakpoint(avail_wh_idx)
        uh_band         % coverage half-band on UH axis [ft/s]
        wh_band         % coverage half-band on WH axis [ft/s]

        brt_tables      % struct array: .uh_idx .wh_idx .value_interp .grad_interp(1x4)
        table_lookup    % numel(avail_uh_idx) x numel(avail_wh_idx) -> index into tables (NaN if missing)
    end

    methods
        function obj = ValueFunction(axis_spec, tables_dir, uh_breakpoint, wh_breakpoint)
            % axis_spec     : FilterConfig.axisSpec(ax)
            % tables_dir    : directory scanned for axis_spec.brt_prefix + '_UH%d_WH%d.mat'
            % uh_breakpoint,
            % wh_breakpoint : UH/WH breakpoint vectors [ft/s] (trim table S.UH, S.WH)
            obj.axis_spec       = axis_spec;
            obj.grid_min        = axis_spec.grid_min(:);
            obj.grid_max        = axis_spec.grid_max(:);
            obj.uh_breakpoint   = uh_breakpoint(:);
            obj.wh_breakpoint   = wh_breakpoint(:);

            % Grid vectors and per-axis spacing (uniform).
            obj.grid_vectors    = cell(1, 4);
            grid_step           = zeros(4, 1);
            for d = 1:4
                obj.grid_vectors{d} = linspace(axis_spec.grid_min(d), axis_spec.grid_max(d), axis_spec.grid_num(d));
                grid_step(d)        = (axis_spec.grid_max(d) - axis_spec.grid_min(d)) / (axis_spec.grid_num(d) - 1);
            end

            % Scan directory for value-function files of this axis.
            % Resolve tables_dir against cwd first, then the MATLAB path (so a
            % relative 'tables/BRT' works when Refactoring/ is on the path but
            % is not the current folder, matching how RSLQR loads tables).
            resolved_dir = ValueFunction.resolve_dir(tables_dir);
            file_pattern = sprintf('%s_UH*_WH*.mat', axis_spec.brt_prefix);
            files        = dir(fullfile(resolved_dir, file_pattern));
            if isempty(files)
                obj.available = false;
                return;
            end

            name_regex = sprintf('^%s_UH(\\d+)_WH(\\d+)\\.mat$', axis_spec.brt_prefix);
            tables = struct('uh_idx', {}, 'wh_idx', {}, 'value_interp', {}, 'grad_interp', {});
            for f = 1:numel(files)
                tokens = regexp( (f).name, name_regex, 'tokens', 'once');
                if isempty(tokens)
                    continue;
                end
                uh_idx = str2double(tokens{1});
                wh_idx = str2double(tokens{2});

                brt_value_table = load(fullfile(resolved_dir, files(f).name));
                assert(isfield(brt_value_table, 'data'), 'ValueFunction:noData', ...
                       '%s has no ''data'' field.', files(f).name);

                data = brt_value_table.data;
                assert(isequal(size(data), axis_spec.grid_num(:)'), ...
                       'ValueFunction:badSize', ...
                       '%s size %s ~= grid_num %s.', files(f).name, ...
                       mat2str(size(data)), mat2str(axis_spec.grid_num(:)'));

                grad_interp = cell(1, 4);
                for d = 1:4
                    grad_data      = ValueFunction.grad_along(data, d, grid_step(d));
                    grad_interp{d} = griddedInterpolant(obj.grid_vectors, grad_data, 'linear', 'nearest');
                end
                k = numel(tables) + 1;
                tables(k).uh_idx       = uh_idx;
                tables(k).wh_idx       = wh_idx;
                tables(k).value_interp = griddedInterpolant(obj.grid_vectors, data, 'linear', 'nearest');
                tables(k).grad_interp  = grad_interp;
            end

            if isempty(tables)
                obj.available = false;
                return;
            end
            obj.brt_tables  = tables;
            obj.available   = true;

            % Available breakpoints and coverage bands.
            obj.avail_uh_idx = unique([tables.uh_idx]);
            obj.avail_wh_idx = unique([tables.wh_idx]);
            obj.avail_uh_vel = obj.uh_breakpoint(obj.avail_uh_idx)';
            obj.avail_wh_vel = obj.wh_breakpoint(obj.avail_wh_idx)';
            obj.uh_band      = ValueFunction.axis_band(obj.avail_uh_vel, obj.uh_breakpoint);
            obj.wh_band      = ValueFunction.axis_band(obj.avail_wh_vel, obj.wh_breakpoint);

            % (avail_uh x avail_wh) -> table index lookup (NaN where absent).
            n_uh = numel(obj.avail_uh_idx);
            n_wh = numel(obj.avail_wh_idx);
            obj.table_lookup = nan(n_uh, n_wh);
            for k = 1:numel(tables)
                uh_pos = find(obj.avail_uh_idx == tables(k).uh_idx, 1);
                wh_pos = find(obj.avail_wh_idx == tables(k).wh_idx, 1);
                obj.table_lookup(uh_pos, wh_pos) = k;
            end
        end

        function [V, gradV, ok] = query(obj, x, uh, wh)
            % x : 4x1 perturbation state [u; w; q; theta] (ft/s, ft/s, rad/s, rad)
            % uh, wh : scheduling velocities [ft/s]
            % Returns interpolated V (scalar), gradV (4x1), ok (coverage flag).
            V = 0; gradV = zeros(4, 1); ok = false;
            if ~obj.available
                return;
            end

            % Bracket scheduling velocities among available breakpoints.
            [uh_lo, uh_hi, uh_frac, uh_ok] = ValueFunction.bracket(obj.avail_uh_vel, uh, obj.uh_band);
            [wh_lo, wh_hi, wh_frac, wh_ok] = ValueFunction.bracket(obj.avail_wh_vel, wh, obj.wh_band);
            if ~(uh_ok && wh_ok)
                return;
            end

            % The four corner tables must all be present.
            table_idx_LL = obj.table_lookup(uh_lo, wh_lo);
            table_idx_HL = obj.table_lookup(uh_hi, wh_lo);
            table_idx_LH = obj.table_lookup(uh_lo, wh_hi);
            table_idx_HH = obj.table_lookup(uh_hi, wh_hi);
            if any(isnan([table_idx_LL, table_idx_HL, table_idx_LH, table_idx_HH]))
                return;
            end

            % Clamp state to grid bounds (nearest-face extrapolation is handled
            % by the interpolants, but clamp keeps queries well-defined).
            x_clamped = min(max(x(:), obj.grid_min), obj.grid_max);

            weight_LL = (1 - uh_frac) * (1 - wh_frac);
            weight_HL =      uh_frac  * (1 - wh_frac);
            weight_LH = (1 - uh_frac) *      wh_frac;
            weight_HH =      uh_frac  *      wh_frac;

            [value_ll, gradV_ll] = obj.eval_table(table_idx_LL, x_clamped);
            [value_hl, gradV_hl] = obj.eval_table(table_idx_HL, x_clamped);
            [value_lh, gradV_lh] = obj.eval_table(table_idx_LH, x_clamped);
            [value_hh, gradV_hh] = obj.eval_table(table_idx_HH, x_clamped);

            V     = weight_LL * value_ll...
                  + weight_HL * value_hl...
                  + weight_LH * value_lh...
                  + weight_HH * value_hh;
            gradV = weight_LL * gradV_ll...
                  + weight_HL * gradV_hl...
                  + weight_LH * gradV_lh...
                  + weight_HH * gradV_hh;
            ok = true;
        end
    end

    methods (Access = private)
        function [V, G] = eval_table(obj, table_idx, x_clamped)
            % Evaluate value and gradient interpolants of table table_idx at clamped x_clamped.
            brt_table   = obj.brt_tables(table_idx);
            V           = brt_table.value_interp(x_clamped(1), x_clamped(2), x_clamped(3), x_clamped(4));
            G           = zeros(4, 1);
            for d = 1:4
                G(d)    = brt_table.grad_interp{d}(x_clamped(1), x_clamped(2), x_clamped(3), x_clamped(4));
            end
        end
    end

    methods (Static, Access = private)
        function d = resolve_dir(tables_dir)
            % Resolve a possibly-relative directory against cwd, then the
            % MATLAB path. Returns tables_dir unchanged if neither resolves
            % (constructor then finds no files and sets available=false).
            if isfolder(tables_dir)
                d = tables_dir;
                return;
            end
            w = what(tables_dir);
            if ~isempty(w)
                d = w(1).path;
                return;
            end
            d = tables_dir;
        end

        function g = grad_along(data, d, h)
            % Central-difference gradient of data along dimension d, uniform
            % spacing h, one-sided at the boundary slices.
            n  = size(data, d);
            g  = (circshift(data, -1, d) - circshift(data, 1, d)) / (2 * h);
            nd = ndims(data);
            idx = repmat({':'}, 1, nd);
            i1 = idx; i1{d} = 1;
            i2 = idx; i2{d} = 2;
            ie = idx; ie{d} = n;
            ip = idx; ip{d} = n - 1;
            g(i1{:}) = (data(i2{:}) - data(i1{:})) / h;
            g(ie{:}) = (data(ie{:}) - data(ip{:})) / h;
        end

        function band = axis_band(avail_vel, full_bp)
            % Coverage half-band for one scheduling axis. For >=2 available
            % breakpoints use half the min spacing between them; for a single
            % available breakpoint use half the min spacing of the full grid.
            if numel(avail_vel) >= 2
                band = 0.5 * min(diff(sort(avail_vel)));
            else
                fb = sort(full_bp(:));
                band = 0.5 * min(diff(fb));
            end
        end

        function [lo, hi, t, ok] = bracket(vals, q, band)
            % Bracket scalar q among sorted available velocities vals. Returns
            % lower/upper indices lo,hi into vals, fraction t in [0,1], and a
            % coverage flag ok. Beyond the range, coverage extends by band
            % (nearest-neighbor), else ok=false.
            n = numel(vals);
            if n == 1
                ok = abs(q - vals(1)) <= band;
                lo = 1; hi = 1; t = 0; return;
            end
            if q <= vals(1)
                ok = (q >= vals(1) - band); lo = 1; hi = 1; t = 0; return;
            end
            if q >= vals(end)
                ok = (q <= vals(end) + band); lo = n; hi = n; t = 0; return;
            end
            lo = find(vals <= q, 1, 'last');
            hi = lo + 1;
            t  = (q - vals(lo)) / (vals(hi) - vals(lo));
            ok = true;
        end
    end
end
