classdef ValueFunctionLUT < handle
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
        spec            % FilterConfig.axisSpec(axis) struct
        available       % true if at least one BRT table was loaded
        grid_vectors    % 1x4 cell of grid vectors [ft/s, ft/s, rad/s, rad]
        grid_min        % 4x1
        grid_max        % 4x1

        uh_breakpoint   % full UH breakpoint vector [ft/s] (trim table S.UH)
        wh_breakpoint   % full WH breakpoint vector [ft/s] (trim table S.WH)

        availUH_idx     % sorted available UH indices (into uh_breakpoint)
        availWH_idx     % sorted available WH indices (into wh_breakpoint)
        availUH_vel     % uh_breakpoint(availUH_idx)  (sorted velocities)
        availWH_vel     % wh_breakpoint(availWH_idx)
        bandUH          % coverage half-band on UH axis [ft/s]
        bandWH          % coverage half-band on WH axis [ft/s]

        brt_tables      % struct array: .uh_idx .wh_idx .Vinterp .GradVinterp(1x4)
        table_lookup    % numel(availUH_idx) x numel(availWH_idx) -> index into tables (NaN if missing)
    end

    methods
        function obj = ValueFunctionLUT(spec, tables_dir, uh_breakpoint, wh_breakpoint)
            % spec          : FilterConfig.axisSpec(ax)
            % tables_dir    : directory scanned for spec.brt_prefix + '_UH%d_WH%d.mat'
            % uh_breakpoint,
            % wh_breakpoint : UH/WH breakpoint vectors [ft/s] (trim table S.UH, S.WH)
            obj.spec            = spec;
            obj.grid_min        = spec.grid_min(:);
            obj.grid_max        = spec.grid_max(:);
            obj.uh_breakpoint   = uh_breakpoint(:);
            obj.wh_breakpoint   = wh_breakpoint(:);

            % Grid vectors and per-axis spacing (uniform).
            obj.grid_vectors    = cell(1, 4);
            grid_step           = zeros(4, 1);
            for d = 1:4
                obj.grid_vectors{d} = linspace(spec.grid_min(d), spec.grid_max(d), spec.grid_num(d));
                grid_step(d)        = (spec.grid_max(d) - spec.grid_min(d)) / (spec.grid_num(d) - 1);
            end

            % Scan directory for value-function files of this axis.
            % Resolve tables_dir against cwd first, then the MATLAB path (so a
            % relative 'tables/BRT' works when Refactoring/ is on the path but
            % is not the current folder, matching how RSLQR loads tables).
            resolved = ValueFunctionLUT.resolve_dir(tables_dir);
            pat   = sprintf('%s_UH*_WH*.mat', spec.brt_prefix);
            files = dir(fullfile(resolved, pat));
            if isempty(files)
                obj.available = false;
                return;
            end

            rex = sprintf('^%s_UH(\\d+)_WH(\\d+)\\.mat$', spec.brt_prefix);
            tbl = struct('uh_idx', {}, 'wh_idx', {}, 'Vinterp', {}, 'GradVinterp', {});
            for f = 1:numel(files)
                tok = regexp(files(f).name, rex, 'tokens', 'once');
                if isempty(tok)
                    continue;
                end
                ui = str2double(tok{1});
                wi = str2double(tok{2});

                brt_value_table = load(fullfile(resolved, files(f).name));
                assert(isfield(brt_value_table, 'data'), 'ValueFunctionLUT:noData', ...
                       '%s has no ''data'' field.', files(f).name);

                data = brt_value_table.data;
                assert(isequal(size(data), spec.grid_num(:)'), ...
                       'ValueFunctionLUT:badSize', ...
                       '%s size %s ~= grid_num %s.', files(f).name, ...
                       mat2str(size(data)), mat2str(spec.grid_num(:)'));

                GradVinterp = cell(1, 4);
                for d = 1:4
                    G = ValueFunctionLUT.grad_along(data, d, grid_step(d));
                    GradVinterp{d} = griddedInterpolant(obj.grid_vectors, G, 'linear', 'nearest');
                end
                k = numel(tbl) + 1;
                tbl(k).uh_idx  = ui;
                tbl(k).wh_idx  = wi;
                tbl(k).Vinterp = griddedInterpolant(obj.grid_vectors, data, 'linear', 'nearest');
                tbl(k).GradVinterp = GradVinterp;
            end

            if isempty(tbl)
                obj.available = false;
                return;
            end
            obj.brt_tables  = tbl;
            obj.available   = true;

            % Available breakpoints and coverage bands.
            obj.availUH_idx = unique([tbl.uh_idx]);
            obj.availWH_idx = unique([tbl.wh_idx]);
            obj.availUH_vel = obj.uh_breakpoint(obj.availUH_idx)';
            obj.availWH_vel = obj.wh_breakpoint(obj.availWH_idx)';
            obj.bandUH      = ValueFunctionLUT.axis_band(obj.availUH_vel, obj.uh_breakpoint);
            obj.bandWH      = ValueFunctionLUT.axis_band(obj.availWH_vel, obj.wh_breakpoint);

            % (availUH x availWH) -> table index lookup (NaN where absent).
            nU = numel(obj.availUH_idx);
            nW = numel(obj.availWH_idx);
            obj.table_lookup = nan(nU, nW);
            for k = 1:numel(tbl)
                iu = find(obj.availUH_idx == tbl(k).uh_idx, 1);
                iw = find(obj.availWH_idx == tbl(k).wh_idx, 1);
                obj.table_lookup(iu, iw) = k;
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
            [uLo, uHi, tU, okU] = ValueFunctionLUT.bracket(obj.availUH_vel, uh, obj.bandUH);
            [wLo, wHi, tW, okW] = ValueFunctionLUT.bracket(obj.availWH_vel, wh, obj.bandWH);
            if ~(okU && okW)
                return;
            end

            % The four corner tables must all be present.
            table_idx_LL = obj.table_lookup(uLo, wLo);
            table_idx_HL = obj.table_lookup(uHi, wLo);
            table_idx_LH = obj.table_lookup(uLo, wHi);
            table_idx_HH = obj.table_lookup(uHi, wHi);
            if any(isnan([table_idx_LL, table_idx_HL, table_idx_LH, table_idx_HH]))
                return;
            end

            % Clamp state to grid bounds (nearest-face extrapolation is handled
            % by the interpolants, but clamp keeps queries well-defined).
            xc = min(max(x(:), obj.grid_min), obj.grid_max);
            
            interp_weight_LLp = (1 - tU) * (1 - tW);
            interp_weight_HLp =      tU  * (1 - tW);
            interp_weight_LHp = (1 - tU) *      tW;
            interp_weight_HHp =      tU  *      tW;

            [value_ll, gradV_ll] = obj.eval_table(table_idx_LL, xc);
            [value_hl, gradV_hl] = obj.eval_table(table_idx_HL, xc);
            [value_lh, gradV_lh] = obj.eval_table(table_idx_LH, xc);
            [value_hh, gradV_hh] = obj.eval_table(table_idx_HH, xc);

            V     = interp_weight_LLp * value_ll...
                  + interp_weight_HLp * value_hl...
                  + interp_weight_LHp * value_lh...
                  + interp_weight_HHp * value_hh;
            gradV = interp_weight_LLp * gradV_ll...
                  + interp_weight_HLp * gradV_hl...
                  + interp_weight_LHp * gradV_lh...
                  + interp_weight_HHp * gradV_hh;
            ok = true;
        end
    end

    methods (Access = private)
        function [V, G] = eval_table(obj, table_idx, xc)
            % Evaluate value and gradient interpolants of table table_idx at clamped xc.
            brt_table   = obj.brt_tables(table_idx);
            V           = brt_table.Vinterp(xc(1), xc(2), xc(3), xc(4));
            G           = zeros(4, 1);
            for d = 1:4
                G(d)    = brt_table.GradVinterp{d}(xc(1), xc(2), xc(3), xc(4));
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
