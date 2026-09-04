classdef BRTValue < handle
    % BRTValue Class
    % Safety value V(x, v) for the longitudinal state x under speed command v.
    % V < 0 means inside the safe set, V >= 0 means outside.
    %
    % A reachable set was computed at each of the 20 anchor speeds only, as a
    % function of the DEVIATION from that anchor's trim state. Between anchors
    % there is no set, so a rule is needed. Two are provided (see `mode`):
    %
    %   'union' (default) : the aircraft is safe if it is inside SOME anchor's
    %       set. Each set is evaluated from its own anchor, which is the only
    %       offset it was ever computed for, so this asserts nothing that was
    %       not certified. The price is that mid-segment the aircraft is half a
    %       step (4.2 ft/s) away from both anchors, and that offset is charged
    %       against both sets.
    %
    %   'blend' : interpolate the two anchors into one moving reference point,
    %       take the deviation from that, and blend the two set values. Cheaper
    %       on the aircraft mid-segment, but the set it tests against was never
    %       computed - it is a blend of two sets, not a reachable set of
    %       anything. This was the original behaviour; every result recorded
    %       before 2026-08-14 (33.7 s band schedule, check_a_law, the governor
    %       work) was measured with it, so set mode='blend' to reproduce those.

    properties
        mode = 'union'  % 'union' | 'blend'  (see class comment)

        gv          % 1x4 cell, BRT grid vectors
        brt         % 1xN cell, 4-D value arrays
        trim_lon    % 4xN anchor trim states [u; w; q; theta]
        u_anchor    % 1xN anchor forward speeds [ft/s] (monotone increasing)
        n_anchor

        g_lo        % 1x4 first grid point per axis   ) cached so the lookup
        g_hi        % 1x4 last  grid point per axis   ) can index the grid by
        g_h         % 1x4 grid spacing per axis       ) formula instead of
        g_n         % 1x4 number of grid points       ) searching it
    end

    methods
        function obj = BRTValue(data_dir, trim_lon, wh_idx)
            % data_dir : folder holding GUAM_LON_BRT_HJIR_UH*_WH*.mat
            % trim_lon : 4 x N anchor trim states
            % wh_idx   : WH breakpoint index of the BRT family
            spec = FilterConfig.channelSpec('lon');
            obj.gv = cell(1, 4);
            for d = 1:4
                obj.gv{d} = linspace(spec.grid_min(d), spec.grid_max(d), ...
                                     spec.grid_num(d));
            end
            obj.g_lo = cellfun(@(g) g(1),   obj.gv);
            obj.g_hi = cellfun(@(g) g(end), obj.gv);
            obj.g_n  = cellfun(@numel,      obj.gv);
            obj.g_h  = (obj.g_hi - obj.g_lo) ./ (obj.g_n - 1);

            obj.trim_lon = trim_lon;
            obj.u_anchor = trim_lon(1, :);
            obj.n_anchor = size(trim_lon, 2);
            obj.brt = cell(1, obj.n_anchor);
            for k = 1:obj.n_anchor
                S = load(fullfile(data_dir, ...
                    sprintf('GUAM_LON_BRT_HJIR_UH%d_WH%d.mat', k, wh_idx)), 'data');
                obj.brt{k} = S.data;
            end
        end

        function xe = trim_at(obj, v)
            % x_e(v): 4x1 trim state whose forward speed is v.
            [k, a] = obj.bracket(v);
            xe = (1 - a) * obj.trim_lon(:, k) + a * obj.trim_lon(:, k + 1);
        end

        function V = value(obj, x_lon, v)
            % V(x, v). x_lon is the ABSOLUTE longitudinal state [u;w;q;theta].
            [k, a] = obj.bracket(v);

            if strcmp(obj.mode, 'union')
                % Inside anchor k's set, or inside anchor k+1's set. Each is
                % read at the deviation from ITS OWN anchor.
                V = min(obj.value_at(x_lon, k), obj.value_at(x_lon, k + 1));
                return;
            end

            % 'blend': one moving reference point, one deviation, blended value
            xe = (1 - a) * obj.trim_lon(:, k) + a * obj.trim_lon(:, k + 1);
            [i0, f] = obj.locate(x_lon(:)' - xe(:)');
            V = (1 - a) * lerp4(obj.brt{k},     i0, f) + ...
                     a  * lerp4(obj.brt{k + 1}, i0, f);
        end

        function V = value_at(obj, x_lon, j)
            % Anchor j's own set, read at the deviation from anchor j.
            [i0, f] = obj.locate(x_lon(:)' - obj.trim_lon(:, j)');
            V = lerp4(obj.brt{j}, i0, f);
        end
    end

    methods (Access = private)
        function [i0, f] = locate(obj, xc)
            % Bracket index and fraction per axis for a deviation xc (1x4).
            % The grids are linspace, so this is a formula; interpn would
            % search them and re-validate the whole call every time. Same
            % arithmetic, agrees with interpn to 4e-15, but ~23x faster - and
            % this runs twice per rollout step, 600 steps per rollout, ~33
            % rollouts per governor solve.
            i0 = [0 0 0 0];  f = [0 0 0 0];
            for d = 1:4
                t = (xc(d) - obj.g_lo(d)) / obj.g_h(d) + 1;
                if t < 1
                    t = 1;
                elseif t > obj.g_n(d)
                    t = obj.g_n(d);
                end
                b = floor(t);
                if b > obj.g_n(d) - 1, b = obj.g_n(d) - 1; end
                i0(d) = b;  f(d) = t - b;
            end
        end

        function [k, a] = bracket(obj, v)
            % Bracket v among the anchor speeds; clamp outside the range.
            u = obj.u_anchor;
            if v <= u(1)
                k = 1;  a = 0;  return;
            end
            if v >= u(end)
                k = obj.n_anchor - 1;  a = 1;  return;
            end
            k = find(u <= v, 1, 'last');
            k = min(k, obj.n_anchor - 1);
            a = (v - u(k)) / (u(k + 1) - u(k));
        end
    end
end

function V = lerp4(D, i0, f)
% 4-D linear interpolation on the 2x2x2x2 block at i0, weights f.
% Collapsed one axis at a time; identical to interpn(...,'linear').
c = D(i0(1):i0(1)+1, i0(2):i0(2)+1, i0(3):i0(3)+1, i0(4):i0(4)+1);
c = c(:,:,:,1) * (1 - f(4)) + c(:,:,:,2) * f(4);
c = c(:,:,1)   * (1 - f(3)) + c(:,:,2)   * f(3);
c = c(:,1)     * (1 - f(2)) + c(:,2)     * f(2);
V = c(1)       * (1 - f(1)) + c(2)       * f(1);
end
