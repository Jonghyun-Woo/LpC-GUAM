classdef TrimScheduleTrajectory
    % TRIMSCHEDULETRAJECTORY guidance-level reference trajectory in trim space.
    %
    % Connects the target-set origins (= trim points) of consecutive LON BRTs
    % with straight lines in the full 25-dim trim vector XU0, sampled on a
    % fixed time grid. The forward-speed schedule fixes the per-segment time
    % T_seg(k) = time allotted to go from the k-th BRT anchor to the (k+1)-th.
    %
    % Structure: [initial hold at trim 1] -> [n_trim-1 linear transition
    % segments] -> [final hold at trim n_trim].
    %
    % Units: ft/s, rad/s, rad (same as trim table / BRT grids).

    methods (Static)
        function traj = build(opts)
            % BUILD trim-schedule reference trajectory.
            %
            % opts fields (all optional):
            %   trim_table : path to trim table .mat
            %                (default '<repo>/tables/trim/trim_table_Poly_ConcatVer4p0.mat')
            %   n_trim     : number of trim points / BRT anchors (default 20)
            %   wh_idx     : WH breakpoint index of the BRT family (default 3)
            %   dt         : sample time [s]                       (default 0.02)
            %   t_hold     : initial/final hold duration [s]       (default 5)
            %   T_seg      : segment times [s]; scalar (uniform) or
            %                1x(n_trim-1) vector                   (default 2.0)
            %   accel      : forward-speed schedule slope [ft/s^2];
            %                sets T_seg(k) = (UH(k+1)-UH(k))/accel when
            %                T_seg is not given
            %
            % traj fields:
            %   time     : 1xN time grid [s]
            %   xu       : 25xN interpolated trim vector [X0(12); U0(13)]
            %   lon      : 4xN absolute LON state [u; w; q; theta]
            %   progress : 1xN trim-index progress in [1, n_trim]
            %              (e.g. 3.4 = 40 % through segment 3->4)
            %   seg      : 1xN active segment index in [1, n_trim-1]
            %              (holds clamp to the first/last segment)
            %   trims    : 25 x n_trim trim vectors (BRT target origins)
            %   trim_lon : 4 x n_trim  [u; w; q; theta] of each trim
            %   UH, WH   : scheduling velocities of the anchors [ft/s]
            %   T_seg    : 1x(n_trim-1) segment durations [s]
            %   t_node   : 1x n_trim time at which each trim point is crossed
            if nargin < 1, opts = struct(); end
            root = fileparts(fileparts(mfilename('fullpath')));
            trim_table = getfield_default(opts, 'trim_table', ...
                fullfile(root, 'tables', 'trim', 'trim_table_Poly_ConcatVer4p0.mat'));
            n_trim = getfield_default(opts, 'n_trim', 20);
            wh_idx = getfield_default(opts, 'wh_idx', 3);
            dt     = getfield_default(opts, 'dt', 0.02);
            t_hold = getfield_default(opts, 't_hold', 5.0);
            accel  = getfield_default(opts, 'accel', []);

            S = load(trim_table, 'XU0_interp', 'UH', 'WH');
            trims = S.XU0_interp(:, 1:n_trim, wh_idx);      % 25 x n_trim
            UH    = S.UH(1:n_trim)';                        % 1 x n_trim
            WH    = S.WH(wh_idx);

            % LON state rows of XU0: u=1, w=3, q=5, theta=11 (RSLQR convention)
            lon_rows = [1, 3, 5, 11];
            trim_lon = trims(lon_rows, :);                  % 4 x n_trim

            if isfield(opts, 'T_seg') && ~isempty(opts.T_seg)
                T_seg = opts.T_seg;
            elseif ~isempty(accel)
                T_seg = diff(UH) / accel;
            else
                T_seg = 2.0;                        % uniform default [s]
            end
            if isscalar(T_seg)
                T_seg = repmat(T_seg, 1, n_trim - 1);
            end
            T_seg = T_seg(:)';
            assert(numel(T_seg) == n_trim - 1, ...
                'TrimScheduleTrajectory:Tseg', ...
                'T_seg must have n_trim-1 = %d entries.', n_trim - 1);

            % --- assemble sample-count layout -----------------------------
            n_hold = round(t_hold / dt);
            n_k    = round(T_seg / dt);                     % samples per segment
            N      = 2 * n_hold + sum(n_k) + 1;             % +1 final sample

            time     = (0:N-1) * dt;
            progress = zeros(1, N);

            i = 1;
            progress(i : i + n_hold - 1) = 1;               % initial hold
            i = i + n_hold;
            for k = 1:n_trim - 1                            % transitions
                frac = (0:n_k(k) - 1) / n_k(k);
                progress(i : i + n_k(k) - 1) = k + frac;
                i = i + n_k(k);
            end
            progress(i:end) = n_trim;                       % final hold

            % --- linear interpolation on the full trim vector -------------
            seg  = min(floor(progress), n_trim - 1);
            frac = progress - seg;
            xu   = trims(:, seg) + (trims(:, seg + 1) - trims(:, seg)) .* frac;

            % time at which each trim node is crossed (for plotting)
            t_node = t_hold + [0, cumsum(n_k)] * dt;

            traj = struct('time', time, 'xu', xu, 'lon', xu(lon_rows, :), ...
                          'progress', progress, 'seg', seg, ...
                          'trims', trims, 'trim_lon', trim_lon, ...
                          'UH', UH, 'WH', WH, 'T_seg', T_seg, ...
                          't_node', t_node, 'dt', dt, 'n_trim', n_trim, ...
                          'wh_idx', wh_idx);
        end
    end
end
