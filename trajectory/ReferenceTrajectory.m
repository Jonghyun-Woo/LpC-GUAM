classdef ReferenceTrajectory
    % REFERENCETRAJECTORY dense per-step reference tables (no Simulink).
    % NED position (Down negative = higher altitude); vel is heading-frame,
    % equal to inertial since chi = 0. All scenarios fly straight north.

    properties (Constant)
        alt = -200;
    end

    methods (Static)
        function ref = build(scenario, dt, T, target_vel)
            % BUILD reference table for a scenario.
            %   scenario  : 'althold' (default) | 'climb' | 'brt_verify'
            %   dt, T     : time grid (from SimConfig)
            %   ref       : struct with time, pos, vel, chi, chidot
            if nargin < 1 || isempty(scenario), scenario = 'althold'; end

            time = 0 : dt : T;
            N    = numel(time);

            % Level forward transition: u ramps 0 -> target_vel while body-w
            % tracks the level-flight trim schedule (vd = -sin(th)*u + cos(th)*w
            % = 0) to hold `alt`; above the trim grid's w range w saturates.
            if strcmp(scenario, 'brt_verify')
                S   = load('trim_table_Poly_ConcatVer4p0.mat');
                THg = squeeze(S.XU0_interp(11, :, :));   % trim pitch [nUH x nWH]
                pos = zeros(3, N);
                vel = zeros(3, N);
                vel(1, :) = linspace(0, target_vel, N);
                vel(3, :) = ReferenceTrajectory.level_body_w(vel(1, :), S.UH(:), S.WH(:), THg);
                pos(1, :) = cumtrapz(time, vel(1, :));
                pos(3, :) = ReferenceTrajectory.alt * ones(1, N);
                ref = struct('time',   time, ...
                             'pos',    pos, ...
                             'vel',    vel, ...
                             'chi',    zeros(1, N), ...
                             'chidot', zeros(1, N));
                return;
            end

            % Hover (0..T/2) then cruise, split by time so the boundary is
            % grid-independent.
            nHover  = nnz(time <= T / 2);
            nCruise = N - nHover;

            pos = zeros(3, N);
            vel = zeros(3, N);

            pos(1, :) = [zeros(1, nHover), linspace(0, 150, nCruise)];
            vel(1, :) = [zeros(1, nHover), linspace(0, 15, nCruise)];

            vel(3, :) = [linspace(-8, 0, nHover), zeros(1, nCruise)];

            % Down position differs only in cruise
            switch scenario
                case 'climb'
                    pos(3, :) = [linspace(0, ReferenceTrajectory.alt, nHover), linspace(ReferenceTrajectory.alt, ReferenceTrajectory.alt - 20, nCruise)];
                case 'althold'
                    pos(3, :) = [linspace(0, ReferenceTrajectory.alt, nHover), ReferenceTrajectory.alt * ones(1, nCruise)];
                otherwise
                    error('ReferenceTrajectory:unknownScenario', ...
                          'Unknown scenario ''%s'' (expected ''althold'' or ''climb'').', ...
                          scenario);
            end

            ref = struct('time',   time, ...
                         'pos',    pos, ...
                         'vel',    vel, ...
                         'chi',    zeros(1, N), ...
                         'chidot', zeros(1, N));
        end

        function w = level_body_w(u, UHb, WHb, THg)
            % Body-w [ft/s] giving level flight (vd = -sin(th)*u + cos(th)*w = 0)
            % at each forward speed u, from the trim pitch THg on the (UH,WH) grid;
            % saturates at the grid's body-w limits (level needs more w than the
            % grid holds at high u, leaving a small residual climb there).
            uq = linspace(min(u), max(u), 201);
            wg = linspace(WHb(1), WHb(end), 201);
            wn = zeros(size(uq));
            for i = 1:numel(uq)
                th = interp2(WHb', UHb, THg, wg, uq(i), 'linear');
                vd = -sin(th) .* uq(i) + cos(th) .* wg;
                if vd(end) <= 0
                    wn(i) = wg(end);
                elseif vd(1) >= 0
                    wn(i) = wg(1);
                else
                    wn(i) = interp1(vd, wg, 0, 'linear');
                end
            end
            w = interp1(uq, wn, u, 'linear');
        end
    end
end
