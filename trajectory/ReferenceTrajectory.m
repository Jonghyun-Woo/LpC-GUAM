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
            
            if strcmp(scenario, 'brt_verify')
                S   = load('trim_table_Poly_ConcatVer4p0.mat');
                THw = squeeze(S.XU0_interp(11, :, 2));       % trim pitch @wh=0 [nUH]
                pos = zeros(3, N);
                vel = zeros(3, N);
                vel(1, :) = linspace(0, target_vel, N);
                vel(3, :) = 0;                                  % wh = 0 (body w)
                u_clamp = min(max(vel(1, :), S.UH(1)), S.UH(end));
                th      = interp1(S.UH(:), THw(:), u_clamp, 'linear');
                vn      =  cos(th) .* vel(1, :);               % NED-north vel at w=0
                vd      = -sin(th) .* vel(1, :);               % NED-down  vel at w=0
                pos(1, :) = cumtrapz(time, vn);
                pos(3, :) = ReferenceTrajectory.alt + cumtrapz(time, vd);
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
    end
end
