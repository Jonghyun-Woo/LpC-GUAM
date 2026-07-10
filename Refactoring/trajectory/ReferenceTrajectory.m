classdef ReferenceTrajectory
    % REFERENCETRAJECTORY m-code-only reference trajectory generator.
    %
    % Gathers every closed-loop reference scenario in one place. Unlike the
    % Exec_Scripts/*.m generators (which build 3-point timeseries for Simulink
    % interpolation), this returns the dense per-step tables that the
    % Refactoring/ port consumes directly. No Simulink and no lib/ dependency.
    %
    % Scenarios (both fly straight north, heading chi = 0):
    %   'althold' (default) - hover climb to 80 ft, then hold 80 ft in cruise
    %   'climb'             - hover climb to 80 ft, then climb to 100 ft (original)
    %
    % Frame/sign conventions: NED position, Down negative = higher altitude.
    % vel is expressed in the heading frame; since chi = 0 here it equals the
    % inertial velocity. A chi ~= 0 scenario would need a QrotZ(chi) transform
    % of the inertial velocity into the heading frame before returning.

    methods (Static)
        function ref = build(scenario, dt, T)
            % BUILD reference trajectory table for a scenario.
            %   scenario : 'althold' (default) | 'climb'
            %   dt, T    : time grid parameters (from SimConfig)
            %   ref      : struct with fields
            %                time   (1xN)  - time grid [s]
            %                pos    (3xN)  - NED position [ft]
            %                vel    (3xN)  - heading-frame velocity [ft/s]
            %                chi    (1xN)  - heading angle [rad]
            %                chidot (1xN)  - heading rate [rad/s]
            if nargin < 1 || isempty(scenario), scenario = 'althold'; end

            time = 0 : dt : T;
            N    = numel(time);

            % Verification-only forward-transition scenario for the LON
            % liveness filter. Mission = LEVEL forward transition: forward
            % speed ramps 17 -> 59 ft/s (crossing several UH trim tables so
            % the filter's UH scheduling is exercised) at CONSTANT altitude.
            % The test harness starts the vehicle off-trim; the WH3
            % (descending) BRT is the safety envelope and the filter, anchored
            % to WH3 (RSLQR safety_wh_anchor), engages only if the transition
            % would leave that envelope. Used by verify_liveness_lon_realbrt.m.
            %
            % Frame rule: vel(3) is the heading/inertial vertical velocity, so
            % level flight is vel(3) = 0. Do NOT put the body-frame WH here (an
            % earlier bug integrated WH3 body-w as an inertial descent rate).
            % The WH3 anchor is pinned by the filter, not encoded here.
            if strcmp(scenario, 'lon_brt_verify')
                pos = zeros(3, N);
                vel = zeros(3, N);
                vel(1, :) = linspace(17, 59, N);    % forward speed sweep across UH tables
                vel(3, :) = 0;                      % level flight (altitude hold)
                pos(1, :) = cumtrapz(time, vel(1, :));
                pos(3, :) = -300 * ones(1, N);      % hold 300 ft up
                ref = struct('time',   time, ...
                             'pos',    pos, ...
                             'vel',    vel, ...
                             'chi',    zeros(1, N), ...
                             'chidot', zeros(1, N));
                return;
            end

            % Split hover (first half, 0..T/2) from cruise (second half).
            % Time-based split so the segment boundary sits at t = T/2
            % regardless of grid length. For T = 40, dt = 0.01 this gives
            % nHover = 2001, nCruise = 2000 (matches the former hardcoded table).
            nHover  = nnz(time <= T / 2);
            nCruise = N - nHover;

            pos = zeros(3, N);
            vel = zeros(3, N);

            % North: hold during hover, accelerate forward to 150 ft in cruise
            pos(1, :) = [zeros(1, nHover), linspace(0, 150, nCruise)];
            vel(1, :) = [zeros(1, nHover), linspace(0, 15, nCruise)];

            % East: unused (straight north)
            % (pos(2,:) and vel(2,:) stay zero)

            % Down velocity: climb (w: -8 -> 0) during hover, zero in cruise
            vel(3, :) = [linspace(-8, 0, nHover), zeros(1, nCruise)];

            % Down position: differs only in the cruise segment
            switch scenario
                case 'climb'
                    % keep climbing 80 -> 100 ft (original, inconsistent with w = 0)
                    pos(3, :) = [linspace(0, -80, nHover), linspace(-80, -100, nCruise)];
                case 'althold'
                    % hold 80 ft through cruise (consistent with w = 0)
                    pos(3, :) = [linspace(0, -80, nHover), -80 * ones(1, nCruise)];
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
