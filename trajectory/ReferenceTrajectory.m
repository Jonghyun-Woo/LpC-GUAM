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
        function topts = trimOpts(dt, overrides)
            % Options handed to TrimScheduleTrajectory.build. Kept public so
            % a governor script can build the SAME schedule the plant sees:
            %   traj = TrimScheduleTrajectory.build( ...
            %              ReferenceTrajectory.trimOpts(dt, params));
            if nargin < 2 || isempty(overrides), overrides = struct(); end
            topts = struct('dt', dt);
            T_seg = getfield_default(overrides, 'T_seg', []);
            if ~isempty(T_seg), topts.T_seg = T_seg; end
        end

        function ref = build(scenario, dt, T, target_vel, overrides)
            % BUILD reference trajectory table for a scenario.
            %   scenario  : 'althold' (default) | 'climb' | 'trim_schedule'
            %   dt, T     : time grid parameters (from SimConfig)
            %   overrides : optional struct; .T_seg sets the per-segment
            %               time of the 'trim_schedule' mission so the
            %               plant-side reference and the governor's own
            %               schedule cannot drift apart
            %   ref      : struct with fields
            %                time   (1xN)  - time grid [s]
            %                pos    (3xN)  - NED position [ft]
            %                vel    (3xN)  - heading-frame velocity [ft/s]
            %                chi    (1xN)  - heading angle [rad]
            %                chidot (1xN)  - heading rate [rad/s]
            if nargin < 1 || isempty(scenario), scenario = 'althold'; end
            if nargin < 5 || isempty(overrides), overrides = struct(); end

            % Trim-schedule mission for the guidance governor: forward speed
            % follows the 20-trim BRT anchor schedule (TrimScheduleTrajectory),
            % passing through every BRT target-set origin in (u, theta) trim
            % space. Sim-side we only need the position/velocity tables; the
            % trim-space bookkeeping (T_seg, anchors, progress) rides along in
            % ref.trim_schedule for the governor/liveness prediction.
            %
            % Frame rule (same as lon_brt_verify): vel(3) is the
            % heading/inertial vertical velocity -> level flight = 0.
            %
            % VERTICAL PROFILE. Two choices, and they are not equivalent.
            %
            %   descend_rate = 0 (default, shipped behaviour)
            %       Level at 80 ft. The 20 trim points are WH3 points, trimmed
            %       at an 11.667 ft/s descent, so flying level leaves the
            %       aircraft a standing -11.7 ft/s off every anchor on the
            %       vertical-speed axis, plus up to +2.2 deg in pitch. The
            %       reachable sets still contain it - V stays near -0.20 - but
            %       most of the room to manoeuvre is spent before the vehicle
            %       accelerates at all. Measured: level flight admits a uniform
            %       1.0 ft/s^2 where sitting on the anchors admits 3.5.
            %
            %   descend_rate = WH (11.667 ft/s, i.e. 700 ft/min)
            %       The condition the anchors and their sets were built for.
            %       Measured against the simulator, the aircraft then settles
            %       within 0.7 ft/s and 0.1 deg of every trim point, and the
            %       standing offset disappears. Costs altitude: descend_rate
            %       times the mission length, so the start altitude is set from
            %       end_alt plus that.
            %
            % Left off by default because every result recorded before
            % 2026-08-20 - the governor studies, the band schedule, the fixed
            % anchor timings - was flown level, and switching silently would
            % invalidate them.
            if strcmp(scenario, 'trim_schedule')
                tst = TrimScheduleTrajectory.build( ...
                          ReferenceTrajectory.trimOpts(dt, overrides));
                N   = numel(tst.time);
                vel = zeros(3, N);
                vel(1, :) = tst.lon(1, :);          % forward speed = trim UH schedule
                pos = zeros(3, N);
                pos(1, :) = cumtrapz(tst.time, vel(1, :));

                wdot = getfield_default(overrides, 'descend_rate', 0);
                if wdot == 0
                    pos(3, :) = -80;                % hold 80 ft up
                else
                    end_alt   = getfield_default(overrides, 'end_alt', 100);
                    vel(3, :) = wdot;               % + is down in NED
                    pos(3, :) = -(end_alt + wdot * (tst.time(end) - tst.time)) ;
                end
                ref = struct('time',   tst.time, ...
                             'pos',    pos, ...
                             'vel',    vel, ...
                             'chi',    zeros(1, N), ...
                             'chidot', zeros(1, N), ...
                             'trim_schedule', tst);
                return;
            end

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
                vel(1, :) = linspace(0, target_vel, N);    % forward speed sweep across UH tables
                vel(3, :) = 0;                      % level flight (altitude hold)
                pos(1, :) = cumtrapz(time, vel(1, :));
                pos(3, :) = -80 * ones(1, N);      % hold 80 ft up
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
