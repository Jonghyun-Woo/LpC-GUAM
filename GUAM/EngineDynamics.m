classdef EngineDynamics < handle
    % ENGINEDYNAMICS First-order rate- and position-limited servo bank for
    % the 9 rotor/propeller speed actuators (8 lift rotors + 1 pusher).
    %
    % Mirrors GUAM's ServosLib "First Order Dynamic Rate and Position
    % Limited Servo" applied to engines in the SurfEng subsystem:
    %   n_dot = sat(wn*(cmd - n), +/-rate_limit), n clamped to limits
    % Parameters from vehicles/Lift+Cruise/setup/setupEngines.m
    % (Polynomial force/moment model limits).
    properties
        dt              % integration time step [s]
        wn              % servo bandwidth [rad/s]
        rate_limit      % rotor speed rate limit [rad/s^2]
        pos_max         % upper speed limits (9x1) [rad/s]
        pos_min         % lower speed limits (9x1) [rad/s]
        pos             % current rotor speeds (9x1) [rad/s]
    end

    methods
        function obj = EngineDynamics(dt)
            obj.dt         = dt;
            obj.wn         = 4 * pi;                  % rad/s (setupEngines)
            obj.rate_limit = 100.0;
            obj.pos_max    = [1600 * ones(8, 1); 2000] .* (2 * pi / 60); % rad/s
            obj.pos_min    = zeros(9, 1);
            obj.pos        = zeros(9, 1);
        end

        function reset(obj, pos0)
            % Initialize rotor speeds (e.g., at trim).
            obj.pos = min(max(pos0(:), obj.pos_min), obj.pos_max);
        end

        function pos = step(obj, cmd)
            % Advance the engine servo bank one time step toward cmd (9x1, rad/s).
            rate    = obj.wn .* (cmd(:) - obj.pos);
            rate    = min(max(rate, -obj.rate_limit), obj.rate_limit);
            obj.pos = min(max(obj.pos + obj.dt .* rate, obj.pos_min), obj.pos_max);
            pos     = obj.pos;
        end
    end
end
