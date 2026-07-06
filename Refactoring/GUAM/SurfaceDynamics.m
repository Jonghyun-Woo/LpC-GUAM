classdef SurfaceDynamics < handle
    % SURFACEDYNAMICS First-order rate- and position-limited servo bank
    % for the 5 control surfaces [LA RA LE RE RUD] (flaperon convention:
    % left/right ailerons carry both flap and aileron deflection).
    %
    % Mirrors GUAM's ServosLib "First Order Dynamic Rate and Position
    % Limited Servo" (default actType = ActuatorEnum.FirstOrder):
    %   pos_dot = sat(wn*(cmd - pos), +/-rate_limit), pos clamped to limits
    % Parameters from vehicles/Lift+Cruise/setup/setupActuators.m.
    properties
        dt              % integration time step [s]
        wn              % servo bandwidth [rad/s]
        rate_limit      % surface rate limit [rad/s]
        pos_max         % upper position limits (5x1) [rad]
        pos_min         % lower position limits (5x1) [rad]
        pos             % current surface positions [LA RA LE RE RUD]' [rad]
    end

    methods
        function obj = SurfaceDynamics(dt)
            obj.dt         = dt;
            obj.wn         = 20.0;                    % rad/s (setupActuators)
            obj.rate_limit = 100.0;                   % rad/s
            obj.pos_max    =  deg2rad(30) * ones(5, 1);
            obj.pos_min    = -deg2rad(30) * ones(5, 1);
            obj.pos        = zeros(5, 1);
        end

        function reset(obj, pos0)
            % Initialize servo states (e.g., at trim surface deflections).
            obj.pos = min(max(pos0(:), obj.pos_min), obj.pos_max);
        end

        function pos = step(obj, cmd)
            % Advance the servo bank one time step toward cmd (5x1, rad).
            rate    = obj.wn .* (cmd(:) - obj.pos);
            rate    = min(max(rate, -obj.rate_limit), obj.rate_limit);
            obj.pos = min(max(obj.pos + obj.dt .* rate, obj.pos_min), obj.pos_max);
            pos     = obj.pos;
        end
    end
end
