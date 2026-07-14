classdef Environment < handle
    % ENVIRONMENT Flat-earth environment model (NED frame only).
    %
    % Replaces GUAM's Environment.slx for the m-code refactor: no ECI/ECEF
    % frames are modeled — position is NED over a flat earth and altitude
    % is h = -rd. Provides:
    %   - US Standard Atmosphere 1976 troposphere (English units), matching
    %     the sea-level constants in setup/setupAtmosphere.m
    %     (rho0 = 0.0023769 slug/ft^3, a0 = 1116.45 ft/s)
    %   - steady wind in the NED frame (default: calm; turbulence not modeled)
    properties
        wind_ned = zeros(3, 1);  % steady wind velocity in NED [ft/s]
        deltaT   = 0;            % temperature deviation from standard [degR]
    end

    properties (Constant)
        T0    = 518.67;      % sea-level temperature [degR]
        P0    = 2116.22;     % sea-level pressure [lbf/ft^2]
        L     = 0.00356616;  % tropospheric lapse rate [degR/ft]
        Rgas  = 1716.49;     % gas constant for air [ft*lbf/(slug*degR)]
        gamma = 1.4;         % ratio of specific heats
        g0    = 32.17405;    % standard gravity [ft/s^2]
    end

    methods
        function [rho, a, T, P] = atmosphere(obj, h)
            % Standard-atmosphere properties at geometric altitude h [ft]
            % (troposphere model; h clamped at sea level).
            h   = max(h, 0);
            T   = obj.T0 - obj.L * h + obj.deltaT;
            P   = obj.P0 * ((obj.T0 - obj.L * h) / obj.T0)^(obj.g0 / (obj.L * obj.Rgas));
            rho = P / (obj.Rgas * T);
            a   = sqrt(obj.gamma * obj.Rgas * T);
        end

        function v_air_b = airspeed_body(obj, v_b, R_i2b)
            % Air-relative velocity in the body frame.
            % v_b   : inertial velocity in body frame (3x1) [ft/s]
            % R_i2b : NED-to-body rotation matrix (3x3)
            v_air_b = v_b - R_i2b * obj.wind_ned;
        end
    end
end
