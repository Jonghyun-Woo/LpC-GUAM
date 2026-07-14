classdef AeroFrame < handle
    % AEROFRAME Computes aero-frame quantities (alpha, beta, airspeed)
    % from body-frame velocity, matching the legacy polynomial aero model
    % convention (vehicles/Lift+Cruise/AeroProp/Polynomial/LpC_aero_p_v2.m).
    properties (Constant, Access = private)
        u_min = 0.01; % below this body-axis u (ft/s), alpha/beta default to zero
    end

    methods
        function [alpha, beta, V] = compute(obj, u, v, w)
            % state: body-frame velocity components u, v, w (ft/s)
            % alpha, beta: angle of attack, angle of sideslip (rad)
            % V: true airspeed (ft/s)
            V = norm([u, v, w]);
            if u > obj.u_min
                alpha = atan(w / u);
                beta = asin(v / max(V, eps));
            else
                alpha = 0;
                beta = 0;
            end
        end
    end
end
