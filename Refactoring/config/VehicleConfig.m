classdef VehicleConfig
    % VEHICLECONFIG Configuration class for vehicle parameters
    properties
        mass        % Mass of the vehicle (lbs)
        inertia     % Inertia matrix (3x3, lbs*ft^2)
        inv_inertia % Inverse of inertia matrix (3x3, 1/(lbs*ft^2))
        gravity     % Gravity acceleration (ft/s^2)

        engine_num      % number of engines
        Prop_location   % propeller/rotor locations (ft)
        Prop_Diameter   % propeller/rotor diameter (ft)
        Prop_angles     % propeller/rotor orientation angles (roll-pitch-yaw) [deg]
    end
    
    methods
        function obj = VehicleConfig()
            obj.mass    = 181.789249;
            obj.inertia = diag([13051.74318, 16660.75897, 24735.13582]);
            obj.inv_inertia = (obj.inertia)\eye(3);
            obj.gravity = 32.174;

            obj.Prop_location = ...
                [ -5.07,  -4.63, -4.63,  -5.07,   -19.2,  -18.76, -18.76,  -19.2, -31.94;...
                -18.750,  -8.45,  8.45, 18.750, -18.750,   -8.45,   8.45, 18.750,  0.000;...
                  -6.73,  -7.04, -7.04,  -6.73,   -9.01,    -9.3,   -9.3,  -9.01,  -7.79];
            obj.Prop_Diameter = [10*ones(1,8),9];
            obj.Prop_angles =...
                [ 0,   0,   0,   0,   0,   0,   0,   0,   0;...
                 90, -90, -90, -90, -90, -90, -90, -90,   0;...
                  0,  +8,  -8,   0,   0,  +8,  -8,   0,   0];
        end
    end
end

