classdef VehicleConfig
    % VEHICLECONFIG Configuration class for vehicle parameters
    properties
        mass        % Mass of the vehicle (slug)
        inertia     % Inertia matrix (3x3, slug*ft^2)
        inv_inertia % Inverse of inertia matrix (3x3, 1/(slug*ft^2))
        gravity     % Gravity acceleration (ft/s^2)

        engine_num      % number of engines
        Prop_location   % propeller/rotor locations (ft)
        Prop_Diameter   % propeller/rotor diameter (ft)
        Prop_inertia    % Propeller/rotor rotational moment of inertia (slug*ft^2)
        Prop_angles     % propeller/rotor orientation angles (roll-pitch-yaw) [deg]
        % rho             % Air density (slug/ft^3)

    end
    
    methods
        function obj = VehicleConfig()
            obj.mass    = 181.789249;                                       % Mass of the LpC GUAM (slug)
            obj.inertia = diag([13051.74318, 16660.75897, 24735.13582]);    % Inertia matrix (3x3, slug*ft^2)
            obj.inv_inertia = (obj.inertia)\eye(3);                         % Inverse of inertia matrix (3x3, 1/(slug*ft^2))
            obj.gravity = 32.174;                                           % Gravity acceleration (ft/s^2)
            % obj.rho = 0.0023769; % Air density at sea level (slug/ft^3)

            obj.engine_num = 9;                                             % number of engines
            obj.Prop_location = ...                                         % propeller/rotor locations (ft)         
                [ -5.07,  -4.63, -4.63,  -5.07,   -19.2,  -18.76, -18.76,  -19.2, -31.94;...
                -18.750,  -8.45,  8.45, 18.750, -18.750,   -8.45,   8.45, 18.750,  0.000;...
                  -6.73,  -7.04, -7.04,  -6.73,   -9.01,    -9.3,   -9.3,  -9.01,  -7.79];
            obj.Prop_Diameter = [10*ones(1,8),9];                           % propeller/rotor diameter (ft)
            obj.Prop_inertia = [13.486*ones(1,8),17.486];                   % propeller/rotor rotational moment of inertia (slug*ft^2)
            obj.Prop_angles =...                                            % propeller/rotor orientation angles (roll-pitch-yaw) [deg] 
                [ 0,   0,   0,   0,   0,   0,   0,   0,   0;...
                 90, -90, -90, -90, -90, -90, -90, -90,   0;...
                  0,  +8,  -8,   0,   0,  +8,  -8,   0,   0];
        end
    end
end

