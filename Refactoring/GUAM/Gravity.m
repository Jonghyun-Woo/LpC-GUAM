classdef Gravity < handle
    properties
        % g = 32.17405284; % Standard gravity (ft/s^2)
        g
    end
    
    methods
        function obj = Gravity(vehicleConfig)
            obj.g = vehicleConfig.gravity;
        end
        
        function gravity_b = compute_gravity_body_frame(obj, state)
            % Compute gravity vector in body frame based on current orientation
            % theta: Pitch angle (radians)
            % phi: Roll angle (radians)

            % theta = state();
            % phi = state();

            % sThe = sin(theta);
            % cThe = cos(theta);
            % sPhi = sin(phi);
            % cPhi = cos(phi);
            
            % gravity_b = obj.g * [-sThe; sPhi*cThe; cPhi*cThe];
            gravity_b = zeros([3, 1]);
        end
    end
end