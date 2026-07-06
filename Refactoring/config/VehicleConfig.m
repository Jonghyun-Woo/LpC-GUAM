classdef VehicleConfig
    % Vehicle configuration parameters for the Lift+Cruise (LpC) GUAM.
    % Source: vehicles/Lift+Cruise/setup/LpC_model_parameters.m
    %         SACD Lift+Cruise reference configuration (publicly releasable)
    %
    % Ref [1]: https://sacd.larc.nasa.gov/uam-refs/  (list-cruiseTE6.list)
    % Ref [2]: lift-cruiseTE6.xlsx

    % -----------------------------------------------------------------------
    % Directly specified parameters — treat as fixed for this configuration.
    % -----------------------------------------------------------------------
    properties (Constant)

        % --- Mass properties ---
        mass  = 181.789249                                          % total vehicle mass [slug]
        I     = diag([13051.74318, 16660.75897, 24735.13582])      % inertia tensor (3x3) [slug·ft²]
        cm_b  = [-13.841959237973750, ...
                  -0.000000144320024, ...
                  -4.603304281304035]                               % CG location (+fwd,+right,+down) [ft]

        % --- Aerodynamic reference geometry ---
        S     = 186.0    % wing reference area [ft²]
        cbar  = 3.18     % mean aerodynamic chord [ft]
        b     = 47.5     % wingspan [ft]

        % --- Propulsion geometry ---
        Prop_location = ...                                         % rotor hub locations [ft]
              [ -5.07,  -4.63,  -4.63,  -5.07, -19.2,  -18.76, -18.76, -19.2,  -31.94; ...
               -18.750,  -8.45,   8.45, 18.750,-18.750,  -8.45,   8.45, 18.750,  0.000; ...
                -6.73,   -7.04,  -7.04,  -6.73,  -9.01,  -9.3,   -9.3,  -9.01,  -7.79]

        Prop_D    = [10*ones(1,8), 9]                              % propeller/rotor diameter [ft]
        prop_spin = [-1, 1, -1, 1, 1, -1, 1, -1, 1]              % rotation direction (CW=+1, CCW=-1)
        Ip        = [13.486*ones(1,8), 17.486]                    % rotational moment of inertia [slug·ft²]

        Prop_angles = [  0,   0,   0,   0,   0,   0,   0,   0,   0; ...
                       -90, -90, -90, -90, -90, -90, -90, -90,   0; ...
                         0,  +8,  -8,   0,   0,  +8,  -8,   0,   0]   % orientation (roll-pitch-yaw) [deg]

        % --- Physical constant ---
        gravity = 32.174    % gravitational acceleration [ft/s²]

    end

    % -----------------------------------------------------------------------
    % Derived quantities — computed once in the constructor from Constant values.
    % -----------------------------------------------------------------------
    properties
        inv_I             % inverse of inertia tensor (3x3) [1/(slug·ft²)]
        Prop_R_BH         % rotation matrix body→hub frame (3x3x9)
        Prop_rot_axis_e   % unit vector of rotor axis of rotation in body frame (3x9)
        p_T_e             % axial thrust direction unit vector in body frame (3x9)
    end

    methods
        function obj = VehicleConfig()
            obj.inv_I = obj.I \ eye(3);

            numEngines          = size(obj.Prop_location, 2);
            obj.Prop_R_BH       = zeros(3, 3, numEngines);
            obj.Prop_rot_axis_e = zeros(3, numEngines);

            for ii = 1:numEngines
                phi   = obj.Prop_angles(1, ii) * pi / 180;
                theta = obj.Prop_angles(2, ii) * pi / 180;
                psi   = obj.Prop_angles(3, ii) * pi / 180;

                % Rotation matrix from body frame to hub frame
                R_HB = [cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta); ...
                        cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi); ...
                       -sin(theta),          cos(theta)*sin(phi),                               cos(phi)*cos(theta)];

                obj.Prop_R_BH(:, :, ii)   = R_HB';
                obj.Prop_rot_axis_e(:, ii) = R_HB' * [1; 0; 0];
            end

            obj.p_T_e = obj.Prop_rot_axis_e;
        end
    end
end