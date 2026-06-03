classdef RigidBody6Dof < handle
    % RIGIDBODY6DOF 6-Degree-of-Freedom Rigid Body Dynamics Model
    % Uses inverse inertia matrix formulation for rotational dynamics.
    properties
        m       % Mass (kg)
        I       % Inertia Matrix (3x3, kg*m^2)
        invI    % Inverse Inertia Matrix (3x3, 1/(kg*m^2))
        g       % Gravity acceleration (m/s^2)
    end
    
    methods
        function obj = RigidBody6Dof(vehicleConfig)
            % Constructor
            if nargin < 3
                gravity = 9.80665;
            end
            obj.m = vehicleConfig.mass;
            obj.I = vehicleConfig.inertia;
            obj.invI = inv(obj.I);
            obj.g = vehicleConfig.gravity;
        end
        
        function dx = calculate_dynamics(obj, x, F, M)
            % x: State vector [rn, re, rd, u, v, w, phi, theta, psi, p, q, r]' (12x1)
            % F: Body forces [Fx, Fy, Fz]' (3x1)
            % M: Body moments [Mx, My, Mz]' (3x1)
            % dx: State derivative (12x1)
            
            % 1. State Extraction
            % Position (not directly needed for dynamics, but for completeness)
            rn = x(1); re = x(2); rd = x(3);
            
            % Linear Velocity (Body frame)
            V_b = x(4:6);
            u = V_b(1); v = V_b(2); w = V_b(3);
            
            % Euler Angles
            phi = x(7); theta = x(8); psi = x(9);
            
            % Angular Velocity (Body frame)
            omega = x(10:12);
            p = omega(1); q = omega(2); r = omega(3);
            
            % Pre-compute trigonometric functions
            cPhi = cos(phi); sPhi = sin(phi);
            cThe = cos(theta); sThe = sin(theta); tThe = tan(theta);
            cPsi = cos(psi); sPsi = sin(psi);
            
            % 2. Translational Kinematics (Position Derivative)
            R_b2i = [
                cThe*cPsi, sPhi*sThe*cPsi - cPhi*sPsi, cPhi*sThe*cPsi + sPhi*sPsi;
                cThe*sPsi, sPhi*sThe*sPsi + cPhi*cPsi, cPhi*sThe*sPsi - sPhi*cPhi;
                -sThe,     sPhi*cThe,                  cPhi*cThe
            ];
            p_dot = R_b2i * V_b;
            
            % 3. Rotational Kinematics (Euler Angle Derivative)
            T_rot = [
                1, sPhi*tThe, cPhi*tThe;
                0, cPhi,      -sPhi;
                0, sPhi/cThe, cPhi/cThe
            ];
            euler_dot = T_rot * omega;
            
            % 4. Translational Dynamics (Linear Velocity Derivative)
            gravity_b = obj.g * [-sThe; sPhi*cThe; cPhi*cThe];
            cross_omega_v = cross(omega, V_b);
            v_dot = -cross_omega_v + gravity_b + (1/obj.m) * F;
            
            % 5. Rotational Dynamics (Angular Velocity Derivative)
            % Form: w_dot = inv(I) * (M - cross(w, I*w))
            % Using mldivide (\) for numerical stability and speed
            cross_w_Iw = cross(omega, obj.I * omega);
            omega_dot = obj.I \ (M - cross_w_Iw);
            
            % 6. Assemble Output
            dx = [p_dot; v_dot; euler_dot; omega_dot];
        end
    end
end