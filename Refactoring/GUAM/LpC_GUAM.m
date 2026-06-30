classdef LpC_GUAM < handle
    % LPC_GUAM Top-level orchestrator for the LpC m-code refactor.
    % Wires VehicleConfig/SimConfig into RigidBody6Dof and integrates the
    % 12-state EOM forward in time using forward Euler.
    properties
        vehicleConfig   % VehicleConfig
        simConfig       % SimConfig
        rigidBody       % RigidBody6Dof
        aeroFrame       % AeroFrame
        engineDynamics  % Engine Dynamics
        surfaceDynamics % Control Surface Dynamics (Aileron, Elevator, Rudder)
        
        state           % Current state vector [rn re rd u v w phi theta psi p q r]' (12x1)
    end

    methods
        function obj = LpC_GUAM()
            obj.vehicleConfig = VehicleConfig();
            obj.simConfig = SimConfig();
            obj.rigidBody = RBD(obj.vehicleConfig);
            obj.aeroFrame = AeroFrame();
            obj.state = zeros(12, 1);
        end

        function step(obj, engine_cmd, surface_cmd)
            
            % engine_out
            % surface_out
            % run_LpC_aero(x,n,d,rho,a,Units,Model);
            run_LpC_aero(x, engine_out, surface_out, rho, a, units, obj.vehicleConfig);
        end
    end
end
