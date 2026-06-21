classdef LpC_GUAM < handle
    % LPC_GUAM Top-level orchestrator for the LpC m-code refactor.
    % Wires VehicleConfig/SimConfig into RigidBody6Dof and integrates the
    % 12-state EOM forward in time using forward Euler.
    properties
        vehicleConfig   % VehicleConfig
        simConfig       % SimConfig
        rigidBody       % RigidBody6Dof
        aeroFrame       % AeroFrame
        state           % Current state vector [rn re rd u v w phi theta psi p q r]' (12x1)
        history         % Struct array log: t, state, alpha, beta
    end

    methods
        function obj = LpC_GUAM()
            obj.vehicleConfig = VehicleConfig();
            obj.simConfig = SimConfig();
            obj.rigidBody = RigidBody6Dof(obj.vehicleConfig);
            obj.aeroFrame = AeroFrame();
            obj.state = zeros(12, 1);
            obj.history = struct('t', {}, 'state', {}, 'alpha', {}, 'beta', {});
        end

        function run(obj)
            t_vec = 0:obj.simConfig.dt:obj.simConfig.T;
            obj.history = struct('t', {}, 'state', {}, 'alpha', {}, 'beta', {});

            for k = 1:length(t_vec)
                u = obj.state(4);
                v = obj.state(5);
                w = obj.state(6);
                [alpha, beta, ~] = obj.aeroFrame.compute(u, v, w);

                % Placeholder body forces/moments until the polynomial
                % aero/propulsion model is ported (follow-up task).
                F = zeros(3, 1);
                M = zeros(3, 1);

                obj.history(end+1) = struct('t', t_vec(k), 'state', obj.state, 'alpha', alpha, 'beta', beta);

                dx = obj.rigidBody.calculate_dynamics(obj.state, F, M);
                obj.state = obj.state + obj.simConfig.dt * dx;
            end
        end
    end
end
