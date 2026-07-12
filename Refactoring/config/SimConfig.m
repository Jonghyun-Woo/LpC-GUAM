classdef SimConfig < handle
    properties
        % Simulation parameters
        dt        % Time step (s)
        T         % Total simulation time (s)
        scenario  % Reference trajectory scenario: 'althold' (default) | 'climb'
                  %   | 'lon_brt_verify' (verification-only descending WH3 scenario)
    end

    methods
        function obj = SimConfig(scenario)
            obj.dt = 0.01;  % Time step (s)
            obj.T = 40;     % Total simulation time (s)
            if nargin < 1 || isempty(scenario)
                scenario = 'althold';
            end
            obj.scenario = scenario;
        end

        function ref = getReferenceTrajectory(obj)
            % Load the reference trajectory table for the configured scenario.
            ref = ReferenceTrajectory.build(obj.scenario, obj.dt, obj.T);
        end

        function display_simulation_config(obj)
            fprintf('Simulation Configuration:\n');
            fprintf('Time step: %.2f s\n', obj.dt);
            fprintf('Total simulation time: %.2f s\n', obj.T);
            fprintf('Scenario: %s\n', obj.scenario);
        end
    end
end