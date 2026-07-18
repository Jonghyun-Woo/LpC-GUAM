classdef SimConfig < handle
    properties
        % Simulation parameters
        steps         % Default number of simulation steps
        dt        % Time step (s)
        T         % Total simulation time (s)
        scenario  % Reference trajectory scenario: 'althold' (default) | 'climb'
        %   | 'lon_brt_verify' (verification-only descending WH3 scenario)
    end

    methods
        function obj = SimConfig(scenario, overrides)
            if nargin < 1 || isempty(scenario), scenario = 'althold'; end
            if nargin < 2 || isempty(overrides), overrides = struct(); end
            obj.dt       = getfield_default(overrides, 'dt', 0.01);
            obj.steps    = getfield_default(overrides, 'steps', 4000);
            obj.T        = obj.steps * obj.dt;
            obj.scenario = scenario;
        end

        function display_simulation_config(obj)
            fprintf('Simulation Configuration:\n');
            fprintf('Number of steps: %d\n', obj.steps);
            fprintf('Time step: %.2f s\n', obj.dt);
            fprintf('Total simulation time: %.2f s\n', obj.T);
            fprintf('Scenario: %s\n', obj.scenario);
        end
    end
end