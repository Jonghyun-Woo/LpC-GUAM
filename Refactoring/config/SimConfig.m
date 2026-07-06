classdef SimConfig < handle
    properties
        % Simulation parameters
        dt  % Time step (s)
        T   % Total simulation time (s)
    end
    
    methods
        function obj = SimConfig()
            obj.dt = 0.01;  % Time step (s)
            obj.T = 40;     % Total simulation time (s)
        end
        
        function display_simulation_config(obj)
            fprintf('Simulation Configuration:\n');
            fprintf('Time step: %.2f s\n', obj.dt);
            fprintf('Total simulation time: %.2f s\n', obj.T);
        end
    end
end