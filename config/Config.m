classdef Config < handle
    % Central configuration hub for the LpC GUAM simulation.
    %
    % Individual config classes stay in their own files for readability, but
    % assembly/ownership/passing is centralized here. The main script builds a
    % single Config and hands it to LpC_GUAM(cfg).
    %
    % Layout (2-tier):
    %   cfg.sim        : SimConfig        (dt, M, T, scenario)
    %   cfg.vehicle    : VehicleConfig    (mass/geometry)
    %   cfg.controller : ControllerConfig (target_vel, refTraj, rslqr, filter)
    %
    % overrides (optional struct) fields and destinations:
    %   .M, .dt          -> SimConfig
    %   .target_vel      -> ControllerConfig
    %   scenario (arg)   -> SimConfig (mission profile)

    properties
        sim         % SimConfig
        vehicle     % VehicleConfig
        controller  % ControllerConfig
    end

    methods
        function obj = Config(scenario, overrides)
            if nargin < 1 || isempty(scenario), scenario = 'althold'; end
            if nargin < 2 || isempty(overrides), overrides = struct(); end

            % Build order matters: controller depends on sim.dt/T/scenario.
            obj.sim        = SimConfig(scenario, overrides);
            obj.vehicle    = VehicleConfig();
            obj.controller = ControllerConfig(obj.sim.scenario, ...
                                              obj.sim.dt, obj.sim.T, overrides);
        end
    end
end
