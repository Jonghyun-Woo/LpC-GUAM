classdef ControllerConfig < handle
    % Controller-side configuration bundle.
    %
    % Holds the mission characteristics common to any controller (reference
    % trajectory and its parameters) plus the per-controller config objects.
    % Future controllers add their own config alongside `rslqr` while sharing
    % `refTraj`/`target_vel`.
    %
    % scenario / dt / T are owned by SimConfig and injected here only as build
    % inputs for the reference trajectory (not stored as properties).

    properties
        target_vel = 15;   % cruise forward speed [ft/s] (mission parameter)
        refTraj            % built reference trajectory table
        rslqr              % RSLQRConfig (gains/limits)
        filter             % FilterConfig (liveness filter config)
    end

    methods
        function obj = ControllerConfig(scenario, dt, T, overrides)
            if nargin < 4 || isempty(overrides), overrides = struct(); end
            obj.target_vel = getfield_default(overrides, 'target_vel', obj.target_vel);
            obj.rslqr      = RSLQRConfig();
            obj.filter     = FilterConfig(overrides);
            obj.refTraj    = ReferenceTrajectory.build(scenario, dt, T, obj.target_vel);
        end

        function ref = getReferenceTrajectory(obj)
            ref = obj.refTraj;
        end
    end
end
