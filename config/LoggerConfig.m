classdef LoggerConfig < handle
    % Logging / plotting / figure-saving options for SimLogger.
    %
    % handle class: assembled on the central hub (cfg.logger) alongside
    % sim/vehicle/controller. dt and step count are NOT stored here — they are
    % owned by SimConfig and read by SimLogger at construction (ADR-0003).
    %
    % Simplifications (fixed, not configurable): figures are always saved as
    % PNG, always created Visible, and every step is logged (no decimation).
    %
    % overrides (optional struct) fields map 1:1 to the properties below and
    % are shared with the other sub-configs through the same overrides struct
    % (keys here do not collide with SimConfig/Controller/Filter keys).

    properties
        enable      = true;   % master logging on/off
        logFilter   = true;   % collect filter/BRT diagnostics (harmless if filter off)
        plotBasic   = true;   % position / velocity / attitude / effector figures
        plotFilter  = true;   % single-run BRT value(t) + nominal vs filtered input
        saveFigures = false;  % save figures (PNG)
        saveData    = false;  % save buffers / exported trace as .mat
        saveDir     = '';     % empty -> SimLogger picks logs/<timestamp> at runtime
    end

    methods
        function obj = LoggerConfig(overrides)
            if nargin < 1 || isempty(overrides), overrides = struct(); end
            obj.enable      = getfield_default(overrides, 'enable',      obj.enable);
            obj.logFilter   = getfield_default(overrides, 'logFilter',   obj.logFilter);
            obj.plotBasic   = getfield_default(overrides, 'plotBasic',   obj.plotBasic);
            obj.plotFilter  = getfield_default(overrides, 'plotFilter',  obj.plotFilter);
            obj.saveFigures = getfield_default(overrides, 'saveFigures', obj.saveFigures);
            obj.saveData    = getfield_default(overrides, 'saveData',    obj.saveData);
            obj.saveDir     = getfield_default(overrides, 'saveDir',     obj.saveDir);
        end
    end
end
