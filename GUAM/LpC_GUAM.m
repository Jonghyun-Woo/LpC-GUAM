classdef LpC_GUAM < handle
    % LPC_GUAM Top-level orchestrator for the LpC m-code refactor.
    %
    % Flat-earth (NED) closed-loop simulation of the GUAM Lift+Cruise
    % hover-to-cruise transition. No ECI/ECEF frames — position is NED,
    % altitude = -rd. The reference trajectory is selected by scenario
    % (see ReferenceTrajectory / SimConfig):
    %   'althold' (default) - hover climb to 80 ft, then hold 80 ft in cruise
    %   'climb'             - hover climb to 80 ft, then climb to 100 ft
    %
    % Per-step pipeline (mirrors GUAM.slx with the default variants:
    % Polynomial aero, FirstOrder actuators, Simple EOM):
    %   RSLQR control -> engine/surface servo dynamics -> standard
    %   atmosphere -> polynomial aero-propulsive forces/moments ->
    %   6-DOF rigid body EOM -> forward-Euler state update.
    properties
        config          % Config (central hub; sub-configs distributed below)
        vehicleConfig   % VehicleConfig (also passed to the aero model as Model)
        simConfig       % SimConfig (owns dt/T and the scenario)
        refTraj         % Reference trajectory table (from ControllerConfig.getReferenceTrajectory)
        units           % Units ('ft','slug') — aero model unit conversions

        rigidBody       % RBD 6-DOF equations of motion
        aeroFrame       % AeroFrame (alpha/beta/airspeed for logging)
        environment     % Environment (flat-earth atmosphere + steady wind)
        controller      % Controller (RSLQR baseline + liveness safety filter)
        engineDynamics  % EngineDynamics (9 rotor speed servos)
        surfaceDynamics % SurfaceDynamics (5 surface servos [LA RA LE RE RUD])

        state           % Current state [rn re rd u v w phi theta psi p q r]' (12x1)
        time            % Current simulation time [s]
    end

    methods
        function obj = LpC_GUAM(cfg)
            % cfg : Config hub (see config/Config.m). Sub-configs are
            % distributed to the owning components below.
            if nargin < 1 || isempty(cfg), cfg = Config('althold'); end

            obj.config          = cfg;
            obj.vehicleConfig   = cfg.vehicle;
            obj.simConfig       = cfg.sim;
            obj.refTraj         = cfg.controller.getReferenceTrajectory();
            obj.units           = Units('ft', 'slug');

            obj.rigidBody       = RBD(cfg.vehicle);
            obj.aeroFrame       = AeroFrame();
            obj.environment     = Environment();
            obj.controller      = Controller(cfg.controller, cfg.sim.dt);
            obj.engineDynamics  = EngineDynamics(cfg.sim.dt);
            obj.surfaceDynamics = SurfaceDynamics(cfg.sim.dt);

            obj.reset();
        end

        function reset(obj)
            % Initialize the vehicle at the trim condition of the first
            % reference-trajectory point (as GUAM's setupTrim does for the
            % initial reference velocity).
            % The trim state and the trim actuator positions both come from the
            % controller, which assembles the surface vector through the same
            % total_cmd() path the control law uses. Doing it here as well used
            % to duplicate that [flap-ail; flap+ail; ele; ele; rud] ordering, so
            % a change on one side would silently disagree with the other.
            [x0, engine0, surface0] = obj.controller.initial_condition(obj.refTraj);

            obj.state = x0;
            obj.engineDynamics.reset(engine0);
            obj.surfaceDynamics.reset(surface0);
            obj.controller.reset();
            obj.time = 0;
        end

        function [engine, surface, Fb, Mb] = step(obj, ref)
            % Advance the closed-loop simulation one time step.
            % ref : struct with fields pos (3x1 NED), vel (3x1 heading
            %       frame), chi, chi_dot

            % 1. Control law -> effector commands
            [eng_cmd, srf_cmd] = obj.controller.control(obj.state, ref);

            % 2. Actuator dynamics
            engine  = obj.engineDynamics.step(eng_cmd);
            surface = obj.surfaceDynamics.step(srf_cmd);

            % 3. Atmosphere at current altitude (flat earth: h = -rd)
            [rho, a] = obj.environment.atmosphere(-obj.state(3));

            % 4. Aero-propulsive forces/moments (air-relative body velocity)
            R_i2b = RSLQR.rotm_i2b(obj.state(7), obj.state(8), obj.state(9));
            v_air = obj.environment.airspeed_body(obj.state(4:6), R_i2b);
            x_aero = [v_air; obj.state(10:12)];
            [Fb, Mb] = run_LpC_aero(x_aero, engine, surface, rho, a, ...
                                    obj.units, obj.vehicleConfig);

            % 5. 6-DOF EOM (gravity added inside RBD) + forward Euler
            dx = obj.rigidBody.calculate_dynamics(obj.state, Fb, Mb);
            obj.state = obj.state + obj.simConfig.dt .* dx;
            obj.time  = obj.time + obj.simConfig.dt;
        end

        function s = saveState(obj)
            % SAVESTATE snapshot of every piece of memory the closed loop
            % carries. Needed by the reference governor, which rolls the
            % whole loop forward as a "what-if" prediction and then rewinds.
            %
            % Missing ANY of these makes the prediction disagree with what
            % the vehicle would actually do:
            %   state      : rigid-body 12 states
            %   eng/srf    : first-order servo positions (actuators lag)
            %   ctrl_*_i   : RSLQR servo-compensator integrators. These are
            %                wound up by past error; starting a rollout with
            %                them zeroed under-predicts the acceleration and
            %                therefore under-predicts the overshoot.
            %   phi_v/theta_v : virtual attitude effectors fed back by the
            %                allocation
            s = struct( ...
                'state',      obj.state, ...
                'time',       obj.time, ...
                'eng_pos',    obj.engineDynamics.pos, ...
                'srf_pos',    obj.surfaceDynamics.pos, ...
                'ctrl_lon_i', obj.controller.baseline_controller.ctrl_lon_i, ...
                'ctrl_lat_i', obj.controller.baseline_controller.ctrl_lat_i, ...
                'phi_v',      obj.controller.baseline_controller.phi_v, ...
                'theta_v',    obj.controller.baseline_controller.theta_v);
        end

        function restoreState(obj, s)
            % RESTORESTATE rewind to a snapshot taken by saveState.
            obj.state                 = s.state;
            obj.time                  = s.time;
            obj.engineDynamics.pos    = s.eng_pos;
            obj.surfaceDynamics.pos   = s.srf_pos;
            bc = obj.controller.baseline_controller;
            bc.ctrl_lon_i = s.ctrl_lon_i;
            bc.ctrl_lat_i = s.ctrl_lat_i;
            bc.phi_v      = s.phi_v;
            bc.theta_v    = s.theta_v;
        end
    end
end
