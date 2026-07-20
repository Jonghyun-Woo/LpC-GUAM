classdef LpC_GUAM < handle
    % LPC_GUAM Flat-earth (NED) 6-DOF plant for the GUAM Lift+Cruise model.
    %
    % Pure plant. Given per-step actuator commands (rotor speeds + surface
    % deflections), it advances its own state through first-order actuator
    % servos, the polynomial aero-propulsive model, and RK4 integration of the
    % rigid-body EOM. No ECI/ECEF frames — position is NED, altitude = -rd.
    %
    % The controller lives outside the plant: the simulation driver calls
    % controller.control() to obtain the effector commands, then passes them to
    % step(). step() returns nothing — its effect is the update of obj.state /
    % obj.time (and the actual actuator positions obj.engine / obj.surface).
    %
    % Per-step pipeline (mirrors GUAM.slx default variants: Polynomial aero,
    % FirstOrder actuators, Simple EOM):
    %   engine/surface servo dynamics -> standard atmosphere ->
    %   polynomial aero-propulsive forces/moments -> 6-DOF rigid body EOM ->
    %   RK4 state update.
    properties
        vehicleConfig   % VehicleConfig (also passed to the aero model as Model)
        simConfig       % SimConfig (owns dt/T and the scenario)
        units           % Units ('ft','slug') — aero model unit conversions

        rigidBody       % RBD 6-DOF equations of motion
        aeroFrame       % AeroFrame (alpha/beta/airspeed for logging)
        environment     % Environment (flat-earth atmosphere + steady wind)
        engineDynamics  % EngineDynamics (9 rotor speed servos)
        surfaceDynamics % SurfaceDynamics (5 surface servos [LA RA LE RE RUD])

        state           % Current state [rn re rd u v w phi theta psi p q r]' (12x1)
        time            % Current simulation time [s]
        engine          % Actual rotor speeds after servo dynamics (9x1) [rad/s]
        surface         % Actual surface deflections after servo dynamics (5x1) [rad]
    end

    methods
        function obj = LpC_GUAM(cfg)
            % cfg : Config hub (see config/Config.m). Only the vehicle/sim
            % sub-configs are used here; the controller is built and owned by
            % the simulation driver, not the plant.
            if nargin < 1 || isempty(cfg), cfg = Config('althold'); end

            obj.vehicleConfig   = cfg.vehicle;
            obj.simConfig       = cfg.sim;
            obj.units           = Units('ft', 'slug');

            obj.rigidBody       = RBD(cfg.vehicle);
            obj.aeroFrame       = AeroFrame();
            obj.environment     = Environment();
            obj.engineDynamics  = EngineDynamics(cfg.sim.dt);
            obj.surfaceDynamics = SurfaceDynamics(cfg.sim.dt);

            obj.state   = zeros(12, 1);
            obj.time    = 0;
            obj.engine  = obj.engineDynamics.pos;
            obj.surface = obj.surfaceDynamics.pos;
        end

        function reset(obj, x0, engine0, surface0)
            % Initialize the plant state and actuator servos. The trim initial
            % condition (rigid state x0 + actuator trims engine0/surface0) is
            % supplied by the driver from the controller (see
            % Controller.initial_condition).
            obj.state = x0(:);
            obj.time  = 0;
            obj.engineDynamics.reset(engine0);
            obj.surfaceDynamics.reset(surface0);
            obj.engine  = obj.engineDynamics.pos;
            obj.surface = obj.surfaceDynamics.pos;
        end

        function step(obj, engine_cmd, surface_cmd)
            % Advance the plant one time step under the given actuator commands.
            % Updates obj.state / obj.time and the actual actuator positions
            % obj.engine / obj.surface. Returns nothing.
            %
            % engine_cmd  : 9x1 commanded rotor speeds [rad/s]
            % surface_cmd : 5x1 commanded surface deflections [rad]
            dt = obj.simConfig.dt;

            % 1. Actuator servo dynamics (advanced once; held over the RK4 stages)
            engine  = obj.engineDynamics.step(engine_cmd);
            surface = obj.surfaceDynamics.step(surface_cmd);
            obj.engine  = engine;
            obj.surface = surface;

            % 2. RK4 integration of the 6-DOF rigid-body state
            x  = obj.state;
            k1 = obj.state_derivative(x,             engine, surface);
            k2 = obj.state_derivative(x + 0.5*dt*k1, engine, surface);
            k3 = obj.state_derivative(x + 0.5*dt*k2, engine, surface);
            k4 = obj.state_derivative(x +     dt*k3, engine, surface);
            obj.state = x + (dt / 6) .* (k1 + 2*k2 + 2*k3 + k4);
            obj.time  = obj.time + dt;
        end

        function dx = state_derivative(obj, x, engine, surface)
            % 6-DOF rigid-body rate at state x for fixed actuator positions.
            % Atmosphere -> air-relative aero-propulsive forces/moments -> EOM
            % (gravity included inside RBD). Used as the RK4 stage function.
            [rho, a] = obj.environment.atmosphere(-x(3));
            R_i2b  = RSLQR.rotm_i2b(x(7), x(8), x(9));
            v_air  = obj.environment.airspeed_body(x(4:6), R_i2b);
            x_aero = [v_air; x(10:12)];
            [Fb, Mb] = run_LpC_aero(x_aero, engine, surface, rho, a, ...
                                    obj.units, obj.vehicleConfig);
            dx = obj.rigidBody.calculate_dynamics(x, Fb, Mb);
        end

        function Xlon_dot = lon_plant_rate(obj, x_full, U0, u_perturb)
            % True longitudinal state rate for the liveness filter's nonlinear
            % dV/dt. Assembles absolute effectors from trim (U0) plus the
            % 11-channel effector perturbation u_perturb = [Pi1..8; Pi9; de; df]
            % (same split/mapping as reset() and RSLQR.srface_cmd), then runs
            % the exact plant path of step() (atmosphere -> aero -> RBD).
            % Returns Xlon_dot = d/dt [u; w; q; theta] (rows [4;6;11;8] of dx).
            engine        = U0(5:13);
            engine(1:9)   = engine(1:9) + u_perturb(1:9);
            surface       = [U0(1) - U0(2); U0(1) + U0(2); U0(3); U0(3); U0(4)];
            surface(3)    = surface(3) + u_perturb(10);   % elevator -> LE
            surface(4)    = surface(4) + u_perturb(10);   % elevator -> RE
            surface(1)    = surface(1) + u_perturb(11);   % flap -> LA
            surface(2)    = surface(2) + u_perturb(11);   % flap -> RA

            [rho, a] = obj.environment.atmosphere(-x_full(3));
            R_i2b = RSLQR.rotm_i2b(x_full(7), x_full(8), x_full(9));
            v_air = obj.environment.airspeed_body(x_full(4:6), R_i2b);
            x_aero = [v_air; x_full(10:12)];
            [Fb, Mb] = run_LpC_aero(x_aero, engine, surface, rho, a, ...
                                    obj.units, obj.vehicleConfig);
            dx = obj.rigidBody.calculate_dynamics(x_full, Fb, Mb);
            Xlon_dot = dx([4; 6; 11; 8]);
        end
    end
end
