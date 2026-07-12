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
        vehicleConfig   % VehicleConfig (also passed to the aero model as Model)
        simConfig       % SimConfig (owns dt/T and the scenario)
        refTraj         % Reference trajectory table (from SimConfig.getReferenceTrajectory)
        units           % Units ('ft','slug') — aero model unit conversions

        rigidBody       % RBD 6-DOF equations of motion
        aeroFrame       % AeroFrame (alpha/beta/airspeed for logging)
        environment     % Environment (flat-earth atmosphere + steady wind)
        controller      % RSLQR gain-scheduled controller
        engineDynamics  % EngineDynamics (9 rotor speed servos)
        surfaceDynamics % SurfaceDynamics (5 surface servos [LA RA LE RE RUD])

        state           % Current state [rn re rd u v w phi theta psi p q r]' (12x1)
        time            % Current simulation time [s]
    end

    methods
        function obj = LpC_GUAM(scenario, target_vel)
            if nargin < 1 || isempty(scenario)
                scenario = 'althold';
            end
            obj.vehicleConfig   = VehicleConfig();
            obj.simConfig       = SimConfig(scenario);
            obj.refTraj         = obj.simConfig.getReferenceTrajectory(target_vel);
            obj.units           = Units('ft', 'slug');

            obj.rigidBody       = RBD(obj.vehicleConfig);
            obj.aeroFrame       = AeroFrame();
            obj.environment     = Environment();
            obj.controller      = RSLQR();
            obj.engineDynamics  = EngineDynamics(obj.simConfig.dt);
            obj.surfaceDynamics = SurfaceDynamics(obj.simConfig.dt);

            obj.reset();
        end

        function reset(obj)
            % Initialize the vehicle at the trim condition of the first
            % reference-trajectory point (as GUAM's setupTrim does for the
            % initial reference velocity).
            uh0 = obj.refTraj.vel(1, 1);
            wh0 = obj.refTraj.vel(3, 1);
            [X0, U0] = obj.controller.interp_xu0(uh0, wh0);

            obj.state         = zeros(12, 1);
            obj.state(1:3)    = obj.refTraj.pos(:, 1);  % NED position
            obj.state(4:6)    = X0(1:3);            % body velocity at trim
            obj.state(7:9)    = X0(10:12);          % Euler angles at trim
            obj.state(10:12)  = X0(4:6);            % body rates at trim

            % Actuators start at trim settings
            obj.engineDynamics.reset(U0(5:13));
            obj.surfaceDynamics.reset([U0(1) - U0(2);...
                                       U0(1) + U0(2);...
                                       U0(3);...
                                       U0(3);...
                                       U0(4)]);
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

        function out = run(obj)
            % Run the full hover-to-cruise transition defined by the
            % reference trajectory table and return logged results.
            rt  = obj.refTraj;
            N   = size(rt.pos, 2);
            dt  = obj.simConfig.dt;

            out.time    = (0:N-1) .* dt;
            out.state   = zeros(12, N);
            out.engine  = zeros(9, N);
            out.surface = zeros(5, N);
            out.alpha   = zeros(1, N);
            out.beta    = zeros(1, N);
            out.V       = zeros(1, N);
            out.ref_pos = rt.pos;
            out.ref_vel = rt.vel;

            obj.reset();
            for k = 1:N
                ref.pos     = rt.pos(:, k);
                ref.vel     = rt.vel(:, k);
                ref.chi     = rt.chi(k);
                ref.chi_dot = rt.chidot(k);

                out.state(:, k) = obj.state;
                [out.alpha(k), out.beta(k), out.V(k)] = ...
                    obj.aeroFrame.compute(obj.state(4), obj.state(5), obj.state(6));

                [engine, surface] = obj.step(ref);
                out.engine(:, k)  = engine;
                out.surface(:, k) = surface;
            end
        end
    end
end
