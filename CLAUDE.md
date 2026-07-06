# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

GUAM (Generic UAM Simulation), NASA LaRC's MATLAB/Simulink simulation of a generic Lift+Cruise eVTOL transition
vehicle. It is a Simulink model (`GUAM.slx`) driven by MATLAB setup scripts, with a 6-DOF rigid body, swappable
aero-propulsive/actuator/sensor/control variant subsystems, and reference-trajectory generation (ramps, timeseries,
piecewise Bezier curves, doublets). No build system, package manager, or test framework is present — this is MATLAB
script/model code intended to be run inside MATLAB/Simulink.

Current branch (`refactoring/mcode`) ports the Simulink simulation into a pure object-oriented MATLAB
implementation under `Refactoring/` (see below). As of 2026-07-06 this port is a **runnable end-to-end
closed-loop simulation** of the hover-to-cruise transition scenario, restricted to a flat-earth formulation
(NED/Body/Aero frames only — no ECI/ECEF, no turbulence). Entry point: `Refactoring/run_transition_sim.m`;
full documentation in `Refactoring/docs/WIKI.md`.

## Running the simulation

There is no CLI; everything runs from the MATLAB command window with the repo root as the working directory.

- `RUNME.m` — top-level interactive entry point. Prompts the user to pick a demo trajectory (1-5), then calls
  `sim(model)` and `simPlots_GUAM` to plot results.
- `setupPath.m` — must run first (directly or via the scripts above); calls `restoredefaultpath` and adds
  `ClassDef/`, `Environment/`, `lib/` (+ subfolders), `utilities/`, `setup/`, `Exec_Scripts/`, `vehicles/`,
  `Bez_Functions/`, `Challenge_Problems/` to the MATLAB path.
- `simSetup.m` — the core setup pipeline: builds `userStruct.variants` → `SimIn` (via `setupTypes`), applies
  switches (`setupSwitches`), resolves the desired reference trajectory into `target.RefInput`, then calls
  `setup(SimIn, target)`, `setupParameters`, `setupVariants`, `setupBuses`, and builds the user-defined SimOut bus.
- Demo trajectory scripts live in `Exec_Scripts/` (files ending `_traj.m`, plus `exam_RAMP.m`, `exam_Bezier.m`).
  Each sets `userStruct`/`target`, calls `simSetup`, then the user runs the open Simulink model.
- `Challenge_Problems/RUNME.m` is a separate demo of an own-ship trajectory + failure scenario for autonomy
  research (collision avoidance, effector failures, etc.); the `Generate_*.m` scripts there regenerate the
  bundled datasets.

To change a variant programmatically without a full setup script, e.g.:
```matlab
userStruct.variants.actType = 3; simSetup
```
Variant enumerations live in `ClassDef/*Enum.m` (e.g. `ActuatorEnum`, `ForceMomentEnum`, `RefInputEnum`). Default
variant choices are vehicle-specific, in `vehicles/Lift+Cruise/setup/setupDefaultChoices.m`.

There are no automated tests or linters in this repo — verification is done by running the Simulink model and
inspecting `logsout{1}.Values` (commonly assigned to `SimOut`) or using `vehicles/Lift+Cruise/Utils/simPlots_GUAM.m`
and `utilities/Animate_SimOut.m`.

## Architecture (current/Simulink simulation)

- **Three core structures** flow through the whole simulation: `SimIn` (fixed simulation inputs / variant
  selections), `SimPar` (tunable simulation parameters), `SimOut` (logged outputs, from `logsout{1}.Values`).
- **Variant system**: `userStruct.variants` selects among swappable subsystem implementations (vehicle type,
  reference-input type, atmosphere, turbulence, controller, actuator, propulsion, force/moment model, EOM,
  sensors, experiment type). Each axis has a `ClassDef/<X>Enum.m` and a `setup/select<X>Type.m` /
  `utilities/variant_definitions/setup<X>Variants.m` pair that wires the Simulink variant subsystem.
- **Two interchangeable aero-propulsive models** (`fmType`): `Polynomial` (default — fast, response-surface/CFD
  fit, found in `vehicles/Lift+Cruise/AeroProp/Polynomial/`, will error outside its valid flight envelope by
  design) and `SFunction` (slower, first-principles strip-theory class/object model under
  `vehicles/Lift+Cruise/AeroProp/SFunction/`, configurable rotor/surface layout).
- **Vehicle definition** lives under `vehicles/Lift+Cruise/` — `ClassDef/` (vehicle-specific classes), `Control/`
  (baseline LQRi gain-scheduled controller, unified across hover/transition/cruise, operating in the
  heading-rotated NED frame), `Subsystems/` (`.slx` blocks: AeroProp, GearSystem, PowerSystem, Sensors, SurfEng,
  VehicleEOM), `Trim/` (offline trim routines; top level `trim_helix.m`, cost `mycost.m`, constraints
  `nlinCon_helix.m`), `setup/` (vehicle-specific setup incl. `User_SimOut/` user output scripts).
- **Reference trajectories**: selected via `userStruct.variants.refInputType` (`FOUR_RAMP`, `ONE_RAMP`,
  `TIMESERIES`, `BEZIER`, `DEFAULT`/doublets). Bezier curve generation/evaluation utilities are in
  `Bez_Functions/` (`genPWCurve.m`, `evalPWCurve.m`, `evalSegments.m`, etc.); a `pwcurve.waypoints` /
  `pwcurve.time_wpts` `.mat` file or a `target.RefInput.Bezier` struct can supply waypoints, expressed in the NED
  frame.
- **Signal buses**: defined in `utilities/bus_definitions/*.m` and assembled by `setup/setupBuses.m` /
  `utilities/buildBusObject.m`; support Earth-Centered Inertial, ECEF, NED, Navigation, Velocity, Wind, Stability,
  and Body reference frames.
- **Quaternion/rotation math** utilities (`Qmult`, `QmultSeq`, `Qinvert`, `Qtrans`, `QrotX/Y/Z`, `qcheck4`) live in
  `lib/utilities/`.
- **Gain scheduling**: `vehicles/Lift+Cruise/Control/ctrl_scheduler_GUAM.m` is the top-level script; it schedules
  longitudinal and lateral axes separately via `get_lin_dynamics_heading.m` (linearizes at a flight condition) and
  axis-specific scripts (`get_lat_dynamics_heading.m`, `ctrl_lat.m`, etc.).

## Architecture (`Refactoring/` m-code port)

A parallel, non-Simulink OOP MATLAB implementation replacing the corresponding Simulink blocks with plain
class-based code. It mirrors GUAM's default variant configuration (Polynomial aero, FirstOrder actuators,
Simple EOM, baseline controller) under a flat-earth assumption. Status: **runnable** — verified headless in
MATLAB R2025a against the Hover2Cruise transition (40 s, stable, tracks the velocity reference).

- `Refactoring/run_transition_sim.m` — entry script: adds paths, runs the closed-loop sim, plots
  position/velocity/attitude/effectors vs. reference. Programmatic use: `guam = LpC_GUAM(); out = guam.run();`.
- `Refactoring/GUAM/LpC_GUAM.m` — top-level orchestrator. `reset()` initializes at the trim of the first
  reference point; `step(ref)` runs one closed-loop step (RSLQR control → engine/surface servos → atmosphere →
  polynomial aero → 6-DOF EOM → forward Euler, Δt = 0.01 s); `run()` executes the whole scenario with logging.
- `Refactoring/GUAM/RBD.m` — flat-earth 6-DOF Euler-angle EOM (`calculate_dynamics`), state
  `[rn re rd u v w phi theta psi p q r]'`; gravity added inside via `Gravity.m` (implemented, body-frame tilt).
- `Refactoring/GUAM/AeroPolynomial/` — copy of the polynomial aero-propulsive model from
  `vehicles/Lift+Cruise/AeroProp/Polynomial/` (entry `run_LpC_aero.m`). Surface input order is
  `[LA RA LE RE RUD]` (flaperon convention); engine vector `n` maps to internal rotors as
  `n(1..9) → N1 N3 N5 N7 N2 N4 N6 N8 N9`. `check_fac_limits.m` here holds boundary values instead of erroring.
- `Refactoring/GUAM/EngineDynamics.m` / `SurfaceDynamics.m` — first-order rate/position-limited servos matching
  ServosLib "First Order Dynamic Rate and Position Limited Servo" (engines: ωn = 4π rad/s, [0,1600/2000 RPM];
  surfaces: ωn = 20 rad/s, ±30°; both RL = 100).
- `Refactoring/GUAM/Environment.m` — flat-earth US Standard Atmosphere 1976 troposphere (`atmosphere(h)`) +
  steady NED wind (`airspeed_body`); no turbulence.
- `Refactoring/GUAM/AeroFrame.m` — α/β/airspeed per the polynomial-model convention (used for logging).
- `Refactoring/controller/RSLQR.m` — port of the baseline unified RSLQR controller (AIAA 2021-0999):
  bilinear gain scheduling over (UH, WH) from `trim_table_Poly_ConcatVer4p0.mat`, lon/lat servo-compensator
  LQR, weighted pseudo-inverse allocation `M = W⁻¹Bᵀ(BW⁻¹Bᵀ)⁻¹` with limit scaling, virtual attitude
  effectors θ_v/φ_v. Trim input ordering: `U0 = [flap; ail; ele; rud; eng1..9]`.
- `Refactoring/config/` — `VehicleConfig.m` (SACD L+C mass/geometry/prop layout, slug/ft), `SimConfig.m`
  (Δt = 0.01 s, T = 40 s), `RSLQRConfig.m` (LQR/allocation weights + embedded Hover2Cruise reference tables
  `ref_pos/ref_vel/ref_chi/ref_chidot`, matching `Exec_Scripts/exam_TS_Hover2Cruise_traj.m`).
- `Refactoring/utils/Units.m` — class wrapper of `lib/utilities/setUnits.m`; the aero model consumes `.deg`
  and `.knot`.
- `Refactoring/docs/WIKI.md` — full documentation: modification history, frames, dynamics/aero/actuator/
  environment models, controller theory, and paper references.

Known characteristics: position tracking shows bounded lag/overshoot (max ~40 ft N) because the original
example's reference position (piecewise-linear) and velocity (ramp) are mutually inconsistent and the outer
position loop gain is weak (0.1) — this matches the source simulation's behavior, not a bug.
`AeroPolynomial/poly_aero_wrapper_Mod.m` is offline linearization tooling, not used at runtime.
