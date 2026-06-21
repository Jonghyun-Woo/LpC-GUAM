# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

GUAM (Generic UAM Simulation), NASA LaRC's MATLAB/Simulink simulation of a generic Lift+Cruise eVTOL transition
vehicle. It is a Simulink model (`GUAM.slx`) driven by MATLAB setup scripts, with a 6-DOF rigid body, swappable
aero-propulsive/actuator/sensor/control variant subsystems, and reference-trajectory generation (ramps, timeseries,
piecewise Bezier curves, doublets). No build system, package manager, or test framework is present — this is MATLAB
script/model code intended to be run inside MATLAB/Simulink.

Current branch (`refactoring/mcode`) is mid-effort porting parts of the Simulink/m-code simulation into a pure
object-oriented MATLAB implementation under `Refactoring/` (see below) — most classes there are skeletons/stubs,
not yet wired into a runnable simulation.

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

## Architecture (in-progress `Refactoring/` port)

A parallel, non-Simulink OOP MATLAB implementation, intended to eventually replace the corresponding Simulink
blocks with plain class-based code:
- `Refactoring/config/VehicleConfig.m` — static vehicle parameters (mass, inertia, gravity, prop layout) in
  slug/ft units, mirroring `vehicles/Lift+Cruise` constants.
- `Refactoring/config/SimConfig.m` — simulation time step/duration.
- `Refactoring/GUAM/RigidBody6Dof.m` — standalone 6-DOF equations of motion (`calculate_dynamics`), state vector
  `[rn re rd u v w phi theta psi p q r]'`, taking body forces/moments and returning the state derivative.
- `Refactoring/GUAM/Gravity.m` — gravity-in-body-frame helper; currently stubbed to return zeros (real
  trig-based computation is commented out — do not assume it's wired up).
- `Refactoring/GUAM/Environment.m` — wind/turbulence (Dryden filter) model; constructor and most methods are
  partially stubbed/commented out, work in progress.
- `Refactoring/GUAM/LpC_GUAM.m` — top-level integration point for the refactor; currently empty.
- `Refactoring/utils/Units.m` — static ft/m, lb/kg unit conversion helpers.

When extending `Refactoring/`, check which methods are still stubs (commented-out bodies, hardcoded zero returns)
before assuming functionality exists — several classes here are skeletons, not working code.
