# LpC-GUAM — MATLAB m-code port of NASA GUAM (Lift+Cruise)

An object-oriented **MATLAB** reimplementation of NASA Langley's Simulink-based **GUAM** (Generic Urban Air Mobility Simulation) for the **Lift+Cruise (LpC)** eVTOL.

It reproduces GUAM's default closed-loop transition demo (`Exec_Scripts/exam_TS_Hover2Cruise_traj.m`: hover → transition → forward flight) as plain MATLAB classes — no Simulink required.

> This is a derivative fork, **not** the official NASA repository. For the original Simulink simulation see [nasa/Generic-Urban-Air-Mobility-GUAM](https://github.com/nasa/Generic-Urban-Air-Mobility-GUAM). This work is distributed under the same NASA Open Source Agreement (see [`license/`](license/) and the notice at the end of this file).

## Scope

To keep the port tractable it targets the low-speed, short-range eVTOL transition regime under a flat-earth assumption:

- **Frames:** NED Inertial / Body / Aero ($\alpha$, $\beta$)
- **Default variants only:** Polynomial aero-propulsive model, first-order rate/position limited actuators, 6-DOF rigid-body-dynamics.
- **Units:** English throughout:  ft, slug, lbf, rad.

## Architecture

Plant and controller are cleanly separated and assembled by a single config hub; an explicit
step loop drives them:

| Layer | Class | Role |
|---|---|---|
| Config hub | `config/Config.m` (+ `SimConfig`, `VehicleConfig`, `ControllerConfig`, `RSLQRConfig`, `FilterConfig`, `LoggerConfig`) | Central assembly of all simulation parameters |
| Controller | `controller/Controller.m` → `RSLQR` + `LivenessFilter` (`ValueFunction`) | Gain-scheduled LQRi control + allocation, optional liveness filter; returns absolute effector commands |
| Plant | `GUAM/LpC_GUAM.m` → `EngineDynamics`, `SurfaceDynamics`, `Environment`, `AeroPolynomial/`, `RBD`, `Gravity` | First-order servos → atmosphere → polynomial aero-propulsion → flat-earth 6-DOF EOM, RK4 integration |
| Trajectory | `trajectory/ReferenceTrajectory.m` | Builds the reference tables per scenario |
| Logging | `logger/SimLogger.m` (+ plotting helpers) | Buffers per-step data, plots, saves figures/mat |

Per step: `Controller.control(state, ref)` produces absolute rotor/surface commands, then `LpC_GUAM.step(...)` integrates the plant one `dt = 0.01 s` step under them.


## Status

- ✅ Full 40 s hover→cruise transition closed loop reproduces the GUAM default configuration.
- 🔬 The longitudinal liveness filter is a work in progress: the linear trim BRT is not yet a valid liveness certificate for the true nonlinear plant. The BRT value-function tables (`tables/`) are not included in this repository, so the filter currently passes through.****
