# LpC-GUAM m-code Refactoring Wiki

Documentation for the `Refactoring/` folder, a pure object-oriented MATLAB (m-code) port
of the Simulink-based Lift+Cruise (LpC) transition simulation from NASA LaRC's **GUAM**
(Generic UAM Simulation).

While the original GUAM supports many reference frames including ECI and ECEF, this
refactor assumes a **flat earth** and uses **only the NED / Body / Aero (α, β) frames**.
The rotating oblate earth, earth rotation, and geodetic navigation are not modeled.

- Original simulation: `GUAM.slx` + `simSetup.m` (repository root), [nasa/Generic-Urban-Air-Mobility-GUAM](https://github.com/nasa/Generic-Urban-Air-Mobility-GUAM)
- Reproduced scenario: `Exec_Scripts/exam_TS_Hover2Cruise_traj.m` (hover → transition → forward flight)
- Entry point: `Refactoring/run_transition_sim.m`

---

## Table of Contents

1. [Modification History](#1-modification-history)
2. [Reference Frames](#2-reference-frames)
3. [6-DOF Rigid Body Dynamics](#3-6-dof-rigid-body-dynamics)
4. [Aero-Propulsive Model (Polynomial)](#4-aero-propulsive-model)
5. [Control Surface and Rotor Actuator Models](#5-control-surface-and-rotor-actuator-models)
6. [Environment Model (Atmosphere/Wind)](#6-environment-model)
7. [Controller (RSLQR)](#7-controller-rslqr)
8. [Simulation Scenario and How to Run](#8-simulation-scenario-and-how-to-run)
9. [References](#9-references)

---

## 1. Modification History

### Commit history

| Date | Commit | Description |
|---|---|---|
| 2026-06-03 | `f7fb782` | Created LpC GUAM skeleton code (`LpC_GUAM`, `RBD`, config classes) |
| 2026-06-20 | `0e514c3` | Added `Environment` and `Gravity` classes |
| 2026-06-21 | `8e75d59` | Intermediate backup of refactoring work |
| 2026-06-30 | `4005894` | Copied polynomial aero model (`AeroPolynomial/`), ported RSLQR controller and trim table |
| 2026-06-30 | `8fff1ac` | Aligned `VehicleConfig`/`Units` with the source simulation conventions (slug/ft) |

### 2026-07-06 — Completed transition simulation (flat-earth alignment)

Completed the full closed-loop simulation matching the original GUAM default variant
configuration (Polynomial aero, FirstOrder actuators, Simple EOM), restricted to a
flat-earth formulation without ECI/ECEF.

| File | Change |
|---|---|
| `GUAM/LpC_GUAM.m` | Completed the orchestrator. Implemented `reset()` (initialize at the trim of the first reference point), `step()` (control → servos → atmosphere → aero → EOM → Euler integration), and `run()` (execute the 40 s transition scenario with logging) |
| `GUAM/EngineDynamics.m` | New — first-order rate/position limited servos for the 9 rotors (parameters from `setupEngines.m`) |
| `GUAM/SurfaceDynamics.m` | Rewritten — first-order rate/position limited servos for the 5 control surfaces (parameters from `setupActuators.m`). Removed the previous `rad2deg`/`deg2rad` swap bug |
| `GUAM/Environment.m` | Rewritten — removed the incomplete Dryden turbulence stub; replaced with US Standard Atmosphere (1976) troposphere + steady NED wind |
| `GUAM/RBD.m` | Bug fixes: nonexistent property `vehicleConfig.inertia` → `.I`, typo in element (2,3) of the body→NED rotation matrix (`sPhi*cPhi` → `sPhi*cPsi`), unit comments lb → slug |
| `controller/RSLQR.m` | Bug fixes: inconsistent call signatures inside `control()` (struct vs. separate arguments), dimension error in the weighted pseudo-inverse formula (`W\eye(3)` → `M = W⁻¹Bᵀ(BW⁻¹Bᵀ)⁻¹`), engine limit check `1:8` → `1:9` (matching the original `Limited_Long_Cont_Alloc.m`). Added `reset()`; corrected the trim input `U0` ordering comment |
| `run_transition_sim.m` | New — run script (path setup, simulation, result plots) |
| `config/SimConfig.m` | Simulation duration 10 s → 40 s (length of the transition scenario) |
| `config/RSLQRConfig.m` | Added Hover→Cruise reference trajectory tables (`ref_pos/ref_vel/ref_chi/ref_chidot`) |

**Verification**: Confirmed a full 40 s run via headless MATLAB R2025a execution.
Final u = 15.2 ft/s (ref 15), θ ≈ 1°, lift rotors stable at about 900 RPM. The position
error (max 40 ft on the North axis) is a lag characteristic caused by the inconsistency
built into the original example's reference itself (piecewise-linear position vs. ramped
velocity), not a divergence.

---

## 2. Reference Frames

The original GUAM provides bus support for ECI/ECEF/NED/Navigation/Velocity/Wind/
Stability/Body frames; this port uses only the four below.

### 2.1 NED inertial frame (flat earth)

- Origin: simulation start point; axes: North–East–Down.
- **Flat-earth assumption**: the NED frame is treated as inertial, ignoring earth
  curvature and rotation (a standard assumption for the low-speed, short-range eVTOL
  transition regime; see Stevens et al. [6], Sec. 1.7).
- Altitude is `h = -r_d`.

### 2.2 Body frame

- Origin: aircraft center of gravity (CG); axes: x-forward, y-right, z-down.
- The NED → Body transformation uses the 3-2-1 (yaw-pitch-roll) Euler angle sequence:

```
R_ib(φ,θ,ψ) = R_x(φ) · R_y(θ) · R_z(ψ)
```

Implementation: `RSLQR.rotm_i2b()` (static method), and `R_b2i = R_ib'` inside
`RBD.calculate_dynamics()`.

### 2.3 Heading (control) frame

The NED frame rotated about its z-axis by the course angle χ. The controller's gain
scheduling variables and velocity commands are defined in this frame (the
"heading-rotated NED frame" of the original GUAM baseline controller, Acheson et al. [1]).

- `UH` = heading-frame x-velocity, `WH` = z-velocity → the two scheduling axes of the
  trim/gain tables.

### 2.4 Aero frame (α, β)

`AeroFrame.m` — identical to the polynomial aero model convention (`LpC_aero_p_v2.m`):

```
V = ‖[u v w]‖,   α = atan(w/u),   β = asin(v/V)    (u > 0.01 ft/s)
α = β = 0                                          (u ≤ 0.01 ft/s, hover)
```

### 2.5 Units

Same English unit system as the original GUAM: length in ft, mass in slug, force in lbf,
angles in rad. `utils/Units.m` wraps the original `lib/utilities/setUnits.m` as a class.

---

## 3. 6-DOF Rigid Body Dynamics

Implementation: `GUAM/RBD.m` (+ `GUAM/Gravity.m`). This reproduces the original GUAM
default EOM variant, the **Simple EOM** (`VehicleEOMRef_Simp.slx`, `EOMEnum.Simple`), in
a flat-earth Euler-angle formulation. The formulation follows the standard flat-earth
6-DOF equations of Stevens, Lewis & Johnson [6].

### 3.1 State vector (12 states)

```
x = [ r_n  r_e  r_d   u  v  w   φ  θ  ψ   p  q  r ]ᵀ
      NED pos [ft]   body vel [ft/s]  Euler [rad]  body rates [rad/s]
```

### 3.2 Equations of motion

**Translational kinematics** (NED position derivative):

```
ṗ = R_bi(φ,θ,ψ) · [u v w]ᵀ        (R_bi = R_ibᵀ, body → NED)
```

**Rotational kinematics** (Euler angle derivatives):

```
[φ̇]   [1   sinφ·tanθ   cosφ·tanθ] [p]
[θ̇] = [0   cosφ        -sinφ    ] [q]
[ψ̇]   [0   sinφ/cosθ   cosφ/cosθ] [r]
```

**Translational dynamics** (body frame):

```
v̇_b = -ω × v_b + g_b(φ,θ) + F_b/m
```

where the body components of gravity (`Gravity.m`) are

```
g_b = g·[-sinθ,  sinφ·cosθ,  cosφ·cosθ]ᵀ,   g = 32.174 ft/s²
```

The original Simple EOM also adds gravity inside the EOM (`SimIn.Units.g0`); the aero
model output `F_b` does not include gravity.

**Rotational dynamics**:

```
ω̇ = I⁻¹ · (M_b - ω × I·ω)
```

### 3.3 Mass properties (`config/VehicleConfig.m`)

Based on the SACD Lift+Cruise reference configuration [4][5] (identical values to
`vehicles/Lift+Cruise/setup/LpC_model_parameters.m`):

| Item | Value |
|---|---|
| Mass m | 181.789 slug |
| Inertia I | diag(13051.7, 16660.8, 24735.1) slug·ft² |
| Wing area S / chord c̄ / span b | 186 ft² / 3.18 ft / 47.5 ft |
| Rotors | 8 lift rotors (D = 10 ft) + 1 pusher (D = 9 ft) |

### 3.4 Numerical integration

Fixed-step **forward Euler**, Δt = 0.01 s (`SimConfig.dt`). The controller, servos, and
EOM all update at the same rate (the original GUAM uses a fixed-step Simulink solver).

---

## 4. Aero-Propulsive Model

Implementation: `GUAM/AeroPolynomial/` — a direct copy of the **polynomial
(response-surface) model** from the original
`vehicles/Lift+Cruise/AeroProp/Polynomial/`.

### 4.1 Model overview (Simmons et al. [2])

A **full-envelope aero-propulsive model** identified at NASA LaRC by polynomial
regression on OVERFLOW CFD computational-experiment (DOE) data. Local models are built
for the hover, transition, and cruise regimes and blended in overlapping regions.

- Entry point: `run_LpC_aero.m` → `LpC_aero_p_v2.m` → `LpC_interp_p_v2.m`
- Components:
  - **Static aerodynamics** — regime-specific polynomial models
    (`LpC_Polynomial_Models_v2p1/MOF_k0p2/`) combined with sigmoid blending
    (blending_method = 2) according to flight condition
  - **Airframe damping** — `LpC_glider_damping.m` (unpowered airframe damping with
    respect to p, q, r)
  - **Isolated rotor/propeller models** — `RotorModels/` (Sec. 4.3)
- When the flight condition leaves the model validity envelope, `check_fac_limits.m`
  holds the boundary value (the original errors out by design; the ported copy is
  relaxed to hold the boundary).

### 4.2 Input/output conventions

```matlab
[Fb, Mb, Validity] = run_LpC_aero(x, n, d, rho, a, Units, Model)
```

| Argument | Meaning |
|---|---|
| `x = [u v w p q r]` | **Air-relative** body velocity [ft/s] and angular rates [rad/s] (after subtracting steady wind) |
| `n` (9×1) | Rotor speeds [rad/s] — simulation ordering (see mapping below) |
| `d` (5×1) | Control surfaces [rad]: **[left aileron, right aileron, left elevator, right elevator, rudder]** |
| `rho`, `a` | Air density [slug/ft³], speed of sound [ft/s] |
| `Fb`, `Mb` | Body-axis aero + propulsive forces [lbf] and moments [ft·lbf] (about the CG, **gravity excluded**) |

The flap is not a separate surface: it is added symmetrically to the left/right
ailerons — a **flaperon** convention (controller `srface_cmd`:
`[flap−ail, flap+ail, ele, ele, rud]`).

**Engine order mapping**: the rotor numbering N1–N9 internal to the polynomial model and
the simulation input vector `n` are related by
`n(1)→N1, n(2)→N3, n(3)→N5, n(4)→N7, n(5)→N2, n(6)→N4, n(7)→N6, n(8)→N8, n(9)→N9`
(`LpC_aero_p_v2.m`, lines 81–89). The trim table and control allocation follow the same
ordering.

### 4.3 Rotor models (`RotorModels/`)

The propeller inflow angle is computed first (`LpC_prop_inflow_angle.m`); then isolated
rotor models identified from wind tunnel/CFD data give each rotor's forces and moments
(`LpC_CCW_RotorModel_J040914.m` — lift rotors, `LpC_PusherProp_Model.m` — pusher). These
are transformed to body-axis forces/moments about the CG using the hub positions and
orientations (`Model.Prop_location`, `Prop_R_BH`) and summed (`LpC_Iso_PropFM.m`,
`LpC_Total_PropFM.m`). Reaction torque signs follow the rotor spin directions
(`prop_spin`).

---

## 5. Control Surface and Rotor Actuator Models

Implementation: `GUAM/SurfaceDynamics.m`, `GUAM/EngineDynamics.m`.
These match the original GUAM default variant `ActuatorEnum.FirstOrder` — the same
structure as the Simulink `ServosLib` block **"First Order Dynamic Rate and Position
Limited Servo"** (ported after inspecting the block structure inside the `.slx` files).

### 5.1 Servo equation

For each channel:

```
δ̇ = sat( ωₙ · (δ_cmd − δ),  ±RL )      (first-order lag + rate limit)
δ  ∈ [δ_min, δ_max]                     (position limit, clamped after integration)
```

Discretization: forward Euler, Δt = 0.01 s.

### 5.2 Parameters

| | Surfaces (5 ch: LA RA LE RE RUD) | Rotors (9 ch: 8 lift + 1 pusher) |
|---|---|---|
| Bandwidth ωₙ | 20 rad/s | 4π ≈ 12.57 rad/s |
| Rate limit RL | 100 rad/s | 100 rad/s² |
| Position limits | ±30° | lift [0, 1600 RPM], pusher [0, 2000 RPM] |
| Source | `setup/setupActuators.m` | `setup/setupEngines.m` (Polynomial fmType limits) |

Initial values are set to the trim values (`U0`) at simulation start (`reset()`).

---

## 6. Environment Model

Implementation: `GUAM/Environment.m`. A flat-earth reduction of the original
`Environment.slx` (atmosphere + wind/turbulence); ECI/ECEF-related elements (e.g.,
latitude-dependent gravity) are excluded.

### 6.1 Standard atmosphere (US Standard Atmosphere 1976 [7], troposphere)

The English-unit constants match the sea-level values of the original
`setup/setupAtmosphere.m` (ρ₀ = 0.0023769 slug/ft³, a₀ = 1116.45 ft/s).

```
T(h) = T₀ − L·h                    T₀ = 518.67 °R,  L = 0.00356616 °R/ft
P(h) = P₀ · (T/T₀)^(g₀/(L·R))      P₀ = 2116.22 lbf/ft², R = 1716.49 ft·lbf/(slug·°R)
ρ    = P/(R·T),   a = √(γ·R·T)     γ = 1.4
```

### 6.2 Wind

A steady NED wind `wind_ned` (default 0) is supported; the aero model receives the
air-relative body velocity

```
v_air_b = v_b − R_ib · v_wind_NED
```

The original Dryden turbulence model is out of scope for this port (still-air
assumption).

---

## 7. Controller (RSLQR)

Implementation: `controller/RSLQR.m`, `config/RSLQRConfig.m`, and the trim/gain table
`controller/trim_table_Poly_ConcatVer4p0.mat`.

This is a port of the original GUAM baseline **unified gain-scheduled LQRi controller**
(`vehicles/Lift+Cruise/Control/`, `CtrlEnum.BASELINE`). The design method is the
**Robust Servomechanism LQR (RSLQR) with generalized control allocation** of Acheson,
Gregory & Cook [1]. The theoretical foundation of the RSLQR follows Davison's robust
servomechanism theory [8] and the formulation in Lavretsky & Wise [9].

### 7.1 Architecture overview

```
Reference trajectory (pos, vel_H, χ, χ̇)
   │
   ▼
Guidance error (guidance_error)          e_pos_b = R_ib·(p_ref − p),  e_χ
   │
   ▼
Perturbation states/commands             deviations from trim (X0, U0), outer-loop gain 0.1
(perturb_lin_ctrl)                       Xlon = [ũ w̃ q̃ θ̃], Xlat = [ṽ p̃ r̃ φ̃]
   │
   ▼
RSLQR lon/lat control (lon_ctrl/lat_ctrl)   m_des = G·r + K_i·∫ − K_x·x
   │
   ▼
Weighted pseudo-inverse allocation       u = M·m_des,  M = W⁻¹Bᵀ(BW⁻¹Bᵀ)⁻¹ + limit scaling
(pseudo_alloc)
   │
   ▼
Effector commands: 9 rotors + 5 surfaces (trim U0 added back)
```

### 7.2 Gain scheduling

- Scheduling variables: heading-frame velocity commands (U_H, W_H).
- Grid: 28 U_H points (0 – 219.4 ft/s) × 3 W_H points (−7.5, 0, 11.67 ft/s).
- At each grid point, the trim state/input `XU0` (25×1) and the lon/lat gain and model
  matrices (K_i, K_x, K_v, F, G, C, C_v, W, B, A_p, B_p) are designed offline and
  stored; at runtime they are obtained by **bilinear interpolation** (`interp_mtrx`,
  `interp_xu0`).
- Trim input ordering: `U0 = [flap; ail; ele; rud; eng1…eng9]`.
- The linear models for the offline design come from numerical differentiation of the
  polynomial aero model (original `ctrl_scheduler_GUAM.m` /
  `get_lin_dynamics_heading.m`; ported copy
  `AeroPolynomial/poly_aero_wrapper_Mod.m` — not used at runtime).

### 7.3 Longitudinal/lateral axis separation (Acheson et al. [1])

| | Longitudinal (lon) | Lateral (lat) |
|---|---|---|
| Perturbation states x | [ũ, w̃, q̃, θ̃] | [ṽ, p̃, r̃, φ̃] |
| Tracking commands r | [u_cmd, w_cmd, 0] | [v_cmd, χ̇_cmd] |
| Effectors incl. virtual u | 9 rotors + elevator + flap + **θ_v** (12) | 8 lift rotors + aileron + rudder + **φ_v** (11) |

**Servo compensator (integrator)**: for each axis,

```
ẋ_i = F·r − C·x + K_v·(η_v − C_v·x)
m_des = G·r + K_i·x_i − K_x·x
```

where η_v is the virtual attitude effector (θ_v or φ_v) fed back from the allocation.
At low speed (hover) acceleration is produced by tilting the attitude, and at high speed
by the aerodynamic effectors — enabling **unified control across the entire transition
envelope**. This is the core of the "generalized control allocation" in [1].

**Allocation limit scaling**: when the rotors/elevator/flap exceed their position
limits, the whole allocation is scaled down while preserving the direction of m_des
(same logic as the original `Limited_Long_Cont_Alloc.m`).

**LQR weights** (`RSLQRConfig`): longitudinal Q = diag(0.01, 0.01, 1000, 0…) —
weighting concentrated on the integral error (especially the q̇ channel); the
allocation weight W penalizes the flap (10⁷) and elevator (10³) so that the allocation
favors rotors/attitude at hover and low speed.

### 7.4 Guidance (outer loop)

- The position error is rotated into the body frame and added to the velocity command
  with a gain of 0.1: `v_cmd += 0.1 · e_pos_b`
- Course error: `χ̇_cmd += 0.1 · e_χ`, with `e_χ` the shortest angle via `atan2`.

### 7.5 Discretization

The controller integrators, servos, and EOM all use Δt = 0.01 s forward Euler
(`RSLQRConfig.dt = SimConfig.dt`).

### 7.6 Longitudinal Liveness Filter (HJ Reachability)

A liveness filter (`LivenessFilter.m` + `ValueFunctionLUT.m`) runs *behind* the
nominal RSLQR allocation, projecting the longitudinal effector perturbation
`act_lon(1:11)` onto a set of liveness-ensuring controls derived from a
Backward Reachable Tube (BRT) value function `V(x)`. Theory and the full
distillation of paper §IV are in
[`docs/refs/liveness-filter.md`](../../docs/refs/liveness-filter.md); the
design decision log is **ADR-0002**.

- **State/control match**: the BRT state `x=[u,w,q,θ]` (perturbation from trim)
  equals RSLQR `Xlon`; the BRT control `[Π1..8,Π9,δe,δf]` equals `act_lon(1:11)`
  (the virtual attitude effector `act_lon(12)=θ_v` is not filtered).
- **Control-affine specialization** (converged `V ⇒ DtV=0`, `d*=0`):
  `dV/dt = α + β'u`, `α = ∇V·(A_lon x)`, `β = (∇V·B_lon)'`, with
  `A_lon,B_lon = ctrl_lon.Ap,Bp`.
- **Filters**: Least-Restrictive (eq.15, banded boundary `V≥−eps_band`) and
  Smooth Blending (eq.20, CBF-like QP `β'u ≤ −γV − α`, `γ=5.0`). The QP is a
  toolbox-free box+halfspace dual bisection.
- **UNITS**: everything is ft/s, rad/s, rad. The `*_scale` (m, deg) factors in
  helperOC `visualize_*_tube.m` are **visualization-only** and are NOT applied
  to the value-function grid.
- **Scheduling / coverage**: `ValueFunctionLUT` auto-parses `(UH,WH)` indices
  from the BRT files present in `tables/BRT/`, interpolates across UH, and
  returns `ok=false` outside coverage (⇒ filter passes through).
- **WH3-only limitation**: the current LON BRT covers only WH3. Note `UH/WH` are
  **body-frame** velocities `u,w` (verified: `WH == XU0(3)=body w`, NOT
  heading/inertial), so `WH3=+11.67` is a vertical descent at hover but ~level
  flight at cruise (trim θ≈5° projects body-w into forward motion). althold
  cruise (WH1–2, w≈−8…0) is outside coverage, so the filter passes through in
  althold (bit-identical to filter-off — verified).
- **Closed-loop status (root cause found, 2026-07-09)**: the filter math and
  integration are verified (Steps 1–6; `tests/verify_liveness_*.m`), and on the
  *linear* model `ẋ=A_lon·x+B_lon·u` the blend filter keeps `V≤0` (works). But
  on the true *nonlinear* plant the filter does not help: the linear
  `dV/dt=α+β'u` disagrees in sign with the actual `dV/dt` under a pitch-rate
  upset — `∇V` is θ-dominated and `dθ/dt=q` is pure kinematics (relative degree
  2, no direct control), so V rises while the filter (deep inside Glive, large
  `−γV` budget) stays inactive, then escapes before control can arrest q. Root
  cause: **the linear trim BRT is not a valid liveness certificate for the
  nonlinear plant** (implementation itself is correct; BRT model = filter model,
  @GUAM_LON.m:40-41). Fixes (separate work): regenerate the BRT with nonlinear
  reachability, or make the filter use the true nonlinear `∇V·f(x,u)`. See
  ADR-0002 and `plans/guam-liveness-filter.md` Step 6.5.

---

## 8. Simulation Scenario and How to Run

### 8.1 Hover → Cruise transition scenario

The same reference trajectory as the original `Exec_Scripts/exam_TS_Hover2Cruise_traj.m`,
embedded as tables in `RSLQRConfig` (Δt = 0.01 s, 4001 points):

| Segment | Maneuver | Reference |
|---|---|---|
| 0 – 20 s | Vertical climb 0 → 80 ft | w_H: −8 → 0 ft/s, position (0,0,0) → (0,0,−80) |
| 20 – 40 s | Forward acceleration + climb | u_H: 0 → 15 ft/s, position → (150,0,−100), χ = 0 |

The initial state is set to the trim of the first reference point (interpolated from the
trim table; W_H clamped to −7.5 ft/s) — corresponding to the initial-condition setup of
the original `setupTrim`.

### 8.2 Running

```matlab
>> run('Refactoring/run_transition_sim.m')
```

Adds paths → constructs `LpC_GUAM` → `run()` → produces four figures
(position/velocity/attitude/effectors). Programmatic use:

```matlab
addpath(genpath('Refactoring'));
guam = LpC_GUAM();
out  = guam.run();   % out.time, out.state(12×N), out.engine(9×N), out.surface(5×N), ...
```

### 8.3 Step pipeline (`LpC_GUAM.step`)

```
1. RSLQR.control(state, ref)        → rotor/surface commands
2. Engine/SurfaceDynamics.step()     → first-order servos (rate/pos limits)
3. Environment.atmosphere(-r_d)      → ρ, a
4. run_LpC_aero(v_air, n, d, ρ, a)   → F_b, M_b  (aero + propulsion)
5. RBD.calculate_dynamics()          → ẋ (gravity included)
6. x ← x + Δt·ẋ                      (forward Euler)
```

---

## 9. References

[1] Acheson, M. J., Gregory, I. M., and Cook, J., "Examination of Unified Control
    Approaches Incorporating Generalized Control Allocation," *AIAA SciTech Forum*,
    AIAA 2021-0999, 2021. — RSLQR unified controller and generalized control allocation
    (source paper for the controller)

[2] Simmons, B. M., Buning, P. G., and Murphy, P. C., "Full-Envelope Aero-Propulsive
    Model Identification for Lift+Cruise Aircraft Using Computational Experiments,"
    *AIAA AVIATION Forum*, AIAA 2021-3170, 2021. https://doi.org/10.2514/6.2021-3170
    — polynomial aero-propulsive model (the aero model used in this port)

[3] Cook, J. W., and Hauser, J., "A Strip Theory Approach to Dynamic Modeling of
    eVTOL Aircraft," *AIAA SciTech Forum*, AIAA 2021-1720, 2021.
    https://doi.org/10.2514/6.2021-1720
    — GUAM's alternative aero model (SFunction strip theory); not included in this port

[4] Silva, C., Johnson, W., Antcliff, K. R., and Patterson, M. D., "VTOL Urban Air
    Mobility Concept Vehicles for Technology Development," *AIAA AVIATION Forum*,
    AIAA 2018-3847, 2018. — NASA Lift+Cruise concept vehicle (origin of the airframe
    configuration)

[5] NASA SACD, "UAM Reference Vehicles," https://sacd.larc.nasa.gov/uam-refs/
    — Lift+Cruise reference configuration data (mass properties, rotor layout; source
    of the `VehicleConfig` values)

[6] Stevens, B. L., Lewis, F. L., and Johnson, E. N., *Aircraft Control and
    Simulation: Dynamics, Controls Design, and Autonomous Systems*, 3rd ed.,
    Wiley, 2015. — standard formulation of the flat-earth 6-DOF rigid body equations
    of motion

[7] NOAA/NASA/USAF, *U.S. Standard Atmosphere, 1976*, NOAA-S/T 76-1562, 1976.
    — standard atmosphere model (`Environment.m`)

[8] Davison, E. J., "The Robust Control of a Servomechanism Problem for Linear
    Time-Invariant Multivariable Systems," *IEEE Transactions on Automatic Control*,
    Vol. 21, No. 1, 1976, pp. 25–34. — foundation of robust servomechanism theory

[9] Lavretsky, E., and Wise, K. A., *Robust and Adaptive Control with Aerospace
    Applications*, Springer, 2013. — RSLQR (servomechanism LQR) design procedure

[10] NASA, "Generic Urban Air Mobility (GUAM) Simulation," GitHub repository,
     https://github.com/nasa/Generic-Urban-Air-Mobility-GUAM — original simulation
