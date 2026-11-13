# Rocket Flight Dynamics & Attitude Control — MATLAB

A personal engineering project exploring **3D rocket flight dynamics, aerodynamic forces, and active pitch–yaw–roll control** in MATLAB.

The main simulation follows a hobby-scale rocket through **powered ascent, coast, apogee, and unpowered ballistic descent**. Four independently actuated, all-moving tail fins share control requests across three axes.

> **Project status:** The current integrated release is **v11**. This is an educational, physics-based simulation with several assumed/illustrative aerodynamic and structural parameters — **not a flight-validated model or deployable flight controller**. Two-way aeroelastic coupling is left for future work.

## Quick start

**Main file:** [`src/rocket_model_v11.m`](src/rocket_model_v11.m)

**Requirements**
- MATLAB (developed with MATLAB R2026a; compatibility with other releases is not independently verified).
- **Control System Toolbox** for transfer functions and discretization (`tf`, `c2d`).
- No Simulink model or external data files are required for the main script.

**Run**

1. Clone the repository:
   ```bash
   git clone https://github.com/sakes-eid/Rocket-dynamic-controller.git
   cd Rocket-dynamic-controller
   ```
2. Open `src/rocket_model_v11.m` in MATLAB.
3. Edit the parameters near the top of the script if desired.
4. Click **Run** (or call `run('src/rocket_model_v11.m')` from the repository root).
5. Inspect the generated figures and the command-window flight/control diagnostics.

The main script clears the workspace and closes existing figures when started.

## What's modeled

| Area | Features in the integrated v11 model |
| --- | --- |
| **3D flight** | World-frame X/Y/Z translation; body-frame forces; quaternion-based attitude integration |
| **Rotational dynamics** | Pitch, yaw and roll rates; distinct polar/transverse inertia; rigid-body gyroscopic coupling |
| **Attitude control** | Pitch/yaw feedback and a roll controller targeting the shortest angular correction toward 0° |
| **Control allocation** | Four independent fin servos with combined pitch/yaw/roll commands, first-order response and **±5°** travel limits |
| **Changing mass** | Propellant depletion, moving center of mass, and varying moments of inertia |
| **Aerodynamics** | Geometry-based drag terms, atmospheric density and temperature, Mach-dependent *illustrative* aerodynamic coefficient/CP schedules, individual fin forces, body crossflow forces |
| **Disturbances** | Configurable horizontal wind gusts in both world X and Y |
| **Structural sensing** | Simplified **4 Hz** bending/sensor-error oscillators, a bending estimator and a notch filter; **not** two-way physical aeroelastic deformation |
| **Descent** | Fin lock at neutral after apogee; passive fin/body aerodynamic loads; no parachute or touchdown-contact dynamics |
| **Diagnostics** | Flight trajectory, 3D ground track, attitude histories, control moments, individual fin motion, aerodynamic and descent plots |

The original 2D numerical stability study is kept in the historical versions; it has **not** been validated or carried over as a full three-axis stability result.

## Default demonstration scenario

The parameters in the current v11 script include:

| Parameter | Default |
| --- | --- |
| Rocket length / body diameter | 2.7 m / 0.12 m |
| Dry mass / initial fuel mass | 12 kg / 5 kg |
| Fuel flow rate / specific impulse | 0.5 kg/s / 180 s |
| Initial pitch / yaw / roll | +5° / +5° / 180° |
| Fin deflection limit | ±5° per fin |
| Servo time constant | 0.04 s |
| X / Y wind-gust peaks | −15 / +12 m/s |
| Gust interval | 7–11 s |
| Simulation timestep | 0.00025 s |

**Illustrative output from a v11 run** (not flight measurements):

| Measurement | Example result |
| --- | ---: |
| Apogee | ~5,226.5 m at 33.0 s |
| Maximum Mach | ~1.122 |
| Maximum dynamic pressure | ~66.8 kPa |
| Burnout pitch / yaw / roll | ~0.798° / 0.083° / −5.197° |
| Touchdown time | ~75.23 s |
| Touchdown speed | ~193.25 m/s |

These figures are **scenario-specific** and depend on the model's assumed parameters, force laws and numerical implementation. They should not be treated as predictions for a real rocket.

## Repository structure

```text
Rocket-dynamic-controller/
├── src/
│   └── rocket_model_v11.m            # Main integrated 3-axis simulation
├── versions/
│   ├── rocket_model_v1_ridgid.m
│   ├── rocket_model_v2_flexible.m
│   ├── rocket_model_v3_aero.m
│   ├── rocket_model_v4_2D.m
│   ├── rocket_model_v5_mach_drag.m
│   ├── rocket_model_v6.m
│   ├── rocket_model_v6_drag.m
│   ├── rocket_model_v7_numerical_stability.m
│   ├── rocket_model_v8_aeroelastic.m
│   ├── rocket_model_v8_ballistic_descent.m
│   ├── rocket_model_v9_descent_diagnostics.m
│   ├── rocket_model_v10_pitch_yaw_3d.m
│   ├── rocket_model_v12_beam_4hz_test.m
│   └── rocket_model_v13_modal_frequency_tests.m
└── results/
    └── MATLAB simulation screenshots
```

[Browse all historical versions](versions/) · [Browse simulation screenshots](results/)

### Development history

- **v1–v4:** Early rigid/flexible and aerodynamic experiments, progressing to 2D flight dynamics.
- **v5–v7:** Expanded drag/Mach studies and pitch-only numerical stability diagnostics.
- **v8–v9:** Experimental aeroelastic branch; separate ballistic-descent implementation and descent diagnostics.
- **v10:** Transition from 2D to 3D pitch/yaw dynamics and trajectory plots.
- **v11 (main):** Three-axis pitch/yaw/roll control, four-fin mixing and full flight/descent visualization.
- **v12–v13 (separate experiments):** Assumed `EI = 40,000 N·m²` beam-bending response, free decay, imposed 4 Hz excitation and frequency sweep. **These are diagnostic tests, not integrated aeroelastic feedback in v11.**

The historical versions document the development process. They may use different assumptions and are not all intended as independently validated releases.

## Parameters worth exploring

In the main script, try adjusting:

```matlab
theta_initial = deg2rad(5);
yaw_initial   = deg2rad(5);
roll_initial  = deg2rad(180);

wind_peak   = -15;  % world-X gust, m/s
wind_peak_y = 12;   % world-Y gust, m/s

roll_tie_direction = +1;  % +1 or -1 when exactly 180 degrees off
delta_max = deg2rad(5);   % individual fin travel limit
sensor_mode = 2;          % 0: ideal, 1: flexible, 2: compensated
```

For the 4 Hz versus physical-body-bending investigation, see [`rocket_model_v12_beam_4hz_test.m`](versions/rocket_model_v12_beam_4hz_test.m) and [`rocket_model_v13_modal_frequency_tests.m`](versions/rocket_model_v13_modal_frequency_tests.m).

## Assumptions and limitations

- This is a **simulation study**, not a verified flight-control system. **Do not use its outputs or controller directly to fly hardware.**
- Many aerodynamic coefficients and their Mach-dependent corrections are **illustrative assumptions**, not calibrated wind-tunnel or CFD data. Transonic and large-angle-of-attack behavior is particularly uncertain.
- The propulsion model uses fixed mass flow and specific impulse during burn rather than a measured thrust curve.
- Fin motion is represented by first-order actuator dynamics and deflection limits; a detailed servo torque/current or hinge-load model is **not** part of the integrated v11 controller.
- The 4 Hz structural component in v11 is a **simplified sensing disturbance/estimation model**. It is not the EI-based airframe bending solution, and it does not deform the aerodynamic geometry.
- Post-apogee aerodynamics are simplified; no parachute, ground contact, impact shock or recovery system is modeled.
- Full three-axis closed-loop eigenvalue/stability validation, experimental parameter identification and timestep/uncertainty studies remain open tasks.

## Future work

- Physically consistent **two-way aeroelastic coupling** (structural deflection ↔ local aerodynamic forces) and validated stiffness/mode shapes.
- Fin torsion, hinge moments and rigorous control-reversal/flutter analysis if supported by suitable aerodynamic data.
- Higher-fidelity atmosphere, propulsion data, unsteady aerodynamics and real actuator limits.
- Three-axis stability analysis, automated regression tests, calibration against measured data and broader scenario sweeps.

## Project purpose

This project grew out of studying control systems and became an exploration of how **dynamics, aerodynamics, actuation, numerical modeling and feedback control** interact. The source code and prior versions are shared for study and experimentation, with limitations stated explicitly.
