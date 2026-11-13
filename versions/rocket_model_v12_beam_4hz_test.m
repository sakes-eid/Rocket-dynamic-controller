clear; clc; close all

%% V12 TEST: V11 3-axis model + one-way free-free beam and 4 Hz disturbance
% World [X,Y,Z(up)]; body [+X forward, +Y, +Z].
% Initial pitch +5 deg, yaw +5 deg, roll 180 deg.
% Four identical all-moving fins share pitch/yaw/roll requests.
% One-way EI-based body bending added; legacy 4 Hz sensor-side mode retained.
% Mach/aerodynamic and descent correlations remain provisional.

%% Parameters

% Integration
dt = 0.00025;
t = 0:dt:120;
n = length(t);

% Horizontal wind gust
wind_peak = -15;       % m/s, world X gust (pitch plane)
wind_peak_y = 12;      % m/s, world Y gust (yaw plane)
gust_start = 7;     % s
gust_duration = 4; % s

% Descent: unpowered, fins locked neutral, no parachute.
% The reverse-flow coefficient is provisional and can be calibrated later.
Cd_tail_first = 1.0;

% Rocket mass
m_dry = 12;
m_fuel0 = 5;
mdot = 0.5;

% Rocket geometry
L_body = 2.7;
D_body = 0.12;
L_nose = 0.45;

x_dry = 1.2;
x_fuel = 2.0;
I_dry = 8;
D_motor = 0.06;  % Motor diameter (m)

% Four all-moving fins
n_fins = 4;
n_pitch_fins = 2;

fin_span = 0.08;
fin_root = 0.22;
fin_tip = 0.10;
fin_thickness = 0.004;

x_fin_root_LE = L_body-fin_root;
fin_tip_LE_offset = 0.06;
fin_pivot_fraction = 0.25;
fin_load_fraction = 0.30;

% Propulsion
Isp = 180;
g = 9.81;

% Atmosphere and drag
rho0 = 1.225;
H = 8500;

% Illustrative Mach-Cd reference curve (plot only; not flight forces)
Cd0 = 0.5;
M_table = [0 0.5 0.7 0.8 0.9 1.0 1.1 1.2 1.4 1.6 2.0];
Cd_ratio = [1 1 1.03 1.10 1.25 1.55 1.70 1.60 1.35 1.25 1.20];

% Fin aerodynamics
CL_alpha = 3;
CL_max = 1;
delta_max = deg2rad(5);

% PD controller
Kp_fin = 0.5;
Kd_fin = 0.2;
theta_initial = deg2rad(5);  % Pitch initial offset
yaw_initial = deg2rad(5);    % Independent yaw initial offset
roll_initial = deg2rad(180); % 180 deg from desired body roll orientation
enable_yaw = true;
enable_roll = true;

% At precisely 180 degrees BOTH roll paths are equally short.
% +1 chooses positive right-hand body-X rotation, -1 the negative path.
roll_tie_direction = +1;

% Roll PD uses a target angular acceleration, converted to fin deflection
% using instantaneous q and the fin roll moment derivative. Tunable.
roll_wn = 1.5;            % rad/s, desired roll correction bandwidth
roll_zeta = 1.0;         % critical damping target
roll_q_enable = 30;      % Pa; fins cannot roll the rocket at zero airspeed
roll_max_contribution = deg2rad(3); % roll allocation before combined clipping
fin_mass_est = 0.14;     % kg/fin; provisional inertia estimate

% Servo
tau_s = 0.04;

% Structural flexibility
f_flex = 4;
zeta_flex = 0.05;
K_flex = 0.001;
wn_flex = 2*pi*f_flex;

% 4 Hz structural notch filter
f_notch = f_flex;
zeta_notch_zero = 0.01;
zeta_notch_pole = 0.15;
wn_notch = 2*pi*f_notch;

Cd_body_cross = 1.1;

e_fin = 0.8;  % Provisional fin efficiency

% Rail buttons / roughness (provisional)
n_rail_buttons = 2;
D_rail_button = 0.008;

roughness_selected = 5e-6;  % m

% Sensor mode: 0 ideal, 1 flexible, 2 compensated flexible
sensor_mode = 2;

% Switch OFF for a direct comparison with the pre-Mach-stability model.
enable_mach_aero = true;

% PRELIMINARY Mach correction tables (scenario assumptions, NOT CFD data).
% 0.25c -> 0.50c shifts the fin aerodynamic center rearward at high Mach.
% Interpolate smoothly across the transonic region: no Prandtl-Glauert
% singularity at Mach 1. Replace with validated CNa/Cm/CP data later.
Ma_nodes = [0 0.3 0.6 0.8 0.9 1.0 1.1 1.2 1.4 1.8 2.0];
nose_gain_nodes = [1 1 1.02 1.06 1.10 1.15 1.18 1.15 1.08 1.00 0.95];
fin_gain_nodes  = [1 1 1.02 1.09 1.16 1.25 1.32 1.30 1.20 1.12 1.08];
ctrl_gain_nodes = [1 1 1.02 1.07 1.12 1.18 1.22 1.22 1.16 1.10 1.05];
fin_ac_nodes    = [0.25 0.25 0.25 0.25 0.26 0.28 0.33 0.37 0.42 0.47 0.50];

if ~enable_mach_aero
    nose_gain_nodes(:) = 1;
    fin_gain_nodes(:) = 1;
    ctrl_gain_nodes(:) = 1;
    fin_ac_nodes(:) = 0.25;
end

pp_nose = pchip(Ma_nodes,nose_gain_nodes);
pp_fin  = pchip(Ma_nodes,fin_gain_nodes);
pp_ctrl = pchip(Ma_nodes,ctrl_gain_nodes);
pp_ac   = pchip(Ma_nodes,fin_ac_nodes);

% Bending estimator with deliberate parameter mismatch
K_flex_est = 0.8*K_flex;
f_flex_est = 1.05*f_flex;
zeta_flex_est = 1.2*zeta_flex;
wn_flex_est = 2*pi*f_flex_est;

%% Geometry and mass properties

% Body geometry
A_ref = pi*D_body^2/4;
R_body = D_body/2;

% Fin geometry
A_fin = 0.5*(fin_root+fin_tip)*fin_span;
S_fin = n_pitch_fins*A_fin;

c_mean = (2/3)*(fin_root^2+fin_root*fin_tip ...
    +fin_tip^2)/(fin_root+fin_tip);
hinge_arm = c_mean*(fin_load_fraction-fin_pivot_fraction);

% Fin aerodynamic center
x_MAC_LE = x_fin_root_LE ...
    +(fin_tip_LE_offset/3) ...
    *(fin_root+2*fin_tip)/(fin_root+fin_tip);

x_fin_CP = x_MAC_LE+0.25*c_mean;

% Mass and inertia
m_fuel = max(0,m_fuel0-mdot*t);
m_total = m_dry+m_fuel;

x_CoM = (m_dry*x_dry+m_fuel*x_fuel)./m_total;

I_total = I_dry ...
    +m_dry*(x_dry-x_CoM).^2 ...
    +m_fuel.*(x_fuel-x_CoM).^2;

% Provisional polar inertia of a cylindrical body plus four fins.
% Fin mass is already included in dry mass: do NOT add it twice.
% Recomputed below after fin_radial_cp is known.

% Geometry-based center of pressure
CN_nose = 2;
x_CP_nose = 2*L_nose/3;

mid_sweep = fin_tip_LE_offset ...
    +(fin_tip-fin_root)/2;

l_mid = hypot(fin_span,mid_sweep);

K_fb = 1+R_body/(R_body+fin_span);

CN_fins = K_fb ...
    *(4*n_fins*(fin_span/D_body)^2) ...
    /(1+sqrt(1+(2*l_mid/(fin_root+fin_tip))^2));

CN_total = CN_nose+CN_fins;

x_CP = (CN_nose*x_CP_nose ...
    +CN_fins*x_fin_CP)/CN_total;

% Flight-dependent CP, fin lever, and restoring stability.
l_fin = zeros(1,n);
static_margin = zeros(1,n);
x_CP_hist = zeros(1,n);
x_fin_CP_hist = zeros(1,n);
CN_nose_hist = zeros(1,n);
CN_fins_hist = zeros(1,n);
CL_alpha_hist = zeros(1,n);

%% Aerodynamic precomputation

drag_L = [L_body-L_nose L_nose c_mean];

drag_S = [ ...
    pi*D_body*(L_body-L_nose), ...
    pi*R_body*hypot(L_nose,R_body), ...
    2*n_fins*A_fin];

K_body = 1+D_body/(2*L_body);
K_fin = 1+2*fin_thickness/c_mean;
drag_K = [K_body K_body K_fin];

Cf_rough0 = 0.032*(roughness_selected./drag_L).^0.2;

% Conical nose coefficients
s_eps = R_body/hypot(L_nose,R_body);

Cd_n0 = 0.8*s_eps^2;
Cd_n1 = s_eps;
dCd_n1 = 4/2.4*(1-0.5*Cd_n1);

M_n2 = 1.3;
Cd_n2 = 2.1*s_eps^2 ...
    +0.5*s_eps/sqrt(M_n2^2-1);

dCd_n2 = -0.5*s_eps*M_n2/(M_n2^2-1)^1.5;
nose_exponent = dCd_n1/(Cd_n1-Cd_n0);

% Fin edges and external components
A_fin_edge = fin_span*fin_thickness;
sweep_LE = atan(fin_tip_LE_offset/fin_span);

A_base = pi*R_body^2;
A_motor = pi*D_motor^2/4;
A_button = pi*D_rail_button^2/4;

% Four-fin positions and orientation
fin_phi = (0:n_fins-1)'*2*pi/n_fins;

fin_span_dir = [zeros(n_fins,1), ...
    cos(fin_phi),sin(fin_phi)];

fin_normal = [zeros(n_fins,1), ...
    -sin(fin_phi),cos(fin_phi)];

% Radial location of fin aerodynamic center
fin_radial_cp = R_body + fin_span ...
    *(fin_root+2*fin_tip)/(3*(fin_root+fin_tip));

% Effective fin aspect ratio
fin_AR = 2*fin_span^2/A_fin;
I_roll = 0.5*(m_total-n_fins*fin_mass_est)*R_body^2 ...
    +n_fins*fin_mass_est*fin_radial_cp^2;
assert(all(I_roll>0),'Invalid polar inertia: check fin_mass_est');

% Four-fin actuator mixer: per-fin [pitch, yaw, roll] decomposition.
% It retains pitch/yaw pair conventions when the other inputs are zero.
fin_mix = [fin_normal(:,2), fin_normal(:,3), ones(n_fins,1)];
fin_mix_inv = pinv(fin_mix);

% Body side-view area
A_side_nose = 0.5*D_body*L_nose;
A_side_cyl = D_body*(L_body-L_nose);
A_side_body = A_side_nose+A_side_cyl;

% Crossflow center of pressure
x_side_nose = 2*L_nose/3;
x_side_cyl = (L_nose+L_body)/2;

x_CP_cross = (A_side_nose*x_side_nose ...
    +A_side_cyl*x_side_cyl)/A_side_body;

%% States and notch filter

q_flex_est = zeros(1,n);
qdot_flex_est = zeros(1,n);
torque_est_history = zeros(1,n);
q_flex_yaw_est = zeros(1,n);
qdot_flex_yaw_est = zeros(1,n);
torque_est_yaw_history = zeros(1,n);

drag_force = zeros(1,n);
wind_history = zeros(1,n);
wind_y_history = zeros(1,n);

% Inertial translational acceleration histories (world frame, m/s^2).
ax_history = zeros(1,n);
ay_history = zeros(1,n);
az_history = zeros(1,n);

% Flight
V = zeros(1,n);
h = zeros(1,n);

% Horizontal motion in both inertial directions
vx = zeros(1,n);
vy = zeros(1,n);
x_pos = zeros(1,n);
y_pos = zeros(1,n);

rho = zeros(1,n);
q_dyn = zeros(1,n);

% Pitch
theta = zeros(1,n);
yaw = zeros(1,n);
omega = zeros(1,n);         % pitch rate, body +Z
omega_yaw = zeros(1,n);     % yaw rate, body -Y
delta = zeros(1,n);         % pitch fin pair
delta_yaw = zeros(1,n);     % achieved yaw mixer component
roll = zeros(1,n);          % roll about body X relative to zero-roll frame
roll_spin = zeros(1,n);     % accumulated body-X rotation (unwrapped)
omega_roll = zeros(1,n);    % true body-X roll rate
omega_roll_dot = zeros(1,n);
delta_roll = zeros(1,n);    % achieved roll mixer component
fin_deflection = zeros(n_fins,n); % four independent servo angles [rad]
fin_command = zeros(n_fins,n);    % clipped per-fin setpoints [rad]
fin_mixer_scale = ones(1,n);
roll_error_history = zeros(1,n);
roll_requested_moment = zeros(1,n);
roll_control_sign = zeros(1,n);
quat_bw = zeros(n,4);       % body-to-world Hamilton quaternion [w x y z]
flow_incidence_history = zeros(1,n);

% Flexibility
q_flex = zeros(1,n);
qdot_flex = zeros(1,n);
q_flex_yaw = zeros(1,n);
qdot_flex_yaw = zeros(1,n);

% Aerodynamic moments
torque_fin = zeros(1,n);
torque_total = zeros(1,n);
torque_fin_yaw = zeros(1,n);
torque_total_yaw = zeros(1,n);
torque_roll = zeros(1,n);          % total body-X aerodynamic moment
active_torque_roll = zeros(1,n);   % fin deflection contribution

theta(1) = theta_initial;
yaw(1) = yaw_initial;
if ~enable_yaw, yaw(1) = 0; end
if ~enable_roll, roll_initial = 0; end

% The upright zero-roll reference aligns +X_body with the rocket nose,
% +Y_body approximately east and +Z_body approximately north.
nose_init = [tan(theta(1));tan(yaw(1));1];
nose_init = nose_init/norm(nose_init);
R_zero = rocketZeroRollFrame(nose_init);
R_spin = [1 0 0;0 cos(roll_initial) -sin(roll_initial); ...
          0 sin(roll_initial) cos(roll_initial)];
quat_bw(1,:) = rocketRotmToQuat(R_zero*R_spin);
roll(1) = roll_initial;
roll_spin(1) = roll_initial;
clear R_zero R_spin nose_init

% The v10 ascent was missing event initialization; carry it explicitly.
descent_active = false;
apogee_lock_idx = NaN;
impact_idx = NaN;
impact_time = NaN;

% Discrete structural notch filter
N_cont = tf( ...
    [1 2*zeta_notch_zero*wn_notch wn_notch^2], ...
    [1 2*zeta_notch_pole*wn_notch wn_notch^2]);

N_disc = c2d(N_cont,dt,'tustin');
[bN,aN] = tfdata(N_disc,'v');

% Filter histories
angle_x = theta_initial*[1 1];
angle_y = theta_initial*[1 1];

rate_x = [0 0];
rate_y = [0 0];
% Independent yaw notch memory (same filter coefficients)
yaw_angle_x = yaw(1)*[1 1];
yaw_angle_y = yaw(1)*[1 1];
yaw_rate_x = [0 0];
yaw_rate_y = [0 0];

%% Numerical stability reminder
% V7's 16-state Jacobian represented only the 2D pitch controller;
% a valid 3-axis Jacobian must contain all four fin servo states and the
% 3-axis rigid-body Euler dynamics. No 3D stability claim is made here.

%% Independent physical beam test (ONE-WAY: does not feed back into flight)
% EI is bending RIGIDITY (N m^2), not the previous K_flex compliance.
% Natural body bending frequency is derived from EI and modal mass.
% The old 4 Hz sensor-side oscillator remains a separate low-frequency DOF.
EI_body = 40000;             % N m^2, assumed structural bending rigidity
zeta_body = 0.03;            % provisional body modal damping ratio
x_sensor_beam = 0.55;        % m aft of nose, assumed sensor location
beam_force_4Hz_N = 15;      % N peak additional transverse fin force (test input)
beam_4Hz_start = 7;         % s, shaker interval start
beam_4Hz_duration = 4;      % s, smooth finite-duration shaker input
beam_4Hz_frequency = 4;     % Hz excitation, NOT airframe natural frequency
beam_yaw_shaker_fraction = 0.5; % second orthogonal forcing component

% First NON-RIGID free-free Euler-Bernoulli mode: cosh(B)*cos(B)=1.
beam_B = 4.730040744862704;
beam_x = linspace(0,L_body,801);
beam_b = beam_B/L_body;
beam_eta = (cosh(beam_B)-cos(beam_B)) ...
    /(sinh(beam_B)-sin(beam_B));
beam_phi_raw = cosh(beam_b*beam_x)+cos(beam_b*beam_x) ...
    -beam_eta*(sinh(beam_b*beam_x)+sin(beam_b*beam_x));
beam_norm = max(abs(beam_phi_raw));
beam_phi = beam_phi_raw/beam_norm; % dimensionless, max |phi| = 1
beam_dphi = beam_b*(sinh(beam_b*beam_x)-sin(beam_b*beam_x) ...
    -beam_eta*(cosh(beam_b*beam_x)+cos(beam_b*beam_x)))/beam_norm;
beam_ddphi = beam_b^2*(cosh(beam_b*beam_x)-cos(beam_b*beam_x) ...
    -beam_eta*(sinh(beam_b*beam_x)-sin(beam_b*beam_x)))/beam_norm;

% Rayleigh approximation for dry distributed mass + fuel point mass.
% The mode shape is held fixed as fuel burns: this is a test approximation.
beam_K = EI_body*trapz(beam_x,beam_ddphi.^2); % N/m generalized
beam_mass_shape = trapz(beam_x,beam_phi.^2);  % m
beam_phi_fuel = interp1(beam_x,beam_phi,x_fuel,'pchip');
beam_M_liftoff = (m_dry/L_body)*beam_mass_shape ...
    +m_fuel0*beam_phi_fuel^2;
beam_M_burnout = (m_dry/L_body)*beam_mass_shape;
beam_f_liftoff = sqrt(beam_K/beam_M_liftoff)/(2*pi);
beam_f_burnout = sqrt(beam_K/beam_M_burnout)/(2*pi);
beam_phi_nose = interp1(beam_x,beam_phi,x_CP_nose,'pchip');
beam_phi_cross = interp1(beam_x,beam_phi,x_CP_cross,'pchip');
beam_phi_tail = beam_phi(end);
beam_phi_fin0 = interp1(beam_x,beam_phi,x_fin_CP,'pchip');
beam_slope_sensor = interp1(beam_x,beam_dphi,x_sensor_beam,'pchip');
beam_Q_pitch = zeros(1,n); % generalized force (N), body-Y plane
beam_Q_yaw = zeros(1,n);   % generalized force (N), body-Z plane

%% Coupled flight and control simulation

for i = 1:n-1
    % All fins lock at zero after apogee; body continues passively in 3D.
    if ~descent_active && i>1 && t(i)>m_fuel0/mdot && ...
            V(i)<=0 && V(i-1)>0 && h(i)>0
        descent_active = true;
        apogee_lock_idx = i;
        fin_deflection(:,i) = 0;
        fin_command(:,i) = 0;
        delta(i) = 0; delta_yaw(i) = 0; delta_roll(i) = 0;
    end
    % Atmosphere and gust

    rho(i) = rho0*exp(-h(i)/H);

    % Two independent smooth horizontal gusts, same timing.
    wind_x = 0; wind_y = 0;
    if t(i)>=gust_start && t(i)<=gust_start+gust_duration
        phase = (t(i)-gust_start)/gust_duration;
        shape = 0.5*(1-cos(2*pi*phase));
        wind_x = wind_peak*shape;
        wind_y = wind_peak_y*shape;
    end
    wind_history(i) = wind_x;
    wind_y_history(i) = wind_y;
    v_air_x = vx(i)-wind_x;
    v_air_y = vy(i)-wind_y;
    speed_air = sqrt(v_air_x^2+v_air_y^2+V(i)^2);
    q_dyn(i) = 0.5*rho(i)*speed_air^2;

    % Body-to-world rotation and body-frame airflow.
    Rbw = rocketQuatToRotm(quat_bw(i,:));
    nose_world = Rbw(:,1);
    theta(i) = atan2(nose_world(1),nose_world(3));
    yaw(i) = atan2(nose_world(2),nose_world(3));
    roll_frame = rocketZeroRollFrame(nose_world);
    roll(i) = atan2(dot(roll_frame(:,3),Rbw(:,2)), ...
                    dot(roll_frame(:,2),Rbw(:,2)));
    v_body = (Rbw.'*[v_air_x;v_air_y;V(i)]).';
    gamma_i = atan2(v_air_x,V(i));
    gamma_y_i = atan2(v_air_y,V(i));
    alpha_i = atan2(-v_body(2),v_body(1));
    alpha_y_i = atan2(-v_body(3),v_body(1));
    flow_incidence_history(i) = atan2(norm(v_body(2:3)),v_body(1));

    % Propulsion and drag

    thrust = Isp*g*mdot*(m_fuel(i)>0);

    T_i = max(216.65,288.15-0.0065*h(i));

    a_i = sqrt(1.4*287.05*T_i);

    Mach_i = speed_air/a_i;

    % Mach-dependent aerodynamic stability (independent of drag model).
    Ma_eval = min(max(Mach_i,Ma_nodes(1)),Ma_nodes(end));
    CN_nose_i = CN_nose*ppval(pp_nose,Ma_eval);
    CN_fins_i = CN_fins*ppval(pp_fin,Ma_eval);
    CL_alpha_i = CL_alpha*ppval(pp_ctrl,Ma_eval);
    x_fin_CP_i = x_MAC_LE+c_mean*ppval(pp_ac,Ma_eval);
    x_CP_i = (CN_nose_i*x_CP_nose+CN_fins_i*x_fin_CP_i) ...
        /(CN_nose_i+CN_fins_i);

    CN_nose_hist(i) = CN_nose_i;
    CN_fins_hist(i) = CN_fins_i;
    CL_alpha_hist(i) = CL_alpha_i;
    x_fin_CP_hist(i) = x_fin_CP_i;
    x_CP_hist(i) = x_CP_i;
    l_fin(i) = x_fin_CP_i-x_CoM(i);
    static_margin(i) = (x_CP_i-x_CoM(i))/D_body;

    % Atmospheric viscosity and Reynolds numbers
    mu_i = 1.716e-5*(T_i/273.15)^1.5 ...
        *(273.15+110.4)/(T_i+110.4);

    Re_i = max(rho(i)*speed_air*drag_L/mu_i,1e4);

    % Skin-friction coefficients
    Cf_s = 1./(1.5*log(Re_i)-5.6).^2;
    Cf_r = Cf_rough0;

    if Mach_i < 1
        factor = 1-0.1*Mach_i^2;
        Cf_s = Cf_s*factor;
        Cf_r = Cf_r*factor;
    else
        Cf_s = Cf_s/(1+0.15*Mach_i^2)^0.58;
        Cf_r = Cf_r/(1+0.18*Mach_i^2);
    end

    Cf_i = max(Cf_s,Cf_r);

    D_friction_i = q_dyn(i)*sum( ...
        Cf_i.*drag_S.*drag_K);

    % Conical nose pressure drag
    if Mach_i <= 1
        Cd_nose_i = Cd_n0 ...
            +(Cd_n1-Cd_n0)*Mach_i^nose_exponent;

    elseif Mach_i < M_n2
        u = (Mach_i-1)/(M_n2-1);
        dM = M_n2-1;

        Cd_nose_i = ...
            (2*u^3-3*u^2+1)*Cd_n1 ...
            +(u^3-2*u^2+u)*dM*dCd_n1 ...
            +(-2*u^3+3*u^2)*Cd_n2 ...
            +(u^3-u^2)*dM*dCd_n2;
    else
        Cd_nose_i = 2.1*s_eps^2 ...
            +0.5*s_eps/sqrt(Mach_i^2-1);
    end

    D_nose_i = q_dyn(i)*A_ref*Cd_nose_i;

    % Four-fin leading-edge pressure drag
    if Mach_i < 0.9
        Cd_LE_i = (1-Mach_i^2)^(-0.417)-1;
    elseif Mach_i < 1
        Cd_LE_i = 1-1.785*(Mach_i-0.9);
    else
        Cd_LE_i = 1.214-0.502/Mach_i^2 ...
            +0.1095/Mach_i^4;
    end

    D_fin_pressure_i = q_dyn(i)*n_fins*A_fin_edge ...
        *Cd_LE_i*cos(sweep_LE)^2;

    % Engine-on / engine-off base drag
    if Mach_i <= 1
        Cd_base_i = 0.12+0.13*Mach_i^2;
    else
        Cd_base_i = 0.25/Mach_i;
    end

    A_base_i = max(A_base ...
        -A_motor*(m_fuel(i)>0),0);

    D_base_i = q_dyn(i)*Cd_base_i*A_base_i;

    % Rail-button drag
    if Mach_i < 1
        stag_i = 1+Mach_i^2/4+Mach_i^4/40;
    else
        stag_i = 1.84-0.76/Mach_i^2 ...
            +0.166/Mach_i^4+0.035/Mach_i^6;
    end

    D_buttons_i = n_rail_buttons*q_dyn(i) ...
        *A_button*0.85*stag_i;

    % Provisional fin-body interference drag
    D_interference_i = q_dyn(i)*n_fins*2 ...
        *Cf_i(3)*K_fin*R_body*fin_root;

    % Total geometry-based axial drag
    D_drag_i = D_friction_i+D_nose_i ...
        +D_fin_pressure_i+D_base_i ...
        +D_buttons_i+D_interference_i;

    % Individual fin/body aerodynamic forces and moments

theta_i = theta(i);
v_body_hat = v_body/max(speed_air,1e-9);
F_aero = -D_drag_i*v_body_hat;
M_body_aero = [0 0 0];
M_active_fin = [0 0 0];
% Reset aeroelastic force inputs each step; avoid reusing prior step values.
F_fin = zeros(n_fins,3);
F_nose = [0 0 0];
F_cross = [0 0 0];

if speed_air > 1 && V(i) > 0

    % Fin aerodynamic-center positions relative to CoM
    r_fin = [ ...
        repmat(x_CoM(i)-x_fin_CP_i,n_fins,1), ...
        fin_radial_cp*cos(fin_phi), ...
        fin_radial_cp*sin(fin_phi)];

    % Full 3D local flow: rolling fins create aerodynamic roll damping.
    w_body = [omega_roll(i) -omega_yaw(i) omega(i)];

    v_fin = repmat(v_body,n_fins,1) ...
        +cross(repmat(w_body,n_fins,1),r_fin,2);

    U_fin = vecnorm(v_fin,2,2);
    v_fin_hat = v_fin./max(U_fin,1e-9);
    q_fin = 0.5*rho(i)*U_fin.^2;
    M_fin = U_fin/a_i;

    % Individual fin incidence and ACTUAL, mixed fin deflections.
    alpha_fin_flow = atan2( ...
        sum(v_fin.*fin_normal,2),v_fin(:,1));
    delta_fin = fin_deflection(:,i);

    CL_passive = (A_ref*CN_fins_i/S_fin)*alpha_fin_flow;
    CL_control = CL_max*tanh( ...
        CL_alpha_i*delta_fin/CL_max);

    CL_fin = CL_passive+CL_control;

    % Lift perpendicular to local airflow
    lift_dir = cross(v_fin_hat,fin_span_dir,2);
    lift_dir = lift_dir./max(vecnorm(lift_dir,2,2),1e-9);

    F_lift = -(q_fin*A_fin.*CL_fin).*lift_dir;
    F_control = -(q_fin*A_fin.*CL_control).*lift_dir;

    % Individual skin friction
    Re_fin = max(rho(i)*U_fin*c_mean/mu_i,1e4);

    Cf_s = 1./(1.5*log(Re_fin)-5.6).^2;
    Cf_r = Cf_rough0(3)*ones(n_fins,1);

    sub = M_fin < 1;

    Cf_s(sub) = Cf_s(sub).*(1-0.1*M_fin(sub).^2);
    Cf_r(sub) = Cf_r(sub).*(1-0.1*M_fin(sub).^2);

    Cf_s(~sub) = Cf_s(~sub) ...
        ./(1+0.15*M_fin(~sub).^2).^0.58;

    Cf_r(~sub) = Cf_r(~sub) ...
        ./(1+0.18*M_fin(~sub).^2);

    Cf_fin = max(Cf_s,Cf_r);

    D_skin = q_fin*(2*A_fin).*Cf_fin*K_fin;

    % Individual leading-edge pressure drag
    Cd_LE = zeros(n_fins,1);

    a = M_fin < 0.9;
    b = M_fin >= 0.9 & M_fin < 1;
    c = M_fin >= 1;

    Cd_LE(a) = (1-M_fin(a).^2).^(-0.417)-1;
    Cd_LE(b) = 1-1.785*(M_fin(b)-0.9);
    Cd_LE(c) = 1.214-0.502./M_fin(c).^2 ...
        +0.1095./M_fin(c).^4;

    D_pressure = q_fin*A_fin_edge ...
        .*Cd_LE*cos(sweep_LE)^2;

    % Induced and interference drag
    D_induced = q_fin*A_fin ...
        .*CL_fin.^2/(pi*e_fin*fin_AR);

    D_interference_fin = q_fin*2 ...
        .*Cf_fin*K_fin*R_body*fin_root;

    % Replace aggregate fin drag to prevent double-counting
    D_fin_old = q_dyn(i)*Cf_i(3)*drag_S(3)*drag_K(3) ...
        +D_fin_pressure_i+D_interference_i;

    D_nonfin = max(D_drag_i-D_fin_old,0);

    D_fin_new = D_skin+D_pressure ...
        +D_induced+D_interference_fin;

    F_fin = F_lift-D_fin_new.*v_fin_hat;

    % Generalized nose normal force, perpendicular to the local airflow.
    % For yaw=0 this reduces to the original 2D small-angle model.
    normal_to_flow = [1 0 0] - v_body_hat(1)*v_body_hat;
    normal_mag = norm(normal_to_flow);
    alpha_total = atan2(norm(v_body(2:3)),v_body(1));
    if normal_mag > 1e-9
        F_nose = q_dyn(i)*A_ref*CN_nose_i*alpha_total ...
            *normal_to_flow/normal_mag;
    else
        F_nose = [0 0 0];
    end
    r_nose = [x_CoM(i)-x_CP_nose,0,0];

    F_cross = [0, ...
        -q_dyn(i)*Cd_body_cross*A_side_body* ...
           v_body_hat(2)*abs(v_body_hat(2)), ...
        -q_dyn(i)*Cd_body_cross*A_side_body* ...
           v_body_hat(3)*abs(v_body_hat(3))];
    r_cross = [x_CoM(i)-x_CP_cross,0,0];

    F_aero = -D_nonfin*v_body_hat ...
        +F_nose+F_cross+sum(F_fin,1);
    M_body_aero = sum(cross(r_fin,F_fin,2),1) ...
        +cross(r_nose,F_nose)+cross(r_cross,F_cross);
    M_active_fin = sum(cross(r_fin,F_control,2),1);

elseif descent_active && speed_air > 1
    % Post-apogee: aerodynamic loading with unpowered, NEUTRAL static fins.
    % Unlike the ascent model, the axial air velocity can be negative
    % (tail-first motion). Avoid the invalid ~180-degree linear nose AoA.
    % The finite-flow normal load opposes the velocity through the fin plane.
    w_body = [omega_roll(i) -omega_yaw(i) omega(i)];
    nose_first_weight = 0.5*(1+tanh(5*v_body_hat(1)));
    x_fin_reverse = x_MAC_LE+0.75*c_mean;
    x_fin_descent = nose_first_weight*x_fin_CP_i ...
        +(1-nose_first_weight)*x_fin_reverse;

    r_fin = [ ...
        repmat(x_CoM(i)-x_fin_descent,n_fins,1), ...
        fin_radial_cp*cos(fin_phi), ...
        fin_radial_cp*sin(fin_phi)];
    v_fin = repmat(v_body,n_fins,1) ...
        +cross(repmat(w_body,n_fins,1),r_fin,2);
    U_fin = vecnorm(v_fin,2,2);
    q_fin = 0.5*rho(i)*U_fin.^2;
    v_fin_hat = v_fin./max(U_fin,1e-9);

    % Passive normal force from each fixed fin (no commanded lift).
    normal_flow = sum(v_fin.*fin_normal,2)./max(U_fin,1e-9);
    CL_static = (A_ref*CN_fins_i/S_fin)*normal_flow;
    CL_static = max(-CL_max,min(CL_max,CL_static));
    F_fin_normal = -(q_fin*A_fin.*CL_static).*fin_normal;
    D_fin_induced = q_fin*A_fin.*CL_static.^2 ...
        /(pi*e_fin*fin_AR);
    F_fin_induced = -D_fin_induced.*v_fin_hat;
    F_fin = F_fin_normal+F_fin_induced;

    % Axial drag varies smoothly between nose-first and tail-first flow.
    % Descent Cd is assumed; it is not the ascent nose-pressure correlation.
    D_reverse = D_friction_i+D_fin_pressure_i+D_buttons_i ...
        +D_interference_i+q_dyn(i)*Cd_tail_first*A_ref;
    D_axial = nose_first_weight*D_drag_i ...
        +(1-nose_first_weight)*D_reverse;

    % Two transverse components of broadside crossflow drag.
    F_cross = [0, ...
        -q_dyn(i)*Cd_body_cross*A_side_body* ...
           v_body_hat(2)*abs(v_body_hat(2)), ...
        -q_dyn(i)*Cd_body_cross*A_side_body* ...
           v_body_hat(3)*abs(v_body_hat(3))];
    r_cross = [x_CoM(i)-x_CP_cross,0,0];

    F_aero = -D_axial*v_body_hat+F_cross+sum(F_fin,1);
    
M_fin_passive = sum(cross(r_fin,F_fin,2),1);
M_cross_passive = cross(r_cross,F_cross);

M_body_aero = M_fin_passive + M_cross_passive;

if abs(t(i)-39.06) <= dt/2
    fprintf('\nDESCENT MOMENT BREAKDOWN AT %.3f s\n',t(i));

    fprintf('Fin pitch / yaw: %.2f / %.2f Nm\n', ...
        M_fin_passive(3),-M_fin_passive(2));

    fprintf('Body pitch / yaw: %.2f / %.2f Nm\n', ...
        M_cross_passive(3),-M_cross_passive(2));

    fprintf('Body-frame velocity X/Y/Z: %.2f / %.2f / %.2f m/s\n', ...
        v_body(1),v_body(2),v_body(3));
end

    M_active_fin = [0 0 0];
end

% Project actual modeled transverse aerodynamic loads onto the bending
% mode. This is ONE-WAY forcing; beam deformation does not change fin AoA.
% Fin longitudinal CP varies with Mach, so its modal lever changes too.
if speed_air > 1
    phi_fin_i = interp1(beam_x,beam_phi,x_fin_CP_i,'pchip');
    if descent_active
        phi_fin_i = interp1(beam_x,beam_phi,x_fin_descent,'pchip');
    end
    beam_Q_pitch(i) = phi_fin_i*sum(F_fin(:,2)) ...
        +beam_phi_nose*F_nose(2)+beam_phi_cross*F_cross(2);
    beam_Q_yaw(i) = phi_fin_i*sum(F_fin(:,3)) ...
        +beam_phi_nose*F_nose(3)+beam_phi_cross*F_cross(3);
end

% Actual drag opposing the air-relative velocity
drag_force(i) = -dot(F_aero,v_body_hat);

    % 3D translational dynamics in the inertial world frame.
    % Body X receives thrust. Aerodynamic forces are body-frame vectors.
    F_world = Rbw*([thrust;0;0]+F_aero.');
    accel_world = F_world/m_total(i)+[0;0;-g];
    ax_history(i) = accel_world(1);
    ay_history(i) = accel_world(2);
    az_history(i) = accel_world(3);

    vx(i+1) = vx(i)+accel_world(1)*dt;
    vy(i+1) = vy(i)+accel_world(2)*dt;
    V(i+1) = V(i)+accel_world(3)*dt;
    x_pos(i+1) = x_pos(i)+vx(i+1)*dt;
    y_pos(i+1) = y_pos(i)+vy(i+1)*dt;
    h_trial = h(i)+V(i+1)*dt;
    h(i+1) = max(0,h_trial);

    if descent_active && h_trial<=0 && h(i)>0
        impact_idx = i+1;
        impact_fraction = min(max(h(i)/(h(i)-h_trial),0),1);
        impact_time = t(i)+impact_fraction*dt;
        x_pos(i+1) = x_pos(i)+impact_fraction*(x_pos(i+1)-x_pos(i));
        y_pos(i+1) = y_pos(i)+impact_fraction*(y_pos(i+1)-y_pos(i));
        vx(i+1) = vx(i)+impact_fraction*(vx(i+1)-vx(i));
        vy(i+1) = vy(i)+impact_fraction*(vy(i+1)-vy(i));
        V(i+1) = V(i)+impact_fraction*(V(i+1)-V(i));
    elseif h_trial<=0 && V(i+1)<0
        V(i+1) = 0; vx(i+1) = 0; vy(i+1) = 0;
    end

    % Two matching bending-compensated sensors and notch filters.
    theta_sensor = theta(i)+q_flex(i);
    omega_sensor = omega(i)+qdot_flex(i);
    yaw_sensor = yaw(i)+q_flex_yaw(i);
    yaw_rate_sensor = omega_yaw(i)+qdot_flex_yaw(i);
    switch sensor_mode
        case 0
            theta_raw = theta(i); omega_raw = omega(i);
            yaw_raw = yaw(i); yaw_rate_raw = omega_yaw(i);
        case 1
            theta_raw = theta_sensor; omega_raw = omega_sensor;
            yaw_raw = yaw_sensor; yaw_rate_raw = yaw_rate_sensor;
        case 2
            theta_raw = theta_sensor-q_flex_est(i);
            omega_raw = omega_sensor-qdot_flex_est(i);
            yaw_raw = yaw_sensor-q_flex_yaw_est(i);
            yaw_rate_raw = yaw_rate_sensor-qdot_flex_yaw_est(i);
    end
    theta_m = bN(1)*theta_raw+bN(2)*angle_x(1)+bN(3)*angle_x(2) ...
        -aN(2)*angle_y(1)-aN(3)*angle_y(2);
    omega_m = bN(1)*omega_raw+bN(2)*rate_x(1)+bN(3)*rate_x(2) ...
        -aN(2)*rate_y(1)-aN(3)*rate_y(2);
    yaw_m = bN(1)*yaw_raw+bN(2)*yaw_angle_x(1)+bN(3)*yaw_angle_x(2) ...
        -aN(2)*yaw_angle_y(1)-aN(3)*yaw_angle_y(2);
    yaw_rate_m = bN(1)*yaw_rate_raw+bN(2)*yaw_rate_x(1)+bN(3)*yaw_rate_x(2) ...
        -aN(2)*yaw_rate_y(1)-aN(3)*yaw_rate_y(2);
    angle_x=[theta_raw angle_x(1)]; angle_y=[theta_m angle_y(1)];
    rate_x=[omega_raw rate_x(1)]; rate_y=[omega_m rate_y(1)];
    yaw_angle_x=[yaw_raw yaw_angle_x(1)];
    yaw_angle_y=[yaw_m yaw_angle_y(1)];
    yaw_rate_x=[yaw_rate_raw yaw_rate_x(1)];
    yaw_rate_y=[yaw_rate_m yaw_rate_y(1)];

    % ROLL-AWARE tilt feedback: world-frame pitch/yaw corrections are
    % rotated into the current BODY fin axes. This is essential at 180 roll;
    % naïvely reusing the v10 commands would steer in the wrong directions.
    % A small-angle world rotation vector [yaw,-pitch,0] corrects the nose.
    world_tilt_correction = [yaw_m;-theta_m;0];
    body_tilt_correction = Rbw.'*world_tilt_correction;
    cmd_pitch = Kp_fin*body_tilt_correction(3)-Kd_fin*omega_m;
    cmd_yaw = -Kp_fin*body_tilt_correction(2)-Kd_fin*yaw_rate_m;

    % Roll coordinate is twist around the rocket nose, not the yaw/pitch
    % Euler angle. Choose the shortest path and break the 180-degree tie.
    roll_error = atan2(sin(roll(i)),cos(roll(i)));
    if abs(abs(roll_error)-pi)<1e-5
        if abs(omega_roll(i))>deg2rad(0.1)
            selected_direction = sign(omega_roll(i));
        else
            selected_direction = roll_tie_direction;
        end
        roll_error = -selected_direction*pi;
    end
    roll_error_history(i) = roll_error;
    roll_control_sign(i) = -sign(roll_error);

    % Required roll moment, then aerodynamic inversion to an equivalent
    % deflection. At low dynamic pressure roll authority is unavailable.
    M_roll_req = I_roll(i)*(-roll_wn^2*roll_error ...
                       -2*roll_zeta*roll_wn*omega_roll(i));
    roll_requested_moment(i) = M_roll_req;
    K_roll_def = -q_dyn(i)*n_fins*A_fin*CL_alpha_i*fin_radial_cp;
    cmd_roll = 0;
    if enable_roll && q_dyn(i)>=roll_q_enable
        cmd_roll = M_roll_req/K_roll_def;
        cmd_roll = max(-roll_max_contribution, ...
                   min(roll_max_contribution,cmd_roll));
    end
    if ~(V(i)>1 && speed_air>1 && ~descent_active)
        cmd_pitch=0; cmd_yaw=0; cmd_roll=0;
    end
    if ~enable_yaw, cmd_yaw=0; end
    if ~enable_roll, cmd_roll=0; end

    % Four-fin mixer: all axes share every servo. Uniform scaling preserves
    % the requested pitch/yaw/roll proportions when any fin reaches 5 deg.
    raw_fin_cmd = fin_mix*[cmd_pitch;cmd_yaw;cmd_roll];
    mix_scale = min(1,delta_max/max(max(abs(raw_fin_cmd)),eps));
    fin_mixer_scale(i) = mix_scale;
    fin_target = raw_fin_cmd*mix_scale;
    if descent_active
        fin_target(:)=0;
        fin_deflection(:,i+1)=0;  % instant mechanical lock at apogee
    else
        fin_deflection(:,i+1)=fin_deflection(:,i) ...
            +(fin_target-fin_deflection(:,i))*dt/tau_s;
        fin_deflection(:,i+1)=max(-delta_max, ...
                              min(delta_max,fin_deflection(:,i+1)));
    end
    fin_command(:,i) = fin_target;
    achieved_mix = fin_mix_inv*fin_deflection(:,i+1);
    delta(i+1)=achieved_mix(1);
    delta_yaw(i+1)=achieved_mix(2);
    delta_roll(i+1)=achieved_mix(3);

    % Keep signed pitch/yaw convention consistent with +5/+5 offsets.
    torque_fin(i)=M_active_fin(3);
    torque_fin_yaw(i)=-M_active_fin(2);
    torque_total(i)=M_body_aero(3);
    torque_total_yaw(i)=-M_body_aero(2);
    torque_roll(i)=M_body_aero(1);
    active_torque_roll(i)=M_active_fin(1);

    % Full axisymmetric Euler rigid-body dynamics, including the gyroscopic
    % coupling between roll and transverse rates; roll inertia is distinct.
    wbody=[omega_roll(i),-omega_yaw(i),omega(i)];
    Ibody=[I_roll(i),I_total(i),I_total(i)];
    gyro=cross(wbody,Ibody.*wbody);
    wdot=(M_body_aero-gyro)./Ibody;
    omega_roll_dot(i)=wdot(1);
    omega_roll(i+1)=omega_roll(i)+wdot(1)*dt;
    omega_yaw(i+1)=omega_yaw(i)-wdot(2)*dt;
    omega(i+1)=omega(i)+wdot(3)*dt;
    wbody_next=[omega_roll(i+1),-omega_yaw(i+1),omega(i+1)];
    roll_spin(i+1)=roll_spin(i)+omega_roll(i+1)*dt;
    dq=0.5*dt*rocketQuatMultiply(quat_bw(i,:),[0 wbody_next]);
    quat_bw(i+1,:)=quat_bw(i,:)+dq;
    quat_bw(i+1,:)=quat_bw(i+1,:)/norm(quat_bw(i+1,:));
    R_next=rocketQuatToRotm(quat_bw(i+1,:));
    nose_next=R_next(:,1);
    theta(i+1)=atan2(nose_next(1),nose_next(3));
    yaw(i+1)=atan2(nose_next(2),nose_next(3));
    Rzero_next=rocketZeroRollFrame(nose_next);
    roll(i+1)=atan2(dot(Rzero_next(:,3),R_next(:,2)), ...
                       dot(Rzero_next(:,2),R_next(:,2)));

    % True and estimated bending

    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;

    % Estimate fin torque from aerodynamic conditions
    % For mixed fins the scalar old fin-pair formula is invalid.
    % Use the modeled active fin torque (computed from individual fin forces).
    torque_est = torque_fin(i);

    torque_est_history(i) = torque_est;

    % Estimated flexible-mode dynamics
    qdd_est = wn_flex_est^2 ...
        *(K_flex_est*torque_est-q_flex_est(i)) ...
        -2*zeta_flex_est*wn_flex_est*qdot_flex_est(i);

    qdot_flex_est(i+1) = qdot_flex_est(i)+qdd_est*dt;
    q_flex_est(i+1) = q_flex_est(i) ...
        +qdot_flex_est(i+1)*dt;

    % Orthogonal yaw bending mode and independent estimator.
    qdd_yaw = wn_flex^2*(K_flex*torque_fin_yaw(i)-q_flex_yaw(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex_yaw(i);
    qdot_flex_yaw(i+1) = qdot_flex_yaw(i)+qdd_yaw*dt;
    q_flex_yaw(i+1) = q_flex_yaw(i)+qdot_flex_yaw(i+1)*dt;
    torque_est_yaw = torque_fin_yaw(i);
    torque_est_yaw_history(i) = torque_est_yaw;
    qdd_yaw_est = wn_flex_est^2 ...
        *(K_flex_est*torque_est_yaw-q_flex_yaw_est(i)) ...
        -2*zeta_flex_est*wn_flex_est*qdot_flex_yaw_est(i);
    qdot_flex_yaw_est(i+1) = qdot_flex_yaw_est(i)+qdd_yaw_est*dt;
    q_flex_yaw_est(i+1) = q_flex_yaw_est(i)+qdot_flex_yaw_est(i+1)*dt;
    if isfinite(impact_idx)
        break
    end
end

% Discard the unused remainder of the 120-second simulation after impact.
if isfinite(impact_idx)
    stop_idx = impact_idx;
    t=t(1:stop_idx); n=numel(t);
    m_fuel=m_fuel(1:n); m_total=m_total(1:n);
    x_CoM=x_CoM(1:n); I_total=I_total(1:n);
    V=V(1:n); h=h(1:n); vx=vx(1:n); x_pos=x_pos(1:n);
    vy=vy(1:n); y_pos=y_pos(1:n);
    theta=theta(1:n); omega=omega(1:n); delta=delta(1:n);
    yaw=yaw(1:n); omega_yaw=omega_yaw(1:n); delta_yaw=delta_yaw(1:n);
    roll=roll(1:n); omega_roll=omega_roll(1:n); delta_roll=delta_roll(1:n);
    roll_spin=roll_spin(1:n); omega_roll_dot=omega_roll_dot(1:n);
    fin_deflection=fin_deflection(:,1:n); fin_command=fin_command(:,1:n);
    fin_mixer_scale=fin_mixer_scale(1:n);
    roll_error_history=roll_error_history(1:n);
    roll_requested_moment=roll_requested_moment(1:n);
    roll_control_sign=roll_control_sign(1:n);
    quat_bw=quat_bw(1:n,:);
    flow_incidence_history=flow_incidence_history(1:n);
    q_flex=q_flex(1:n); qdot_flex=qdot_flex(1:n);
    q_flex_est=q_flex_est(1:n); qdot_flex_est=qdot_flex_est(1:n);
    q_flex_yaw=q_flex_yaw(1:n); qdot_flex_yaw=qdot_flex_yaw(1:n);
    q_flex_yaw_est=q_flex_yaw_est(1:n);
    qdot_flex_yaw_est=qdot_flex_yaw_est(1:n);
    rho=rho(1:n); q_dyn=q_dyn(1:n);
    drag_force=drag_force(1:n); wind_history=wind_history(1:n);
    wind_y_history=wind_y_history(1:n);
    ax_history=ax_history(1:n); ay_history=ay_history(1:n);
    az_history=az_history(1:n);
    torque_fin=torque_fin(1:n); torque_total=torque_total(1:n);
    torque_fin_yaw=torque_fin_yaw(1:n);
    torque_total_yaw=torque_total_yaw(1:n);
    torque_roll=torque_roll(1:n);
    active_torque_roll=active_torque_roll(1:n);
    torque_est_history=torque_est_history(1:n);
    torque_est_yaw_history=torque_est_yaw_history(1:n);
    l_fin=l_fin(1:n); static_margin=static_margin(1:n);
    x_CP_hist=x_CP_hist(1:n); x_fin_CP_hist=x_fin_CP_hist(1:n);
    CN_nose_hist=CN_nose_hist(1:n); CN_fins_hist=CN_fins_hist(1:n);
    CL_alpha_hist=CL_alpha_hist(1:n);
end

% Last sample has no integration step: carry over stability diagnostics.
x_CP_hist(end) = x_CP_hist(end-1);
x_fin_CP_hist(end) = x_fin_CP_hist(end-1);
CN_nose_hist(end) = CN_nose_hist(end-1);
CN_fins_hist(end) = CN_fins_hist(end-1);
CL_alpha_hist(end) = CL_alpha_hist(end-1);
l_fin(end) = x_fin_CP_hist(end)-x_CoM(end);
static_margin(end) = (x_CP_hist(end)-x_CoM(end))/D_body;

% Complete the final atmosphere/dynamic-pressure sample.
rho(end) = rho0*exp(-h(end)/H);
q_dyn(end) = 0.5*rho(end)*( ...
    (vx(end)-wind_history(end))^2+ ...
    (vy(end)-wind_y_history(end))^2+V(end)^2);

% The terminal sample is the interpolated pre-impact state; no contact
% acceleration/impact impulse is simulated. Carry the last pre-impact
% acceleration and drag sample into plots (no further state integration).
ax_history(end) = ax_history(end-1);
ay_history(end) = ay_history(end-1);
az_history(end) = az_history(end-1);
drag_force(end) = drag_force(end-1);

%% Summary diagnostics

burn_time = m_fuel0/mdot;
burn_idx = find(t <= burn_time);
[apogee,idx_apogee] = max(h);
[maxQ,k] = max(q_dyn(burn_idx));
idx_maxQ = burn_idx(k);

air_speed = sqrt((vx-wind_history).^2+ ...
    (vy-wind_y_history).^2+V.^2);
air_temp = max(216.65,288.15-0.0065*h);
Mach = air_speed./sqrt(1.4*287.05*air_temp);
[M_max,idx_M] = max(Mach(1:idx_apogee));

gamma = atan2(vx-wind_history,V);
gamma_y = atan2(vy-wind_y_history,V);
alpha_AoA = flow_incidence_history; % 3D unsigned total incidence
alpha_AoA(end) = alpha_AoA(end-1);
valid_AoA = t <= burn_time & V > 0 & air_speed > 20;
gust_idx = t >= gust_start & t <= gust_start+gust_duration;

% World-frame net (kinematic) acceleration, and proper acceleration
% from thrust + aerodynamics only: what an ideal onboard IMU would read.
accel_vertical_g = az_history/g;
accel_net_g = sqrt(ax_history.^2+ay_history.^2+az_history.^2)/g;
accel_proper_g = sqrt(ax_history.^2+ay_history.^2 ...
    +(az_history+g).^2)/g;

% Air-relative rocket incidence. 0 deg = nose-first along airflow,
% 180 deg = tail-first, and 90 deg = broadside. This remains valid when
% descending; small-angle C_N-alpha assumptions do not.
flow_incidence_deg = rad2deg(alpha_AoA);

% Effective drag referenced to frontal body area, not a drag correlation.
% Suppress near-apogee low-q division because Cd becomes ill-conditioned.
Cd_effective = NaN(1,n);
valid_cd = q_dyn >= 100;            % Pa threshold for diagnostics
Cd_effective(valid_cd) = ...
    drag_force(valid_cd)./(q_dyn(valid_cd)*A_ref);

fprintf('\nFLIGHT AND CONTROL SUMMARY\n');
fprintf('Apogee: %.1f m at %.1f s\n',apogee,t(idx_apogee));
fprintf('Maximum dynamic pressure: %.1f kPa at %.2f s\n',maxQ/1000,t(idx_maxQ));
fprintf('Maximum Mach: %.3f at %.2f s\n',M_max,t(idx_M));
fprintf('Downrange at apogee X/Y: %.2f / %.2f m\n', ...
    x_pos(idx_apogee),y_pos(idx_apogee));
fprintf('Horizontal downrange magnitude: %.2f m\n', ...
    hypot(x_pos(idx_apogee),y_pos(idx_apogee)));
fprintf('\nBALLISTIC DESCENT (FINS LOCKED NEUTRAL, NO PARACHUTE)\n');
if isfinite(apogee_lock_idx)
    fprintf('Fin lock / apogee time: %.3f s\n',t(apogee_lock_idx));
end
if isfinite(impact_idx)
    fprintf('Touchdown time: %.3f s\n',impact_time);
    fprintf('Time from apogee to touchdown: %.2f s\n', ...
        impact_time-t(idx_apogee));
    fprintf('Touchdown position X/Y: %.2f / %.2f m\n',x_pos(end),y_pos(end));
    fprintf('Touchdown vertical velocity: %.2f m/s\n',V(end));
    fprintf('Touchdown velocity X/Y: %.2f / %.2f m/s\n',vx(end),vy(end));
    fprintf('Touchdown speed: %.2f m/s\n', ...
        sqrt(vx(end)^2+vy(end)^2+V(end)^2));
    fprintf('Peak descent dynamic pressure: %.1f kPa\n', ...
        max(q_dyn(idx_apogee:end))/1000);
else
    fprintf('No touchdown within simulated %.1f s; extend t if needed.\n',t(end));
end
if isfinite(apogee_lock_idx)
    fprintf('Max fin deflection after lock (pitch/yaw): %.6f / %.6f deg\n', ...
        max(abs(rad2deg(delta(apogee_lock_idx:end)))), ...
        max(abs(rad2deg(delta_yaw(apogee_lock_idx:end)))));
end
fprintf('Maximum powered AoA: %.3f deg\n',max(abs(rad2deg(alpha_AoA(valid_AoA)))));
fprintf('Pitch / yaw at burnout: %.3f / %.3f deg\n', ...
    rad2deg(theta(burn_idx(end))),rad2deg(yaw(burn_idx(end))));
fprintf('Maximum pitch / yaw during gust: %.3f / %.3f deg\n', ...
    max(abs(rad2deg(theta(gust_idx)))), ...
    max(abs(rad2deg(yaw(gust_idx)))));
fprintf('Peak powered fin deflection pitch / yaw: %.3f / %.3f deg\n', ...
    max(abs(rad2deg(delta(burn_idx)))), ...
    max(abs(rad2deg(delta_yaw(burn_idx)))));
fprintf('Maximum fin deflection during gust: %.3f deg\n',max(abs(rad2deg(delta(gust_idx)))));
fprintf('Peak control torque during gust: %.2f N m\n',max(abs(torque_fin(gust_idx))));
fprintf('Maximum structural bending during gust: %.3f deg\n',max(abs(rad2deg(q_flex(gust_idx)))));
fprintf('RMS bending estimation error during gust: %.3f deg\n', ...
    rad2deg(sqrt(mean((q_flex(gust_idx)-q_flex_est(gust_idx)).^2))));
fprintf('Sensor mode: %d, Mach-dependent stability: %d\n', ...
    sensor_mode,enable_mach_aero);
fprintf('Drag at maximum Q: %.2f N, effective Cd: %.4f\n', ...
    drag_force(idx_maxQ),drag_force(idx_maxQ)/max(maxQ*A_ref,1e-9));
fprintf('CP: %.3f m, initial/burnout static margin: %.2f / %.2f calibers\n', ...
    x_CP_hist(idx_maxQ),static_margin(1),static_margin(burn_idx(end)));

% Descent-only aerodynamic outputs; coefficients remain provisional.
d_idx = idx_apogee:n;
d_cd = d_idx(isfinite(Cd_effective(d_idx)));
fprintf('\nDESCENT AERODYNAMIC DIAGNOSTIC\n');
fprintf('Maximum descent Mach: %.3f\n',max(Mach(d_idx)));
if ~isempty(d_cd)
    fprintf('Mean reference-area Cd for q >= 100 Pa: %.3f\n', ...
        mean(Cd_effective(d_cd)));
    fprintf('Effective Cd immediately before impact: %.3f\n', ...
        Cd_effective(d_cd(end)));
end
fprintf('Flow incidence before impact: %.2f deg (0 nose-first, 180 tail-first)\n', ...
    flow_incidence_deg(end));
fprintf('Peak passive pitch / yaw / roll moment during descent: %.2f / %.2f / %.2f N m\n', ...
    max(abs(torque_total(d_idx))),max(abs(torque_total_yaw(d_idx))), ...
    max(abs(torque_roll(d_idx))));
fprintf('Peak absolute pitch / yaw / roll rate during descent: %.2f / %.2f / %.2f deg/s\n', ...
    max(abs(rad2deg(omega(d_idx)))), ...
    max(abs(rad2deg(omega_yaw(d_idx)))), ...
    max(abs(rad2deg(omega_roll(d_idx)))));
fprintf('Peak sensed/proper acceleration in descent: %.2f g\n', ...
    max(accel_proper_g(d_idx)));
fprintf('Impact shock/contact forces NOT modeled.\n');

%% Full-flight acceleration (FIRST graph)
figure('Name','Full-Flight Acceleration in g')
plot(t,accel_vertical_g,'LineWidth',1.2)
hold on
plot(t,accel_net_g,'LineWidth',1.2)
plot(t,accel_proper_g,'LineWidth',1.2)
xline(t(idx_apogee),':','Apogee')
if isfinite(impact_idx)
    xline(impact_time,':','Touchdown')
end
legend('Vertical net acceleration (signed)', ...
    'Magnitude of net acceleration', ...
    'Proper g-load: thrust + aerodynamics', 'Location','best')
xlabel('Flight time (s)')
ylabel('Acceleration (g, g = 9.81 m/s^2)')
title('Rocket Acceleration from Liftoff to Touchdown (No Impact Shock)')
grid on

%% Original graphs (all 22 retained)
% These plots reproduce the original diagnostic figures. Reference-only
% correlations and the legacy linear Bode approximation are not used by
% the nonlinear flight/control model.

% 1. Rocket Altitude
figure
plot(t,h)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('Rocket Altitude')
grid on

% 2. Rocket Velocity
figure
plot(t,V)
xlabel('Time (s)')
ylabel('Vertical Velocity (m/s)')
title('Rocket Velocity')
grid on

% 3. Dynamic Pressure
figure
plot(t,q_dyn/1000)
xlabel('Time (s)')
ylabel('Dynamic Pressure (kPa)')
title('Dynamic Pressure')
grid on

% 4. Downrange position
figure
plot(t,x_pos)
xlabel('Time (s)')
ylabel('Horizontal Position (m)')
title('Rocket Downrange X Position')
xlim([0 t(end)])
grid on

% 5. Powered-flight AoA
figure
plot(t(valid_AoA),rad2deg(alpha_AoA(valid_AoA)))
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Angle of Attack During Powered Flight')
xlim([0 burn_time])
grid on

% 6. Attitude and air-relative flight path
figure
plot(t(burn_idx),rad2deg(theta(burn_idx)), ...
    t(burn_idx),rad2deg(yaw(burn_idx)), ...
    t(burn_idx),rad2deg(gamma(burn_idx)), ...
    t(burn_idx),rad2deg(gamma_y(burn_idx)), ...
    'LineWidth',1.2)
legend('Rocket pitch','Rocket yaw','Flow angle X','Flow angle Y')
xlabel('Time (s)')
ylabel('Angle from vertical (deg)')
title('Powered Ascent: Attitude vs Trajectory')
grid on

% 7. True rocket pitch, yaw and roll (roll shown wrapped)
figure
plot(t,rad2deg(theta),t,rad2deg(yaw),t,rad2deg(roll),'LineWidth',1.2)
legend('Pitch','Yaw','Roll (wrapped)')
xlabel('Time (s)')
ylabel('Angle (deg)')
title('Rocket Pitch, Yaw and Roll')
xlim([0 15])
grid on

% 8. Achieved mixer components (four actual fins plotted separately below)
figure
plot(t,rad2deg(delta),t,rad2deg(delta_yaw), ...
     t,rad2deg(delta_roll),'LineWidth',1.2)
legend('Pitch mixer','Yaw mixer','Roll mixer')
xlabel('Time (s)')
ylabel('Fin Deflection (deg)')
title('Pitch, Yaw and Roll Mixer Components')
xlim([0 15])
grid on

% 9. Bending state
figure
plot(t,rad2deg(q_flex),t,rad2deg(q_flex_yaw),'LineWidth',1.2)
legend('Pitch bending','Yaw bending')
xlabel('Time (s)')
ylabel('Structural Bending (deg)')
title('Structural Bending')
xlim([0 15])
grid on

% 10. Pitch moment components
figure
plot(t,torque_fin,t,torque_fin_yaw,t,active_torque_roll, ...
    t,torque_total,t,torque_total_yaw,t,torque_roll)
legend('Active pitch','Active yaw','Active roll', ...
       'Total pitch','Total yaw','Total roll')
xlabel('Time (s)')
ylabel('Moment (N m)')
title('3-Axis Aerodynamic Moments')
xlim([0 5])
grid on

% 11. Frozen-time reduced-order stability estimate at Max Q.
% Captures Mach-varying stiffness, fin authority, fin pitch-rate damping,
% actuator, structural mode, estimator and measurement notch. Still an
% approximation: not a full numerical Jacobian of the nonlinear model.
s = tf('s');
C = Kp_fin+Kd_fin*s;
Gs = 1/(tau_s*s+1);
Gflex_true = K_flex*wn_flex^2 / ...
    (s^2+2*zeta_flex*wn_flex*s+wn_flex^2);
Gflex_est = K_flex_est*wn_flex_est^2 / ...
    (s^2+2*zeta_flex_est*wn_flex_est*s+wn_flex_est^2);

if sensor_mode==0
    Gflex_residual = 0;
elseif sensor_mode==1
    Gflex_residual = Gflex_true;
else
    Gflex_residual = Gflex_true-Gflex_est;
end

Ka_hist = q_dyn.*S_fin.*CL_alpha_hist.*l_fin;
Kpass_hist = q_dyn.*A_ref.*(CN_nose_hist+CN_fins_hist) ...
    .*(x_CP_hist-x_CoM);
Kq_hist = q_dyn.*A_ref.*CN_fins_hist.*l_fin.^2 ...
    ./max(air_speed,1);

Gr = 1/(I_total(idx_maxQ)*s^2 + Kq_hist(idx_maxQ)*s ...
    +Kpass_hist(idx_maxQ));
L_flexible = C*Gs*Ka_hist(idx_maxQ)*(Gr+Gflex_residual)*N_cont;
figure
margin(L_flexible)
title('Mach-Varying Reduced-Order Stability Estimate at Max Q')
grid on

% 12. Projected area
A_projected = A_ref*abs(cos(alpha_AoA)) ...
    +L_body*D_body*abs(sin(alpha_AoA));
A_fins_projected = n_fins*(fin_span*fin_thickness) ...
    *abs(cos(alpha_AoA)) + 2*A_fin*abs(sin(alpha_AoA));
figure
plot(t,A_projected,t,A_projected+A_fins_projected)
legend('Body only','Body + four fins')
xlabel('Time (s)')
ylabel('Projected Area (m^2)')
title('Effect of Fins on Projected Area')
xlim([0 5])
grid on

% 13. Static margin
figure
plot(t,static_margin,'LineWidth',1.3)
hold on
plot(t,(x_CP-x_CoM)/D_body,'--')
yline(0,':')
legend('Mach-dependent','Original constant-C_P','Zero margin')
xlabel('Time (s)')
ylabel('Static Margin (calibers)')
title('Rocket Static Stability')
xlim([0 15])
grid on

% 14. Passive aerodynamic pitching moment
figure
plot(t,torque_total-torque_fin)
xlabel('Time (s)')
ylabel('Pitch Moment (N m)')
title('Passive Aerodynamic Restoring Moment')
xlim([0 5])
grid on

% 15. All-moving fin hinge torque estimate
F_fin_each = q_dyn .* ( ...
    A_ref*(CN_fins_hist/n_fins).*alpha_AoA ...
    - A_fin*CL_max.*tanh(CL_alpha_hist.*delta/CL_max));
torque_hinge = F_fin_each*hinge_arm;
torque_hinge(V<=0) = 0;
figure
plot(t,torque_hinge)
xlabel('Time (s)')
ylabel('Hinge Torque per Fin (N m)')
title('Aerodynamic Hinge Torque - All-Moving Fin')
xlim([0 15])
grid on

% 16. Mach during ascent
figure
plot(t,Mach,'LineWidth',1.3)
yline(1,'--','Mach 1')
xlabel('Time (s)')
ylabel('Mach number')
title('Rocket Mach Number')
xlim([0 t(idx_apogee)])
grid on

% 17. Empirical Cd reference, not the geometry-based flight Cd
M_plot = linspace(0,2,500);
Cd_plot = Cd0*interp1(M_table,Cd_ratio,M_plot,'pchip');
figure
plot(M_plot,Cd_plot,'LineWidth',1.5)
xline(1,'--','Mach 1')
xlabel('Mach number')
ylabel('Drag coefficient Cd')
title('Mach-Dependent Rocket Drag Coefficient (Reference Curve)')
grid on

% 18. Conical nose pressure/wave-drag coefficient
M_eval = linspace(0,2,401);
Cd_nose_pressure = zeros(size(M_eval));
sub = M_eval <= 1;
Cd_nose_pressure(sub) = Cd_n0 ...
    +(Cd_n1-Cd_n0)*M_eval(sub).^nose_exponent;
trans = M_eval > 1 & M_eval < M_n2;
u = (M_eval(trans)-1)/(M_n2-1);
dM = M_n2-1;
Cd_nose_pressure(trans) = ...
    (2*u.^3-3*u.^2+1)*Cd_n1 ...
    +(u.^3-2*u.^2+u)*dM*dCd_n1 ...
    +(-2*u.^3+3*u.^2)*Cd_n2 ...
    +(u.^3-u.^2)*dM*dCd_n2;
sup = M_eval >= M_n2;
Cd_nose_pressure(sup) = 2.1*s_eps^2 ...
    +0.5*s_eps./sqrt(M_eval(sup).^2-1);
figure
plot(M_eval,Cd_nose_pressure,'LineWidth',1.5)
xline(1,'--','Mach 1')
xlabel('Mach number')
ylabel('Nose pressure drag coefficient')
title('Conical Nose Pressure Drag vs Mach')
grid on

% 19. Component drag curves (diagnostic decomposition)
L_cyl = L_body-L_nose;
S_wet_cyl = pi*D_body*L_cyl;
S_wet_nose = pi*R_body*hypot(L_nose,R_body);
S_wet_fins = 2*n_fins*A_fin;
L_parts = [L_cyl L_nose c_mean];
S_parts = [S_wet_cyl S_wet_nose S_wet_fins];
K_parts = [K_body K_body K_fin];
T_drag = max(216.65,288.15-0.0065*h);
mu_drag = 1.716e-5*(T_drag/273.15).^1.5 ...
    .*(273.15+110.4)./(T_drag+110.4);
D_friction_series = zeros(size(t));
for p = 1:3
    Re = max(rho.*air_speed*L_parts(p)./mu_drag,1e4);
    Cf_smooth = 1./(1.5*log(Re)-5.6).^2;
    Cf_rough = 0.032*(roughness_selected/L_parts(p))^0.2;
    smooth_factor = zeros(size(t));
    rough_factor = zeros(size(t));
    sub = Mach < 1;
    smooth_factor(sub) = 1-0.1*Mach(sub).^2;
    rough_factor(sub) = smooth_factor(sub);
    smooth_factor(~sub) = 1./(1+0.15*Mach(~sub).^2).^0.58;
    rough_factor(~sub) = 1./(1+0.18*Mach(~sub).^2);
    Cf = max(Cf_smooth.*smooth_factor,Cf_rough.*rough_factor);
    D_friction_series = D_friction_series ...
        +q_dyn.*Cf*S_parts(p)*K_parts(p);
end
Cd_nose_series = interp1(M_eval,Cd_nose_pressure, ...
    min(max(Mach,0),2),'pchip');
D_nose_series = q_dyn.*A_ref.*Cd_nose_series;
Cd_LE_series = zeros(size(Mach));
a = Mach < 0.9;
b = Mach >= 0.9 & Mach < 1;
c = Mach >= 1;
Cd_LE_series(a) = (1-Mach(a).^2).^(-0.417)-1;
Cd_LE_series(b) = 1-1.785*(Mach(b)-0.9);
Cd_LE_series(c) = 1.214-0.502./Mach(c).^2 ...
    +0.1095./Mach(c).^4;
D_fin_series = q_dyn.*n_fins*A_fin_edge ...
    .*Cd_LE_series*cos(sweep_LE)^2;
Cd_base_series = zeros(size(Mach));
sub = Mach <= 1;
Cd_base_series(sub) = 0.12+0.13*Mach(sub).^2;
Cd_base_series(~sub) = 0.25./Mach(~sub);
D_base_series = q_dyn.*Cd_base_series ...
    .*(A_base-A_motor*(m_fuel>0));
D_model = D_friction_series+D_nose_series+D_fin_series+D_base_series;
idx_plot = t <= burn_time+5;
figure
plot(t(idx_plot),D_model(idx_plot),'LineWidth',1.5)
hold on
plot(t(idx_plot),D_friction_series(idx_plot))
plot(t(idx_plot),D_nose_series(idx_plot))
plot(t(idx_plot),D_fin_series(idx_plot))
plot(t(idx_plot),D_base_series(idx_plot))
legend('Total','Friction','Nose','Fins','Base')
xlabel('Time (s)')
ylabel('Drag force (N)')
title('Geometry-Based Drag Throughout Ascent')
grid on

% 20. Wind gust attitude
figure
plot(t,rad2deg(theta),t,rad2deg(yaw),'LineWidth',1.3)
legend('Pitch','Yaw')
xlabel('Time (s)')
ylabel('Angle relative to vertical (deg)')
title('Rocket Attitude Under Wind Gust')
xlim([5 15])
grid on

% 21. Structural sensor measurement: raw, not feedback after compensation
figure
plot(t,rad2deg(theta),t,rad2deg(yaw),'LineWidth',1.3)
hold on
plot(t,rad2deg(theta+q_flex),'--',t,rad2deg(yaw+q_flex_yaw),'--','LineWidth',1.3)
legend('Pitch true','Yaw true','Pitch sensor raw','Yaw sensor raw')
xlabel('Time (s)')
ylabel('Pitch angle (deg)')
title('Structural Bending and Attitude Measurement Error')
xlim([5 15])
grid on

% 22. Torque-estimator verification
figure
plot(t(gust_idx),torque_fin(gust_idx),t(gust_idx),torque_fin_yaw(gust_idx),'LineWidth',1.3)
hold on
plot(t(gust_idx),torque_est_history(gust_idx),'--', ...
    t(gust_idx),torque_est_yaw_history(gust_idx),'--','LineWidth',1.3)
legend('Pitch actual','Yaw actual','Pitch estimated','Yaw estimated')
xlabel('Time (s)')
ylabel('Torque (N m)')
title('Fin Torque Estimation During Gust')
grid on


%% Mach-dependent aerodynamic stability plots (5 new figures)
Ma_sweep = linspace(0,2,401);
CN_nose_sweep = CN_nose*ppval(pp_nose,Ma_sweep);
CN_fin_sweep = CN_fins*ppval(pp_fin,Ma_sweep);
CL_control_sweep = CL_alpha*ppval(pp_ctrl,Ma_sweep);
x_fin_sweep = x_MAC_LE+c_mean*ppval(pp_ac,Ma_sweep);
x_CP_sweep = (CN_nose_sweep*x_CP_nose+CN_fin_sweep.*x_fin_sweep) ...
    ./(CN_nose_sweep+CN_fin_sweep);

% 23. Mach-dependent aerodynamic derivatives
figure
plot(Ma_sweep,CN_nose_sweep,'LineWidth',1.3)
hold on
plot(Ma_sweep,CN_fin_sweep,'LineWidth',1.3)
plot(Ma_sweep,CL_control_sweep,'LineWidth',1.3)
xline(1,':','Mach 1')
legend('Nose C_{N\alpha}','Four-fin C_{N\alpha}', ...
    'Control-fin C_{L\alpha}','Location','best')
xlabel('Mach number')
ylabel('Slope (per radian)')
title('Assumed Mach-Dependent Normal-Force Slopes')
grid on

% 24. Aerodynamic center shift with Mach
figure
plot(Ma_sweep,x_CP_sweep,'LineWidth',1.4)
hold on
plot(Ma_sweep,x_fin_sweep,'--','LineWidth',1.2)
yline(x_CP,':','Low-Mach CP')
xline(1,':','Mach 1')
legend('Overall CP','Fin AC','Original constant CP','Location','best')
xlabel('Mach number')
ylabel('Distance from nose (m)')
title('Mach-Dependent Aerodynamic Center (Provisional)')
grid on

% 25. Restoring stiffness along actual flight
figure
plot(t,Kpass_hist,'LineWidth',1.3)
hold on
yline(0,':')
xlabel('Time (s)')
ylabel('Linear restoring stiffness (N m/rad)')
title('Estimated Aerodynamic Pitch Restoring Stiffness')
xlim([0 15])
grid on

% 26. Control authority at a fixed dynamic pressure, isolating Mach effect
figure
Ka_sweep = maxQ*S_fin.*CL_control_sweep ...
    .*(x_fin_sweep-x_CoM(idx_maxQ));
plot(Ma_sweep,Ka_sweep,'LineWidth',1.4)
hold on
xline(1,':','Mach 1')
xlabel('Mach number')
ylabel('dM_{control}/d\delta (N m/rad)')
title('Small-Angle Fin Control Authority vs Mach (Fixed Max Q)')
grid on

% 27. Approximate closed-loop poles sampled along powered flight/coast
sample_times = 2:0.5:15;
max_real_pole = NaN(size(sample_times));
for k = 1:numel(sample_times)
    [~,j] = min(abs(t-sample_times(k)));
    if air_speed(j)<20 || V(j)<=1
        continue
    end
    Gr_j = 1/(I_total(j)*s^2 + Kq_hist(j)*s+Kpass_hist(j));
    Loop_j = C*Gs*Ka_hist(j)*(Gr_j+Gflex_residual)*N_cont;
    CL_j = feedback(Loop_j,1);
    p_j = pole(CL_j);
    max_real_pole(k) = max(real(p_j));
end
figure
plot(sample_times,max_real_pole,'o-','LineWidth',1.2)
hold on
yline(0,':','Stability boundary')
xlabel('Flight time (s)')
ylabel('Largest real part of pole (1/s)')
title('Legacy Pitch-Only Reduced-Order Poles (NOT 3D Stability)')
grid on

fprintf('\nLEGACY PITCH-ONLY REDUCED-ORDER ESTIMATE (not 3D validation)\n');
fprintf('Mach at max Q: %.3f\n',Mach(idx_maxQ));
fprintf('CP at max Q: %.3f m (low-Mach reference %.3f m)\n', ...
    x_CP_hist(idx_maxQ),x_CP);
fprintf('Static margin at max Q: %.2f calibers\n',static_margin(idx_maxQ));
fprintf('Restoring stiffness at max Q: %.1f N m/rad\n',Kpass_hist(idx_maxQ));
fprintf('Control authority at max Q: %.1f N m/rad\n',Ka_hist(idx_maxQ));
idx_valid = isfinite(max_real_pole);
if any(idx_valid)
    fprintf('Largest reduced-order pole real part: %+.3f 1/s\n', ...
        max(max_real_pole(idx_valid)));
end
fprintf('These transonic curves are illustrative, NOT validated rocket data.\n');

%% Full-model linear stability analysis status
% The v7 16-state 2D Jacobian is deliberately not reused for three axes.
% Rebuild around quaternion and all fin servo/filter and roll states
% before interpreting new eigenvalues. Existing Mach-coefficient graphs are
% retained as static aerodynamic diagnostics only.

%% Ballistic descent plots (new; existing figures preserved)

% Horizontal/vertical plane trajectory, including descent to touchdown.
figure
plot(x_pos(1:idx_apogee),h(1:idx_apogee),'LineWidth',1.5)
hold on
plot(x_pos(idx_apogee:end),h(idx_apogee:end),'LineWidth',1.5)
plot(x_pos(idx_apogee),h(idx_apogee),'o','MarkerSize',7)
if isfinite(impact_idx)
    plot(x_pos(end),h(end),'x','MarkerSize',9)
    legend('Ascent','Ballistic descent','Apogee','Touchdown','Location','best')
else
    legend('Ascent','Ballistic descent','Apogee','Location','best')
end
xlabel('Downrange (m)')
ylabel('Altitude (m)')
title('Full 2D Trajectory: Ascent, Apogee and Unpowered Descent')
grid on

% Descent velocity: vertical velocity must become negative.
figure
plot(t(idx_apogee:end),V(idx_apogee:end),'LineWidth',1.3)
hold on
plot(t(idx_apogee:end),vx(idx_apogee:end),'LineWidth',1.3)
plot(t(idx_apogee:end),vy(idx_apogee:end),'LineWidth',1.3)
yline(0,':')
legend('Vertical velocity','X velocity','Y velocity','Zero','Location','best')
xlabel('Time (s)')
ylabel('Ground-relative velocity (m/s)')
title('Ballistic Descent Velocities')
grid on

% Passive dynamics, locked fins throughout descent.
figure
subplot(2,1,1)
plot(t(idx_apogee:end),rad2deg(theta(idx_apogee:end)), ...
    t(idx_apogee:end),rad2deg(yaw(idx_apogee:end)))
legend('Pitch','Yaw')
ylabel('Body attitude (deg)')
title('Post-Apogee Passive Attitude and Fin Lock')
grid on
subplot(2,1,2)
plot(t(idx_apogee:end),rad2deg(delta(idx_apogee:end)), ...
    t(idx_apogee:end),rad2deg(delta_yaw(idx_apogee:end)))
legend('Pitch fins','Yaw fins')
ylabel('Fin deflection (deg)')
xlabel('Time (s)')
grid on

% Dynamic pressure from apogee to touchdown, with no parachute.
figure
plot(t(idx_apogee:end),q_dyn(idx_apogee:end)/1000,'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Dynamic pressure (kPa)')
title('Post-Apogee Ballistic Dynamic Pressure')
grid on

% Passive aerodynamic pitching moment with no fin actuation.
figure
plot(t(idx_apogee:end),torque_total(idx_apogee:end), ...
    t(idx_apogee:end),torque_total_yaw(idx_apogee:end),'LineWidth',1.3)
legend('Pitch moment','Yaw moment')
yline(0,':')
xlabel('Time (s)')
ylabel('Passive pitch moment (N m)')
title('Post-Apogee Static-Fin and Body Pitching Moment')
grid on

%% Descent aerodynamic diagnostics — all five use simulated descent states
figure('Name','Post-Apogee Aerodynamic Diagnostics', ...
    'Position',[100 80 1150 780])
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(t(d_idx),Mach(d_idx),'LineWidth',1.2)
ylabel('Mach')
title('Air-relative Mach')
grid on

nexttile
plot(t(d_idx),Cd_effective(d_idx),'LineWidth',1.2)
ylabel('Effective C_D')
title('Aerodynamic Drag Coefficient (q >= 100 Pa)')
grid on

nexttile
plot(t(d_idx),flow_incidence_deg(d_idx),'LineWidth',1.2)
ylabel('Incidence (deg)')
title('0° Nose-first, 180° Tail-first')
ylim([0 180])
grid on

nexttile
plot(t(d_idx),torque_total(d_idx),t(d_idx),torque_total_yaw(d_idx), ...
     t(d_idx),torque_roll(d_idx),'LineWidth',1.2)
legend('Pitch','Yaw','Roll')
ylabel('Moment (N m)')
title('Passive Aerodynamic Pitch/Yaw/Roll Moments')
grid on

nexttile
plot(t(d_idx),rad2deg(omega(d_idx)), ...
    t(d_idx),rad2deg(omega_yaw(d_idx)), ...
    t(d_idx),rad2deg(omega_roll(d_idx)),'LineWidth',1.2)
legend('Pitch','Yaw','Roll')
ylabel('Rate (deg/s)')
title('Rigid-body Angular Velocity')
grid on

nexttile
plot(t(d_idx),accel_proper_g(d_idx),'LineWidth',1.2)
ylabel('Proper acceleration (g)')
title('Aerodynamic Load During Descent')
grid on

xlabel(tl,'Flight time (s)')
sgtitle('Ballistic Descent: Mach, Drag, Incidence, Moments and Rotation')


%% 3D flight trajectory and ground track (yaw +5 deg)
figure('Name','3D Rocket Trajectory')
plot3(x_pos(1:idx_apogee),y_pos(1:idx_apogee), ...
      h(1:idx_apogee),'LineWidth',1.6)
hold on
plot3(x_pos(idx_apogee:end),y_pos(idx_apogee:end), ...
      h(idx_apogee:end),'LineWidth',1.6)
plot3(x_pos(idx_apogee),y_pos(idx_apogee),h(idx_apogee), ...
      'o','MarkerSize',8,'LineWidth',1.5)
plot3(x_pos(end),y_pos(end),h(end),'x','MarkerSize',9,'LineWidth',1.5)
legend('Ascent','Ballistic descent','Apogee','Touchdown','Location','best')
xlabel('X downrange (m)');ylabel('Y downrange (m)');zlabel('Altitude (m)')
title('Full 3D Trajectory - Pitch +5, Yaw +5, Roll +180 deg')
grid on; axis tight; view(35,22)

figure('Name','Horizontal X-Y Ground Track')
plot(x_pos(1:idx_apogee),y_pos(1:idx_apogee),'LineWidth',1.4)
hold on
plot(x_pos(idx_apogee:end),y_pos(idx_apogee:end),'LineWidth',1.4)
plot(x_pos(idx_apogee),y_pos(idx_apogee),'o','MarkerSize',7)
plot(x_pos(end),y_pos(end),'x','MarkerSize',8)
axis equal; grid on
legend('Ascent','Descent','Apogee','Touchdown','Location','best')
xlabel('X downrange (m)');ylabel('Y downrange (m)')
title('Horizontal Landing Footprint')

figure('Name','Pitch Yaw Roll Control Comparison')
tiledlayout(4,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,rad2deg(theta),t,rad2deg(yaw),t,rad2deg(roll),'LineWidth',1.2)
legend('Pitch','Yaw','Roll');ylabel('Attitude (deg)');grid on
xlim([0 min(15,t(end))]);title('Pitch +5, Yaw +5, Roll +180 Degree Offsets')
nexttile
plot(t,rad2deg(delta),t,rad2deg(delta_yaw),t,rad2deg(delta_roll),'LineWidth',1.2)
legend('Pitch mix','Yaw mix','Roll mix');ylabel('Deflection (deg)');grid on
xlim([0 min(15,t(end))])
nexttile
plot(t,rad2deg(omega),t,rad2deg(omega_yaw),t,rad2deg(omega_roll),'LineWidth',1.2)
legend('Pitch','Yaw','Roll');ylabel('Rate (deg/s)');grid on
xlim([0 min(15,t(end))])
nexttile
plot(t,torque_total,t,torque_total_yaw,t,torque_roll,'LineWidth',1.2)
legend('Pitch','Yaw','Roll');ylabel('Aero moment (N m)')
xlabel('Time (s)');grid on;xlim([0 min(15,t(end))])

fprintf('\n3D PITCH + YAW + ROLL CHECKPOINT\n')
fprintf('Initial pitch/yaw: %.3f / %.3f deg\n', ...
    rad2deg(theta(1)),rad2deg(yaw(1)))
fprintf('Burnout pitch/yaw: %.3f / %.3f deg\n', ...
    rad2deg(theta(burn_idx(end))),rad2deg(yaw(burn_idx(end))))
fprintf('Apogee position [X Y Z]: %.2f %.2f %.2f m\n', ...
    x_pos(idx_apogee),y_pos(idx_apogee),h(idx_apogee))
fprintf('Touchdown [X Y Z]: %.2f %.2f %.2f m\n', ...
    x_pos(end),y_pos(end),h(end))
fprintf('Maximum yaw fin deflection FROM LOCK: %.6f deg\n', ...
    max(abs(rad2deg(delta_yaw(apogee_lock_idx:end)))))
fprintf('Initial roll: %.3f deg, wrapped burnout roll: %.3f deg\n', ...
    rad2deg(roll_initial),rad2deg(roll(burn_idx(end))))
fprintf('Roll tie preference: %+.0f body-X direction\n',roll_tie_direction)
fprintf('Roll error at burnout: %.3f deg\n', ...
    rad2deg(atan2(sin(roll(burn_idx(end))),cos(roll(burn_idx(end))))))
fprintf('Peak roll rate in powered flight: %.2f deg/s\n', ...
    max(abs(rad2deg(omega_roll(burn_idx)))))
fprintf('Peak roll moment in powered flight: %.2f N m\n', ...
    max(abs(torque_roll(burn_idx))))
fprintf('Peak absolute per-fin deflection: %.3f deg (limit %.1f)\n', ...
    max(abs(rad2deg(fin_deflection(:)))),rad2deg(delta_max))
fprintf('Mixer saturation fraction (powered): %.2f%%\n', ...
    100*mean(fin_mixer_scale(burn_idx)<1))
fprintf('Maximum fin deflection FROM LOCK: %.6f deg\n', ...
    max(abs(rad2deg(fin_deflection(:,apogee_lock_idx:end))),[],'all'))
fprintf('V7 pitch-only frozen Jacobian intentionally omitted in 3D.\n')


%% Roll response and per-fin actuation (new)
figure('Name','Roll Correction and Direction Choice')
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,rad2deg(roll),'LineWidth',1.25)
hold on
yline(0,':','Target');xline(t(idx_apogee),':','Apogee')
ylabel('Wrapped roll (deg)')
title('Roll: 180 deg Initial Offset, Shortest-Angle Correction')
xlim([0 min(15,t(end))]);grid on
nexttile
plot(t,rad2deg(omega_roll),'LineWidth',1.25)
ylabel('Body X roll rate (deg/s)');xlim([0 min(15,t(end))]);grid on
nexttile
plot(t,torque_roll,t,roll_requested_moment,'LineWidth',1.2)
legend('Actual aerodynamic roll moment','Requested roll moment')
ylabel('Roll torque (N m)');xlabel('Flight time (s)')
xlim([0 min(15,t(end))]);grid on

figure('Name','Four Fin Mixed Actuation')
plot(t,rad2deg(fin_deflection.'),'LineWidth',1.1)
hold on
yline(rad2deg(delta_max),'--','+5 deg');
yline(-rad2deg(delta_max),'--','-5 deg')
xline(t(idx_apogee),':','Lock at apogee')
legend('Fin 1','Fin 2','Fin 3','Fin 4','Location','best')
xlabel('Flight time (s)');ylabel('Actual fin angle (deg)')
title('Four Independent Servos: Mixed Pitch + Yaw + Roll')
xlim([0 min(15,t(end))]);grid on

figure('Name','Full-Flight Roll Spin and Yaw Pitch')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,rad2deg(theta),t,rad2deg(yaw),t,rad2deg(roll),'LineWidth',1.1)
legend('Pitch','Yaw','Roll (wrapped)','Location','best')
ylabel('Orientation (deg)');title('3-Axis Body Orientation')
grid on
nexttile
plot(t,rad2deg(roll_spin-roll_spin(1)),'LineWidth',1.1)
hold on
xline(t(idx_apogee),':','Apogee')
ylabel('Integrated axial spin (deg)');xlabel('Flight time (s)')
title('Accumulated Body-X Rotation (Not Euler Roll)');grid on

fprintf('\n3-AXIS MIXER / ROLL DIAGNOSTIC\n')
fprintf('Initial pitch / yaw / roll: %.1f / %.1f / %.1f deg\n', ...
    rad2deg(theta_initial),rad2deg(yaw_initial),rad2deg(roll_initial))
fprintf('Burnout pitch / yaw / roll: %.3f / %.3f / %.3f deg\n', ...
    rad2deg(theta(burn_idx(end))),rad2deg(yaw(burn_idx(end))), ...
    rad2deg(roll(burn_idx(end))))
fprintf('Max yaw/pitch angles during first 10s: %.3f / %.3f deg\n', ...
    max(abs(rad2deg(yaw(burn_idx)))),max(abs(rad2deg(theta(burn_idx)))))
fprintf('Axial spin during powered flight: %.2f deg\n', ...
    rad2deg(roll_spin(burn_idx(end))-roll_spin(1)))
fprintf('Roll inertia at liftoff/burnout: %.5f / %.5f kg m^2\n', ...
    I_roll(1),I_roll(burn_idx(end)))
fprintf('Any mixer individual fin violations: %d\n', ...
    any(abs(fin_deflection(:))>delta_max+1e-10))
fprintf('Maximum actual per-fin deflection FROM APOGEE LOCK: %.9f deg\n', ...
    max(abs(rad2deg(fin_deflection(:,apogee_lock_idx:end))),[],'all'))

%% EI-based bending versus the original 4 Hz sensor-side oscillator
% The main flight/controls above are identical to v11. This is a one-way
% structural experiment, NOT a closed-loop aeroelastic validation.
beam_Q_pitch = beam_Q_pitch(1:n);
beam_Q_yaw = beam_Q_yaw(1:n);
beam_modal_mass = (m_dry/L_body)*beam_mass_shape ...
    +m_fuel*beam_phi_fuel^2;
beam_frequency = sqrt(beam_K./beam_modal_mass)/(2*pi);

% A separate smooth 4 Hz shaker force, applied at the fin station.
% A 4 Hz INPUT need not make the BODY's natural frequency 4 Hz.
beam_4hz_force = zeros(1,n);
beam_force_mask = t>=beam_4Hz_start & ...
    t<=beam_4Hz_start+beam_4Hz_duration;
beam_phase = (t(beam_force_mask)-beam_4Hz_start)/beam_4Hz_duration;
beam_4hz_force(beam_force_mask) = beam_force_4Hz_N ...
    .*sin(pi*beam_phase).^2 ...
    .*sin(2*pi*beam_4Hz_frequency* ...
        (t(beam_force_mask)-beam_4Hz_start));
beam_Q4_pitch = beam_phi_fin0*beam_4hz_force;
beam_Q4_yaw = beam_yaw_shaker_fraction*beam_Q4_pitch;

% Four parallel modal channels: pitch/aero-only, pitch/aero+4Hz,
% yaw/aero-only, yaw/aero+4Hz. Coordinates q are physical meters.
beam_q = zeros(4,n);
beam_qd = zeros(4,n);
for j = 1:n-1
    Mm = beam_modal_mass(j);
    Cm = 2*zeta_body*sqrt(beam_K*Mm);
    generalized_F = [beam_Q_pitch(j); ...
        beam_Q_pitch(j)+beam_Q4_pitch(j); ...
        beam_Q_yaw(j); ...
        beam_Q_yaw(j)+beam_Q4_yaw(j)];
    beam_accel = (generalized_F-Cm*beam_qd(:,j) ...
        -beam_K*beam_q(:,j))/Mm;
    beam_qd(:,j+1) = beam_qd(:,j)+dt*beam_accel;
    beam_q(:,j+1) = beam_q(:,j)+dt*beam_qd(:,j+1);
end

% Local sensor angular error = derivative of mode shape * displacement.
beam_sensor_pitch_aero = beam_slope_sensor*beam_q(1,:);
beam_sensor_pitch_shaken = beam_slope_sensor*beam_q(2,:);
beam_sensor_yaw_aero = beam_slope_sensor*beam_q(3,:);
beam_sensor_yaw_shaken = beam_slope_sensor*beam_q(4,:);
% q_flex in V11 remains the SEPARATE 4 Hz sensor-side/mount DOF.
beam_sensor_pitch_combined = beam_sensor_pitch_shaken+q_flex;
beam_sensor_yaw_combined = beam_sensor_yaw_shaken+q_flex_yaw;

beam_active_idx = t<=burn_time;
beam_test_idx = t>=beam_4Hz_start & ...
    t<=beam_4Hz_start+beam_4Hz_duration;

fprintf('\nPHYSICAL BEAM + SEPARATE 4 Hz TEST (ONE-WAY)\n')
fprintf('Assumed EI: %.0f N m^2\n',EI_body)
fprintf('Body first-mode frequency liftoff/burnout: %.2f / %.2f Hz\n', ...
    beam_f_liftoff,beam_f_burnout)
fprintf('Generalized modal stiffness: %.1f N/m\n',beam_K)
fprintf('Modal mass liftoff/burnout: %.3f / %.3f kg\n', ...
    beam_M_liftoff,beam_M_burnout)
fprintf('Maximum physical nose bending (aero + 4 Hz input): %.4f mm\n', ...
    1000*max(abs(beam_phi(1)*beam_q(2,beam_active_idx))))
fprintf('Maximum physical tail bending (aero + 4 Hz input): %.4f mm\n', ...
    1000*max(abs(beam_phi_tail*beam_q(2,beam_active_idx))))
fprintf('Peak physical sensor rotation (pitch): %.5f deg\n', ...
    max(abs(rad2deg(beam_sensor_pitch_shaken(beam_active_idx)))))
fprintf('Peak separate 4 Hz sensor-side rotation (pitch): %.3f deg\n', ...
    max(abs(rad2deg(q_flex(beam_active_idx)))))
fprintf('Additional body sensor rotation from applied 4 Hz force: %.6f deg\n', ...
    max(abs(rad2deg(beam_sensor_pitch_shaken-beam_sensor_pitch_aero))))
fprintf('Applied 4 Hz peak shaker force: %.1f N (assumed experiment)\n', ...
    beam_force_4Hz_N)
fprintf('This test DOES NOT apply deformation back into aerodynamic forces.\n')

figure('Name','Beam Mode Shape and Physical Frequency')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(beam_x,beam_phi,'LineWidth',1.3)
hold on
xline(x_sensor_beam,':','Sensor'); xline(x_fin_CP,':','Fin station')
ylabel('Normalized first mode');xlabel('Position aft of nose (m)')
title('Free-Free First Elastic Beam Mode (Physical Shape)');grid on
nexttile
plot(t,beam_frequency,'LineWidth',1.2)
hold on
yline(beam_4Hz_frequency,'--','4 Hz forcing')
xline(burn_time,':','Burnout')
xlabel('Flight time (s)');ylabel('Frequency (Hz)')
title('Body Natural Frequency Calculated from EI and Mass');grid on

figure('Name','Beam Aerodynamic Modal Loads')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,beam_Q_pitch,t,beam_Q_yaw,'LineWidth',1.05)
legend('Body-Y / pitch plane','Body-Z / yaw plane','Location','best')
ylabel('Generalized aero force (N)')
title('Actual V11 Nose, Body and Fin Loads Projected onto Beam');grid on
xlim([0 min(15,t(end))])
nexttile
plot(t,beam_Q4_pitch,t,beam_Q4_yaw,'LineWidth',1.05)
legend('Additional pitch-plane 4 Hz force','Additional yaw-plane 4 Hz force')
ylabel('Generalized force (N)');xlabel('Flight time (s)')
grid on;xlim([0 min(15,t(end))])

figure('Name','Physical Beam Bending in mm')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,1000*beam_phi_tail*beam_q(2,:), ...
     t,1000*beam_phi_tail*beam_q(4,:),'LineWidth',1.2)
legend('Tail transverse displacement: pitch plane', ...
       'Tail transverse displacement: yaw plane')
ylabel('Physical deflection (mm)')
title('EI-Based First-Mode Deformation, Including 4 Hz Force')
xlim([0 min(15,t(end))]);grid on
nexttile
plot(t,rad2deg(beam_sensor_pitch_shaken),'LineWidth',1.2)
hold on
plot(t,rad2deg(beam_sensor_yaw_shaken),'LineWidth',1.2)
legend('Physical pitch-plane slope','Physical yaw-plane slope')
ylabel('Sensor-location slope (deg)');xlabel('Flight time (s)')
xlim([0 min(15,t(end))]);grid on

figure('Name','Separate 4Hz Mount vs Physical Beam')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,rad2deg(beam_sensor_pitch_shaken),'LineWidth',1.25)
hold on
plot(t,rad2deg(q_flex),'LineWidth',1.25)
plot(t,rad2deg(beam_sensor_pitch_combined),'--','LineWidth',1.1)
legend('Physical beam slope','Separate 4 Hz legacy mount', ...
       'Combined pitch sensor perturbation','Location','best')
ylabel('Pitch sensor error (deg)');xlim([0 min(15,t(end))]);grid on
nexttile
plot(t,rad2deg(beam_sensor_yaw_shaken),'LineWidth',1.25)
hold on
plot(t,rad2deg(q_flex_yaw),'LineWidth',1.25)
plot(t,rad2deg(beam_sensor_yaw_combined),'--','LineWidth',1.1)
legend('Physical beam slope','Separate 4 Hz legacy mount', ...
       'Combined yaw sensor perturbation','Location','best')
ylabel('Yaw sensor error (deg)');xlabel('Flight time (s)')
xlim([0 min(15,t(end))]);grid on

figure('Name','4Hz Shaker Effect on Stiff Beam')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')
nexttile
plot(t,1000*beam_phi_tail*beam_q(1,:), ...
     t,1000*beam_phi_tail*beam_q(2,:),'LineWidth',1.1)
legend('Actual V11 aero loads only','Aero loads + 4 Hz test force')
ylabel('Pitch-plane tail displacement (mm)')
xlim([beam_4Hz_start-1 beam_4Hz_start+beam_4Hz_duration+1]);grid on
nexttile
plot(t,rad2deg(beam_sensor_pitch_shaken-beam_sensor_pitch_aero), ...
     'LineWidth',1.1)
ylabel('Extra 4 Hz sensor slope (deg)')
xlabel('Flight time (s)')
xlim([beam_4Hz_start-1 beam_4Hz_start+beam_4Hz_duration+1]);grid on

function R = rocketQuatToRotm(q)
% Hamilton quaternion [w x y z]; body-to-world DCM.
q=q(:).'/norm(q);
w=q(1); x=q(2); y=q(3); z=q(4);
R=[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y); ...
   2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x); ...
   2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)];
end

function q = rocketQuatMultiply(a,b)
% Hamilton product; a,b are row quaternions.
q=[a(1)*b(1)-dot(a(2:4),b(2:4)), ...
   a(1)*b(2)+b(1)*a(2)+a(3)*b(4)-a(4)*b(3), ...
   a(1)*b(3)+b(1)*a(3)+a(4)*b(2)-a(2)*b(4), ...
   a(1)*b(4)+b(1)*a(4)+a(2)*b(3)-a(3)*b(2)];
end

function q = rocketRotmToQuat(R)
% Construct the body-to-world initial quaternion (no extra toolbox).
t=trace(R);
if t>0
    S=sqrt(t+1)*2;
    q=[S/4,(R(3,2)-R(2,3))/S, ...
       (R(1,3)-R(3,1))/S,(R(2,1)-R(1,2))/S];
else
    [~,k]=max(diag(R));
    if k==1
        S=sqrt(1+R(1,1)-R(2,2)-R(3,3))*2;
        q=[(R(3,2)-R(2,3))/S,S/4, ...
           (R(1,2)+R(2,1))/S,(R(1,3)+R(3,1))/S];
    elseif k==2
        S=sqrt(1+R(2,2)-R(1,1)-R(3,3))*2;
        q=[(R(1,3)-R(3,1))/S,(R(1,2)+R(2,1))/S, ...
           S/4,(R(2,3)+R(3,2))/S];
    else
        S=sqrt(1+R(3,3)-R(1,1)-R(2,2))*2;
        q=[(R(2,1)-R(1,2))/S,(R(1,3)+R(3,1))/S, ...
           (R(2,3)+R(3,2))/S,S/4];
    end
end
q=q/norm(q);
end

function R = rocketZeroRollFrame(nose_world)
% Body X points along nose. Body Y is the projected world X direction.
% Body Z completes the right-handed triad. The reference is unique except
% when nose is nearly parallel to world X, handled with a world-Y fallback.
    e1 = nose_world(:)/norm(nose_world);
    ref=[1;0;0];
    e2=ref-dot(ref,e1)*e1;
    if norm(e2)<0.05
        ref=[0;1;0];
        e2=ref-dot(ref,e1)*e1;
    end
    e2=e2/norm(e2);
    e3=cross(e1,e2); e3=e3/norm(e3);
    R=[e1,e2,e3];
end
