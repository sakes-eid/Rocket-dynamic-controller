clear; clc; close all

%% Parameters

% Integration
dt = 0.00025;
t = 0:dt:120;
n = length(t);

% Horizontal wind gust
wind_peak = -15;       % m/s, horizontal gust
gust_start = 7;     % s
gust_duration = 4; % s

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
theta_initial = deg2rad(5);

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
enable_mach_aero = false;

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

drag_force = zeros(1,n);
wind_history = zeros(1,n);

% Flight
V = zeros(1,n);
h = zeros(1,n);

% Horizontal motion
vx = zeros(1,n);
x_pos = zeros(1,n);

rho = zeros(1,n);
q_dyn = zeros(1,n);

% Pitch
theta = zeros(1,n);
omega = zeros(1,n);
delta = zeros(1,n);

% Flexibility
q_flex = zeros(1,n);
qdot_flex = zeros(1,n);

% Aerodynamic moments
torque_fin = zeros(1,n);
torque_total = zeros(1,n);

theta(1) = theta_initial;

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

%% Coupled flight and control simulation

for i = 1:n-1
    % Atmosphere and gust

    rho(i) = rho0*exp(-h(i)/H);

    % Smooth horizontal gust
    wind_x = 0;

    if t(i) >= gust_start && ...
            t(i) <= gust_start+gust_duration

        phase = (t(i)-gust_start)/gust_duration;

        wind_x = 0.5*wind_peak ...
            *(1-cos(2*pi*phase));
    end

    wind_history(i) = wind_x;

    % Air-relative horizontal velocity
    v_air_x = vx(i)-wind_x;

    speed_air = hypot(v_air_x,V(i));

    q_dyn(i) = 0.5*rho(i)*speed_air^2;

    % Air-relative flight path and AoA

    gamma_i = atan2(v_air_x,V(i));

    alpha_i = atan2( ...
        sin(theta(i)-gamma_i), ...
        cos(theta(i)-gamma_i));

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

% Air-relative velocity in rocket-body coordinates
v_body = [ ...
    v_air_x*sin(theta_i)+V(i)*cos(theta_i), ...
    v_air_x*cos(theta_i)-V(i)*sin(theta_i),0];

v_body_hat = v_body/max(speed_air,1e-9);

F_aero = -D_drag_i*v_body_hat;
M_pitch_aero = 0;
M_active_fin = 0;

if speed_air > 1 && V(i) > 0

    % Fin aerodynamic-center positions relative to CoM
    r_fin = [ ...
        repmat(x_CoM(i)-x_fin_CP_i,n_fins,1), ...
        fin_radial_cp*cos(fin_phi), ...
        fin_radial_cp*sin(fin_phi)];

    % Individual local airflow, including pitch rotation
    w_body = [0 0 omega(i)];

    v_fin = repmat(v_body,n_fins,1) ...
        +cross(repmat(w_body,n_fins,1),r_fin,2);

    U_fin = vecnorm(v_fin,2,2);
    v_fin_hat = v_fin./max(U_fin,1e-9);
    q_fin = 0.5*rho(i)*U_fin.^2;
    M_fin = U_fin/a_i;

    % Individual fin incidence and commanded deflection
    alpha_fin_flow = atan2( ...
        sum(v_fin.*fin_normal,2),v_fin(:,1));

    delta_fin = delta(i)*fin_normal(:,2);

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

    % Nose normal force
    lift_body = [-v_body_hat(2),v_body_hat(1),0];

    F_nose = q_dyn(i)*A_ref*CN_nose_i ...
        *alpha_i*lift_body;

    r_nose = [x_CoM(i)-x_CP_nose,0,0];

    % Body crossflow force
    sa = sin(alpha_i);

    F_cross = [0, ...
        q_dyn(i)*Cd_body_cross*A_side_body ...
        *sa*abs(sa),0];

    r_cross = [x_CoM(i)-x_CP_cross,0,0];

    % Total aerodynamic force
    F_aero = -D_nonfin*v_body_hat ...
        +F_nose+F_cross+sum(F_fin,1);

    % Total aerodynamic pitching moment
    M_fin_total = sum(cross(r_fin,F_fin,2),1);

    M_nose = cross(r_nose,F_nose);
    M_cross = cross(r_cross,F_cross);

    M_pitch_aero = M_fin_total(3) ...
        +M_nose(3)+M_cross(3);

    % Active control-fin moment for structural excitation
    M_control_fin = sum(cross(r_fin,F_control,2),1);
    M_active_fin = M_control_fin(3);
end

% Actual drag opposing the air-relative velocity
drag_force(i) = -dot(F_aero,v_body_hat);

    % Translational dynamics

    ax = (thrust*sin(theta_i) ...
        +F_aero(1)*sin(theta_i) ...
        +F_aero(2)*cos(theta_i))/m_total(i);

    az = (thrust*cos(theta_i) ...
        +F_aero(1)*cos(theta_i) ...
        -F_aero(2)*sin(theta_i))/m_total(i)-g;

    vx(i+1) = vx(i)+ax*dt;
    V(i+1) = V(i)+az*dt;

    x_pos(i+1) = x_pos(i)+vx(i+1)*dt;
    h(i+1) = max(0,h(i)+V(i+1)*dt);

    if h(i+1)==0 && V(i+1)<0
        V(i+1) = 0;
        vx(i+1) = 0;
    end

    % Sensor model and notch filter

    theta_sensor = theta(i)+q_flex(i);
    omega_sensor = omega(i)+qdot_flex(i);

    switch sensor_mode
        case 0  % Ideal reference sensor
            theta_raw = theta(i);
            omega_raw = omega(i);

        case 1  % Flexible sensor without compensation
            theta_raw = theta_sensor;
            omega_raw = omega_sensor;

        case 2  % Flexible sensor with compensation
            theta_raw = theta_sensor-q_flex_est(i);
            omega_raw = omega_sensor-qdot_flex_est(i);
    end

    theta_m = bN(1)*theta_raw ...
        +bN(2)*angle_x(1)+bN(3)*angle_x(2) ...
        -aN(2)*angle_y(1)-aN(3)*angle_y(2);

    omega_m = bN(1)*omega_raw ...
        +bN(2)*rate_x(1)+bN(3)*rate_x(2) ...
        -aN(2)*rate_y(1)-aN(3)*rate_y(2);

    angle_x = [theta_raw angle_x(1)];
    angle_y = [theta_m angle_y(1)];

    rate_x = [omega_raw rate_x(1)];
    rate_y = [omega_m rate_y(1)];

    % True-vertical PD controller

    % Control during powered ascent and coast
    if V(i) > 1 && speed_air > 1
        theta_ref = 0;

        theta_error = atan2( ...
            sin(theta_m-theta_ref), ...
            cos(theta_m-theta_ref));

        cmd = -Kp_fin*theta_error-Kd_fin*omega_m;
    else
        cmd = 0;
    end

    % Fin actuator dynamics

    cmd = max(-delta_max,min(delta_max,cmd));

    delta(i+1) = delta(i) ...
        +(cmd-delta(i))*dt/tau_s;

    % Pitch moments

    torque_fin(i) = M_active_fin;
    torque_total(i) = M_pitch_aero;

    % Rigid-body pitch

    pitch_accel = torque_total(i)/I_total(i);

    omega(i+1) = omega(i)+pitch_accel*dt;
    theta(i+1) = theta(i)+omega(i+1)*dt;

    % True and estimated bending

    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;

    % Estimate fin torque from aerodynamic conditions
    torque_est = q_dyn(i)*S_fin*CL_max ...
        *tanh(CL_alpha_i*delta(i)/CL_max)*l_fin(i);

    torque_est_history(i) = torque_est;

    % Estimated flexible-mode dynamics
    qdd_est = wn_flex_est^2 ...
        *(K_flex_est*torque_est-q_flex_est(i)) ...
        -2*zeta_flex_est*wn_flex_est*qdot_flex_est(i);

    qdot_flex_est(i+1) = qdot_flex_est(i)+qdd_est*dt;
    q_flex_est(i+1) = q_flex_est(i) ...
        +qdot_flex_est(i+1)*dt;

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
q_dyn(end) = 0.5*rho(end)*hypot(vx(end),V(end))^2;

%% Summary diagnostics

burn_time = m_fuel0/mdot;
burn_idx = find(t <= burn_time);
[apogee,idx_apogee] = max(h);
[maxQ,k] = max(q_dyn(burn_idx));
idx_maxQ = burn_idx(k);

air_speed = hypot(vx-wind_history,V);
air_temp = max(216.65,288.15-0.0065*h);
Mach = air_speed./sqrt(1.4*287.05*air_temp);
[M_max,idx_M] = max(Mach(1:idx_apogee));

gamma = atan2(vx-wind_history,V);
alpha_AoA = atan2(sin(theta-gamma),cos(theta-gamma));
valid_AoA = t <= burn_time & V > 0 & air_speed > 20;
gust_idx = t >= gust_start & t <= gust_start+gust_duration;


fprintf('\nFLIGHT AND CONTROL SUMMARY\n');
fprintf('Apogee: %.1f m at %.1f s\n',apogee,t(idx_apogee));
fprintf('Maximum dynamic pressure: %.1f kPa at %.2f s\n',maxQ/1000,t(idx_maxQ));
fprintf('Maximum Mach: %.3f at %.2f s\n',M_max,t(idx_M));
fprintf('Downrange at apogee: %.2f m\n',x_pos(idx_apogee));
fprintf('Maximum powered AoA: %.3f deg\n',max(abs(rad2deg(alpha_AoA(valid_AoA)))));
fprintf('Pitch at burnout: %.3f deg\n',rad2deg(theta(burn_idx(end))));
fprintf('Maximum pitch during gust: %.3f deg\n',max(abs(rad2deg(theta(gust_idx)))));
fprintf('Peak powered fin deflection: %.3f deg\n',max(abs(rad2deg(delta(burn_idx)))));
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
title('Rocket Downrange Position')
xlim([0 40])
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
    t(burn_idx),rad2deg(gamma(burn_idx)), ...
    'LineWidth',1.2)
legend('Rocket pitch','Flight-path angle')
xlabel('Time (s)')
ylabel('Angle from vertical (deg)')
title('Powered Ascent: Attitude vs Trajectory')
grid on

% 7. True rocket pitch
figure
plot(t,rad2deg(theta))
xlabel('Time (s)')
ylabel('Pitch Angle (deg)')
title('Rocket Pitch')
xlim([0 15])
grid on

% 8. Fin deflection
figure
plot(t,rad2deg(delta))
xlabel('Time (s)')
ylabel('Fin Deflection (deg)')
title('Servo Fin Movement')
xlim([0 15])
grid on

% 9. Bending state
figure
plot(t,rad2deg(q_flex))
xlabel('Time (s)')
ylabel('Structural Bending (deg)')
title('Structural Bending')
xlim([0 15])
grid on

% 10. Pitch moment components
figure
plot(t,torque_fin,t,torque_total-torque_fin,t,torque_total)
legend('Active fin','Passive restoring','Total')
xlabel('Time (s)')
ylabel('Pitch Moment (N m)')
title('Aerodynamic Pitch Moments')
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
plot(t,rad2deg(theta),'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Pitch relative to vertical (deg)')
title('Rocket Attitude Under Wind Gust')
xlim([5 15])
grid on

% 21. Structural sensor measurement: raw, not feedback after compensation
figure
plot(t,rad2deg(theta),'LineWidth',1.3)
hold on
plot(t,rad2deg(theta+q_flex),'--','LineWidth',1.3)
legend('True rocket pitch','Raw sensor pitch')
xlabel('Time (s)')
ylabel('Pitch angle (deg)')
title('Structural Bending and Attitude Measurement Error')
xlim([5 15])
grid on

% 22. Torque-estimator verification
figure
plot(t(gust_idx),torque_fin(gust_idx),'LineWidth',1.3)
hold on
plot(t(gust_idx),torque_est_history(gust_idx),'--','LineWidth',1.3)
legend('Actual active fin torque','Estimated torque')
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
title('Reduced-Order Closed-Loop Poles (Not Full-Model Proof)')
grid on

fprintf('\nMACH-DEPENDENT STABILITY (provisional coefficients)\n');
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
