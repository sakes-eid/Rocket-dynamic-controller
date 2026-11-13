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

%% Frozen-trajectory numerical stability snapshots
% Capture the actual controller/filter states, not inferred initial conditions.
lin_times = unique([2:0.5:15 9.99 10.01]);
lin_indices = 1+round(lin_times/dt);
lin_n = numel(lin_times);
lin_z = zeros(16,lin_n);
lin_znext = zeros(16,lin_n);
lin_snap = repmat(struct('t',0,'V',0,'h',0,'mfuel',0,...
    'mass',0,'com',0,'inertia',0),1,lin_n);
lin_slot = 1;

% Event flags. Descent begins at the first downward velocity after apogee.
descent_active = false;
apogee_lock_idx = NaN;
impact_idx = NaN;
impact_time = NaN;

%% Coupled flight and control simulation

for i = 1:n-1
    % Lock the two actuated fins at neutral once the ascent reaches apogee.
    if ~descent_active && t(i)>m_fuel0/mdot && ...
            V(i)<=0 && V(i-1)>0 && h(i)>0
        descent_active = true;
        apogee_lock_idx = i;
        delta(i) = 0;
    end

    % Frozen-time Jacobian uses the pre-update flight and filter states.
    lin_capture = lin_slot <= lin_n && i == lin_indices(lin_slot);
    if lin_capture
        lin_z(:,lin_slot) = [vx(i);theta(i);omega(i);delta(i);...
            q_flex(i);qdot_flex(i);q_flex_est(i);qdot_flex_est(i);...
            angle_x(:);angle_y(:);rate_x(:);rate_y(:)];
        lin_snap(lin_slot) = struct('t',t(i),'V',V(i),'h',h(i),...
            'mfuel',m_fuel(i),'mass',m_total(i),...
            'com',x_CoM(i),'inertia',I_total(i));
    end
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

elseif descent_active && speed_air > 1
    % Post-apogee: aerodynamic loading with unpowered, NEUTRAL static fins.
    % Unlike the ascent model, the axial air velocity can be negative
    % (tail-first motion). Avoid the invalid ~180-degree linear nose AoA.
    % The finite-flow normal load opposes the velocity through the fin plane.
    w_body = [0 0 omega(i)];
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

    % Broadside body crossflow acts against local normal airflow.
    v_side = v_body_hat(2);
    F_cross = [0, -q_dyn(i)*Cd_body_cross*A_side_body ...
        *v_side*abs(v_side), 0];
    r_cross = [x_CoM(i)-x_CP_cross,0,0];

    F_aero = -D_axial*v_body_hat+F_cross+sum(F_fin,1);
    M_fin_total = sum(cross(r_fin,F_fin,2),1);
    M_cross = cross(r_cross,F_cross);
    M_pitch_aero = M_fin_total(3)+M_cross(3);
    M_active_fin = 0;
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
    h_trial = h(i)+V(i+1)*dt;
    h(i+1) = max(0,h_trial);

    if descent_active && h_trial <= 0 && h(i)>0
        impact_idx = i+1;
        % Interpolate the crossing within this timestep. Preserve the
        % impact velocity; do not artificially set it to zero.
        impact_fraction = min(max(h(i)/(h(i)-h_trial),0),1);
        impact_time = t(i)+impact_fraction*dt;
        x_pos(i+1) = x_pos(i)+impact_fraction ...
            *(x_pos(i+1)-x_pos(i));
        vx(i+1) = vx(i)+impact_fraction*(vx(i+1)-vx(i));
        V(i+1) = V(i)+impact_fraction*(V(i+1)-V(i));
    elseif h_trial <= 0 && V(i+1)<0
        % Retain v7's launch-pad floor behavior before first ascent.
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

    if descent_active
        % Immediate mechanical lock at zero, with no active fin command.
        delta(i+1) = 0;
    else
        delta(i+1) = delta(i) ...
            +(cmd-delta(i))*dt/tau_s;
    end

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

    if lin_capture
        lin_znext(:,lin_slot) = [vx(i+1);theta(i+1);omega(i+1);...
            delta(i+1);q_flex(i+1);qdot_flex(i+1);...
            q_flex_est(i+1);qdot_flex_est(i+1);...
            angle_x(:);angle_y(:);rate_x(:);rate_y(:)];
        lin_slot = lin_slot+1;
    end
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
    theta=theta(1:n); omega=omega(1:n); delta=delta(1:n);
    q_flex=q_flex(1:n); qdot_flex=qdot_flex(1:n);
    q_flex_est=q_flex_est(1:n); qdot_flex_est=qdot_flex_est(1:n);
    rho=rho(1:n); q_dyn=q_dyn(1:n);
    drag_force=drag_force(1:n); wind_history=wind_history(1:n);
    torque_fin=torque_fin(1:n); torque_total=torque_total(1:n);
    torque_est_history=torque_est_history(1:n);
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
fprintf('\nBALLISTIC DESCENT (FINS LOCKED NEUTRAL, NO PARACHUTE)\n');
if isfinite(apogee_lock_idx)
    fprintf('Fin lock / apogee time: %.3f s\n',t(apogee_lock_idx));
end
if isfinite(impact_idx)
    fprintf('Touchdown time: %.3f s\n',impact_time);
    fprintf('Time from apogee to touchdown: %.2f s\n', ...
        impact_time-t(idx_apogee));
    fprintf('Touchdown horizontal position: %.2f m\n',x_pos(end));
    fprintf('Touchdown vertical velocity: %.2f m/s\n',V(end));
    fprintf('Touchdown horizontal velocity: %.2f m/s\n',vx(end));
    fprintf('Touchdown speed: %.2f m/s\n',hypot(vx(end),V(end)));
    fprintf('Peak descent dynamic pressure: %.1f kPa\n', ...
        max(q_dyn(idx_apogee:end))/1000);
else
    fprintf('No touchdown within simulated %.1f s; extend t if needed.\n',t(end));
end
if isfinite(apogee_lock_idx)
    fprintf('Maximum fin deflection after lock: %.6f deg\n', ...
        max(abs(rad2deg(delta(apogee_lock_idx:end)))));
end
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

%% Numerical frozen-trajectory Jacobian (actual discrete control implementation)
% Linearization includes vx, theta, omega, fin servo, true+estimated flex,
% and both complete IIR notch-filter memories (16 states).
% Altitude/vertical speed/mass and gust profile are frozen at each snapshot.
% This is a local perturbation-growth test, NOT a proof for the changing
% nonlinear trajectory or for unmodeled aerodynamic/structural uncertainty.
lin_P = struct();
lin_P.dt = dt;
lin_P.wind_peak = wind_peak;
lin_P.gust_start = gust_start;
lin_P.gust_duration = gust_duration;
lin_P.rho0 = rho0;
lin_P.H = H;
lin_P.Isp = Isp;
lin_P.g = g;
lin_P.mdot = mdot;
lin_P.drag_L = drag_L;
lin_P.drag_S = drag_S;
lin_P.drag_K = drag_K;
lin_P.Cf_rough0 = Cf_rough0;
lin_P.Cd_n0 = Cd_n0;
lin_P.Cd_n1 = Cd_n1;
lin_P.dCd_n1 = dCd_n1;
lin_P.Cd_n2 = Cd_n2;
lin_P.dCd_n2 = dCd_n2;
lin_P.nose_exponent = nose_exponent;
lin_P.M_n2 = M_n2;
lin_P.s_eps = s_eps;
lin_P.n_fins = n_fins;
lin_P.A_fin = A_fin;
lin_P.A_ref = A_ref;
lin_P.A_fin_edge = A_fin_edge;
lin_P.sweep_LE = sweep_LE;
lin_P.A_base = A_base;
lin_P.A_motor = A_motor;
lin_P.n_rail_buttons = n_rail_buttons;
lin_P.A_button = A_button;
lin_P.K_fin = K_fin;
lin_P.R_body = R_body;
lin_P.fin_root = fin_root;
lin_P.fin_radial_cp = fin_radial_cp;
lin_P.fin_phi = fin_phi;
lin_P.fin_span_dir = fin_span_dir;
lin_P.fin_normal = fin_normal;
lin_P.fin_AR = fin_AR;
lin_P.e_fin = e_fin;
lin_P.c_mean = c_mean;
lin_P.D_body = D_body;
lin_P.x_MAC_LE = x_MAC_LE;
lin_P.x_CP_nose = x_CP_nose;
lin_P.x_CP_cross = x_CP_cross;
lin_P.Cd_body_cross = Cd_body_cross;
lin_P.A_side_body = A_side_body;
lin_P.CN_nose = CN_nose;
lin_P.CN_fins = CN_fins;
lin_P.CL_alpha = CL_alpha;
lin_P.CL_max = CL_max;
lin_P.S_fin = S_fin;
lin_P.Ma_nodes = Ma_nodes;
lin_P.pp_nose = pp_nose;
lin_P.pp_fin = pp_fin;
lin_P.pp_ctrl = pp_ctrl;
lin_P.pp_ac = pp_ac;
lin_P.tau_s = tau_s;
lin_P.Kp_fin = Kp_fin;
lin_P.Kd_fin = Kd_fin;
lin_P.delta_max = delta_max;
lin_P.sensor_mode = sensor_mode;
lin_P.wn_flex = wn_flex;
lin_P.K_flex = K_flex;
lin_P.zeta_flex = zeta_flex;
lin_P.wn_flex_est = wn_flex_est;
lin_P.K_flex_est = K_flex_est;
lin_P.zeta_flex_est = zeta_flex_est;
lin_P.bN = bN;
lin_P.aN = aN;

lin_growth = NaN(1,lin_n);
lin_match_error = NaN(1,lin_n);
lin_spectral_radius = NaN(1,lin_n);
lin_eps_sensitivity = NaN(1,lin_n);
lin_eigs = cell(1,lin_n);
lin_steps = [0.01 1e-5 1e-4 1e-5 1e-5 1e-4 1e-5 1e-4 ...
             1e-5 1e-5 1e-5 1e-5 1e-4 1e-4 1e-4 1e-4]';

for k = 1:lin_n
    z0 = lin_z(:,k);
    f0 = rocketFrozenStep(z0,lin_snap(k),lin_P);
    lin_match_error(k) = max(abs(f0-lin_znext(:,k)));
    if lin_match_error(k)>1e-7
        warning('Snapshot t=%.3f s: copied step does not match flight loop (%.3g).', ...
            lin_times(k),lin_match_error(k));
        continue
    end

    A_discrete = zeros(16);
    A_half = zeros(16);
    for c = 1:16
        e = zeros(16,1);
        e(c) = lin_steps(c);
        plus_step = rocketFrozenStep(z0+e,lin_snap(k),lin_P);
        minus_step = rocketFrozenStep(z0-e,lin_snap(k),lin_P);
        A_discrete(:,c) = (plus_step-minus_step)/(2*lin_steps(c));
        e(c) = lin_steps(c)/2;
        plus_half = rocketFrozenStep(z0+e,lin_snap(k),lin_P);
        minus_half = rocketFrozenStep(z0-e,lin_snap(k),lin_P);
        A_half(:,c) = (plus_half-minus_half)/lin_steps(c);
    end
    lin_eps_sensitivity(k) = abs(log(max(abs(eig(A_discrete)))) ...
        -log(max(abs(eig(A_half)))))/dt;
    poles_d = eig(A_discrete);
    lin_eigs{k} = poles_d;
    lin_spectral_radius(k) = max(abs(poles_d));
    lin_growth(k) = log(lin_spectral_radius(k))/dt;
end

figure
plot(lin_times,lin_growth,'o-','LineWidth',1.2)
hold on
plot(sample_times,max_real_pole,'--','LineWidth',1.1)
yline(0,':','Neutral growth')
legend('16-state numerical Jacobian','Earlier reduced-order model',...
    'Location','best')
xlabel('Flight time (s)')
ylabel('Largest local growth exponent (1/s)')
title('Frozen-Flight Closed-Loop Stability: Numerical vs Reduced-Order')
grid on

fprintf('\nNUMERICAL FROZEN-FLIGHT STABILITY (16-state)\n');
fprintf('V, h, mass and external wind frozen; longitudinal/lateral motion beyond vx excluded.\n');
fprintf('Time (s) | Mach | Growth (1/s) | Step mismatch | Eps-sensitivity\n');
for k = 1:lin_n
    if ismember(round(lin_times(k)*100), ...
            round([2 3 4 5 7 8 9 9.99 10.01 11 12 15]*100))
        j = lin_indices(k);
        fprintf('%8.2f | %5.3f | %+13.4f | %11.2g | %11.2g\n', ...
            lin_times(k),Mach(j),lin_growth(k),...
            lin_match_error(k),lin_eps_sensitivity(k));
    end
end
fprintf('Local positive growth does not establish full-trajectory instability.\n');
fprintf('The Jacobian is valid only where the controls and aerodynamic models are differentiable.\n');




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
yline(0,':')
legend('Vertical velocity','Horizontal velocity','Zero','Location','best')
xlabel('Time (s)')
ylabel('Ground-relative velocity (m/s)')
title('Ballistic Descent Velocities')
grid on

% Passive dynamics, locked fins throughout descent.
figure
subplot(2,1,1)
plot(t(idx_apogee:end),rad2deg(theta(idx_apogee:end)))
ylabel('Body pitch (deg)')
title('Post-Apogee Passive Attitude and Fin Lock')
grid on
subplot(2,1,2)
plot(t(idx_apogee:end),rad2deg(delta(idx_apogee:end)))
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
plot(t(idx_apogee:end),torque_total(idx_apogee:end),'LineWidth',1.3)
yline(0,':')
xlabel('Time (s)')
ylabel('Passive pitch moment (N m)')
title('Post-Apogee Static-Fin and Body Pitching Moment')
grid on

function z1 = rocketFrozenStep(z0,snap,P)
% Exact one-step replay of the active flight-loop equations with frozen
% altitude, vertical speed, mass properties and wind-time input.
% z = [vx theta omega delta flex flexrate flexest flexestrate ...
%      angle_x(2) angle_y(2) rate_x(2) rate_y(2)]'.
    dt = P.dt;
    wind_peak = P.wind_peak;
    gust_start = P.gust_start;
    gust_duration = P.gust_duration;
    rho0 = P.rho0;
    H = P.H;
    Isp = P.Isp;
    g = P.g;
    mdot = P.mdot;
    drag_L = P.drag_L;
    drag_S = P.drag_S;
    drag_K = P.drag_K;
    Cf_rough0 = P.Cf_rough0;
    Cd_n0 = P.Cd_n0;
    Cd_n1 = P.Cd_n1;
    dCd_n1 = P.dCd_n1;
    Cd_n2 = P.Cd_n2;
    dCd_n2 = P.dCd_n2;
    nose_exponent = P.nose_exponent;
    M_n2 = P.M_n2;
    s_eps = P.s_eps;
    n_fins = P.n_fins;
    A_fin = P.A_fin;
    A_ref = P.A_ref;
    A_fin_edge = P.A_fin_edge;
    sweep_LE = P.sweep_LE;
    A_base = P.A_base;
    A_motor = P.A_motor;
    n_rail_buttons = P.n_rail_buttons;
    A_button = P.A_button;
    K_fin = P.K_fin;
    R_body = P.R_body;
    fin_root = P.fin_root;
    fin_radial_cp = P.fin_radial_cp;
    fin_phi = P.fin_phi;
    fin_span_dir = P.fin_span_dir;
    fin_normal = P.fin_normal;
    fin_AR = P.fin_AR;
    e_fin = P.e_fin;
    c_mean = P.c_mean;
    D_body = P.D_body;
    x_MAC_LE = P.x_MAC_LE;
    x_CP_nose = P.x_CP_nose;
    x_CP_cross = P.x_CP_cross;
    Cd_body_cross = P.Cd_body_cross;
    A_side_body = P.A_side_body;
    CN_nose = P.CN_nose;
    CN_fins = P.CN_fins;
    CL_alpha = P.CL_alpha;
    CL_max = P.CL_max;
    S_fin = P.S_fin;
    Ma_nodes = P.Ma_nodes;
    pp_nose = P.pp_nose;
    pp_fin = P.pp_fin;
    pp_ctrl = P.pp_ctrl;
    pp_ac = P.pp_ac;
    tau_s = P.tau_s;
    Kp_fin = P.Kp_fin;
    Kd_fin = P.Kd_fin;
    delta_max = P.delta_max;
    sensor_mode = P.sensor_mode;
    wn_flex = P.wn_flex;
    K_flex = P.K_flex;
    zeta_flex = P.zeta_flex;
    wn_flex_est = P.wn_flex_est;
    K_flex_est = P.K_flex_est;
    zeta_flex_est = P.zeta_flex_est;
    bN = P.bN;
    aN = P.aN;
    i = 1;
    t = [snap.t snap.t+dt];
    h = [snap.h 0];
    V = [snap.V 0];
    vx = [z0(1) 0];
    theta = [z0(2) 0];
    omega = [z0(3) 0];
    delta = [z0(4) 0];
    q_flex = [z0(5) 0];
    qdot_flex = [z0(6) 0];
    q_flex_est = [z0(7) 0];
    qdot_flex_est = [z0(8) 0];
    angle_x = z0(9:10).';
    angle_y = z0(11:12).';
    rate_x = z0(13:14).';
    rate_y = z0(15:16).';
    m_fuel = [snap.mfuel snap.mfuel];
    m_total = [snap.mass snap.mass];
    x_CoM = [snap.com snap.com];
    I_total = [snap.inertia snap.inertia];
    q_dyn = zeros(1,2);
    rho = zeros(1,2);
    x_pos = zeros(1,2);
    wind_history = zeros(1,2);
    CN_nose_hist = zeros(1,2);
    CN_fins_hist = zeros(1,2);
    CL_alpha_hist = zeros(1,2);
    x_fin_CP_hist = zeros(1,2);
    x_CP_hist = zeros(1,2);
    l_fin = zeros(1,2);
    static_margin = zeros(1,2);
    drag_force = zeros(1,2);
    torque_fin = zeros(1,2);
    torque_total = zeros(1,2);
    torque_est_history = zeros(1,2);

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


    z1 = [vx(i+1);theta(i+1);omega(i+1);delta(i+1);...
        q_flex(i+1);qdot_flex(i+1);q_flex_est(i+1);...
        qdot_flex_est(i+1);angle_x(:);angle_y(:);rate_x(:);rate_y(:)];
end
