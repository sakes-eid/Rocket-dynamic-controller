
clear; clc; close all

%% 1. ADJUSTABLE PARAMETERS

% Simulation
dt = 0.00025;
t = 0:dt:120;
n = length(t);

% Horizontal wind gust
wind_peak = -15;       % m/s
gust_start = 7;      % s
gust_duration = 4;   % s

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
Cd0 = 0.5;

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
M_servo_available = 1.96;

% Structural flexibility
f_flex = 4;
zeta_flex = 0.05;
K_flex = 0.001;
wn_flex = 2*pi*f_flex;

% Structural notch filter
f_notch = f_flex;
zeta_notch_zero = 0.01;
zeta_notch_pole = 0.15;
wn_notch = 2*pi*f_notch;

% Mach-dependent drag coefficient
M_table = [0 0.5 0.7 0.8 0.9 1.0 ...
    1.1 1.2 1.4 1.6 2.0];

Cd_ratio = [1 1 1.03 1.10 1.25 1.55 ...
    1.70 1.60 1.35 1.25 1.20];

Cd_curve = pchip(M_table,Cd_ratio);

Cd_body_cross = 1.1;

e_fin = 0.8;  % Provisional fin efficiency

% External protrusions (provisional)
n_rail_buttons = 2;
D_rail_button = 0.008;

roughness_selected = 5e-6;  % m

% 0 = ideal sensor
% 1 = flexible sensor
% 2 = flexible sensor with bending compensation
sensor_mode = 2;



% Working bending estimator with model uncertainty
K_flex_est = 0.8*K_flex;
f_flex_est = 1.05*f_flex;
zeta_flex_est = 1.2*zeta_flex;
wn_flex_est = 2*pi*f_flex_est;


%% 2. GEOMETRY AND MASS PROPERTIES

% Body geometry
A_ref = pi*D_body^2/4;
R_body = D_body/2;

% Fin geometry
A_fin = 0.5*(fin_root+fin_tip)*fin_span;
A_fins_total = n_fins*A_fin;
S_fin = n_pitch_fins*A_fin;

c_mean = (2/3)*(fin_root^2+fin_root*fin_tip ...
    +fin_tip^2)/(fin_root+fin_tip);

x_pivot_chord = fin_pivot_fraction*c_mean;
hinge_arm = c_mean*(fin_load_fraction-fin_pivot_fraction);

A_fin_frontal = fin_span*fin_thickness;
A_front_total = A_ref+n_fins*A_fin_frontal;

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

CN_alpha = CN_total;

% Aerodynamic moment arms
l_fin = x_fin_CP-x_CoM;

% Static stability
static_margin = (x_CP-x_CoM)/D_body;


%% DRAG MODEL PRECOMPUTATION

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

%% 3. INITIALIZE STATES

torque_est_history = zeros(1,n);

q_flex_est = zeros(1,n);
qdot_flex_est = zeros(1,n);

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
torque_passive = zeros(1,n);
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

%% 4. COUPLED FLIGHT AND CONTROL SIMULATION

for i = 1:n-1
    %% 4.1 ATMOSPHERE AND AIRSPEED

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

    %% 4.2 FLIGHT-PATH ANGLE AND ANGLE OF ATTACK

    gamma_i = atan2(v_air_x,V(i));

    alpha_i = atan2( ...
        sin(theta(i)-gamma_i), ...
        cos(theta(i)-gamma_i));


    %% 4.4 PROPULSION AND DRAG

    thrust = Isp*g*mdot*(m_fuel(i)>0);

    T_i = max(216.65,288.15-0.0065*h(i));

    a_i = sqrt(1.4*287.05*T_i);

    Mach_i = speed_air/a_i;

   
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
    
    drag_force(i) = D_drag_i;
    
    % Preserve existing flight-equation interface
    drag_factor = D_drag_i/max(speed_air,1e-9);



%% 4.5 INDIVIDUAL AERODYNAMIC FORCES AND MOMENTS

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
        repmat(x_CoM(i)-x_fin_CP,n_fins,1), ...
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

    CL_passive = (A_ref*CN_fins/S_fin)*alpha_fin_flow;
    CL_control = CL_max*tanh( ...
        CL_alpha*delta_fin/CL_max);

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

    F_nose = q_dyn(i)*A_ref*CN_nose ...
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

    %% 4.6 TWO-DIMENSIONAL FLIGHT DYNAMICS

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

    %% 4.7 SENSOR MEASUREMENTS
    
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

    %% 4.8 PD CONTROLLER - TRUE VERTICAL ATTITUDE

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

    %% 4.9 FIN DEFLECTION LIMIT AND SERVO

    cmd = max(-delta_max,min(delta_max,cmd));

    delta(i+1) = delta(i) ...
        +(cmd-delta(i))*dt/tau_s;

    %% 4.10 AERODYNAMIC PITCH MOMENTS
    
    torque_fin(i) = M_active_fin;
    torque_total(i) = M_pitch_aero;
    torque_passive(i) = torque_total(i)-torque_fin(i);

    %% 4.11 RIGID-BODY ROTATION

    pitch_accel = torque_total(i)/I_total(i);

    omega(i+1) = omega(i)+pitch_accel*dt;
    theta(i+1) = theta(i)+omega(i+1)*dt;

    %% 4.12 STRUCTURAL FLEXIBILITY

    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;


    % Estimate fin torque from aerodynamic conditions
    torque_est = q_dyn(i)*S_fin*CL_max ...
        *tanh(CL_alpha*delta(i)/CL_max)*l_fin(i);
    
    torque_est_history(i) = torque_est;

    % Estimated flexible-mode dynamics
    qdd_est = wn_flex_est^2 ...
        *(K_flex_est*torque_est-q_flex_est(i)) ...
        -2*zeta_flex_est*wn_flex_est*qdot_flex_est(i);

    qdot_flex_est(i+1) = qdot_flex_est(i)+qdd_est*dt;
    q_flex_est(i+1) = q_flex_est(i) ...
        +qdot_flex_est(i+1)*dt;

end

%% 4.13 FINAL ATMOSPHERIC CONDITIONS

rho(end) = rho0*exp(-h(end)/H);

q_dyn(end) = 0.5*rho(end) ...
    *hypot(vx(end),V(end))^2;

%% 5. FLIGHT RESULTS

[apogee,idx_apogee] = max(h);

burn_time = m_fuel0/mdot;
burn_idx = find(t<=burn_time);

[maxQ,k] = max(q_dyn(burn_idx));
idx_maxQ = burn_idx(k);

fprintf('\nFLIGHT RESULTS\n')
fprintf('Apogee: %.1f m at %.1f s\n', ...
    apogee,t(idx_apogee))
fprintf('Max Q: %.1f kPa at %.1f s\n', ...
    maxQ/1000,t(idx_maxQ))

figure
plot(t,h)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('Rocket Altitude')
grid on

figure
plot(t,V)
xlabel('Time (s)')
ylabel('Vertical Velocity (m/s)')
title('Rocket Velocity')
grid on

figure
plot(t,q_dyn/1000)
xlabel('Time (s)')
ylabel('Dynamic Pressure (kPa)')
title('Dynamic Pressure')
grid on

figure
plot(t,x_pos)
xlabel('Time (s)')
ylabel('Horizontal Position (m)')
title('Rocket Downrange Position')
xlim([0 40])
grid on

fprintf('Downrange at apogee: %.2f m\n', ...
    x_pos(idx_apogee))

% Air-relative flight-path angle
gamma = atan2(vx-wind_history,V);
alpha_AoA = atan2(sin(theta-gamma),cos(theta-gamma));

% Powered flight, excluding very low speeds
speed_air = hypot(vx-wind_history,V);
valid = t <= burn_time & V > 0 & speed_air > 20;

figure
plot(t(valid),rad2deg(alpha_AoA(valid)))
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Angle of Attack During Powered Flight')
xlim([0 burn_time])
grid on

fprintf('Maximum powered-flight AoA: %.3f deg\n', ...
    max(abs(rad2deg(alpha_AoA(valid)))));

fprintf('Pitch at burnout: %.3f deg\n',rad2deg(theta(burn_idx(end))))
fprintf('Pitch at apogee: %.3f deg\n',rad2deg(theta(idx_apogee)))

figure
plot(t(burn_idx),rad2deg(theta(burn_idx)), ...
    t(burn_idx),rad2deg(gamma(burn_idx)), ...
    'LineWidth',1.2)

legend('Rocket pitch','Flight-path angle')
xlabel('Time (s)')
ylabel('Angle from vertical (deg)')
title('Powered Ascent: Attitude vs Trajectory')
grid on


% Four-fin body-coordinate geometry
fin_phi = (0:n_fins-1)'*2*pi/n_fins;

fin_span_dir = [zeros(n_fins,1), ...
    cos(fin_phi),sin(fin_phi)];

fin_normal = [zeros(n_fins,1), ...
    -sin(fin_phi),cos(fin_phi)];

fin_radial_cp = R_body + fin_span ...
    *(fin_root+2*fin_tip)/(3*(fin_root+fin_tip));

fin_AR = 2*fin_span^2/A_fin;

% Body crossflow geometry
A_side_nose = 0.5*D_body*L_nose;
A_side_cyl = D_body*(L_body-L_nose);
A_side_body = A_side_nose+A_side_cyl;

x_CP_cross = (A_side_nose*(2*L_nose/3) ...
    +A_side_cyl*(L_nose+L_body)/2)/A_side_body;


%% 6. CONTROL RESULTS

figure
plot(t,rad2deg(theta))
xlabel('Time (s)')
ylabel('Pitch Angle (deg)')
title('Rocket Pitch')
xlim([0 15])
grid on

figure
plot(t,rad2deg(delta))
xlabel('Time (s)')
ylabel('Fin Deflection (deg)')
title('Servo Fin Movement')
xlim([0 15])
grid on

figure
plot(t,rad2deg(q_flex))
xlabel('Time (s)')
ylabel('Structural Bending (deg)')
title('Structural Bending')
xlim([0 15])
grid on

figure
plot(t,torque_fin,t,torque_passive,t,torque_total)
legend('Active fin','Passive restoring','Total')
xlabel('Time (s)')
ylabel('Pitch Moment (N m)')
title('Aerodynamic Pitch Moments')
xlim([0 5])
grid on


fprintf('\nINTEGRATED DRAG MODEL\n')
fprintf('Drag at Max Q: %.2f N\n',drag_force(idx_maxQ))
fprintf('Effective Cd at Max Q: %.4f\n', ...
    drag_force(idx_maxQ)/max(q_dyn(idx_maxQ)*A_ref,1e-9))


%% 7. FROZEN-TIME BODE ANALYSIS AT MAX Q

s = tf('s');

I0 = I_total(idx_maxQ);

Ka = maxQ*S_fin*CL_alpha*l_fin(idx_maxQ);

Kpass = maxQ*A_ref*CN_alpha ...
    *(x_CP-x_CoM(idx_maxQ));

C = Kp_fin+Kd_fin*s;

% Servo
Gs = 1/(tau_s*s+1);

V0 = max(V(idx_maxQ),1);
arm0 = x_CP-x_CoM(idx_maxQ);

Kq0 = maxQ*A_ref*CN_alpha*arm0^2/V0;

Gr = 1/(I0*s^2+Kq0*s+Kpass);

Gf = K_flex*wn_flex^2 ...
    /(s^2+2*zeta_flex*wn_flex*s+wn_flex^2);

L_rigid = C*Gs*Ka*Gr*N_cont;
L_flexible = C*Gs*Ka*(Gr+Gf)*N_cont;

figure
margin(L_flexible)
title('Flexible Rocket - Stability Margins at Max Q')
grid on

fprintf('\nFROZEN-TIME STABILITY\n')
fprintf('Rigid stable: %d\n', ...
    isstable(feedback(L_rigid,1)))
fprintf('Flexible stable: %d\n', ...
    isstable(feedback(L_flexible,1)))

%% 8. FROZEN-TIME STABILITY DURING ASCENT

for tc = [5 6 6.5 7 8 9 10]
    [~,i] = min(abs(t-tc));

    Ka_i = q_dyn(i)*S_fin*CL_alpha*l_fin(i);

    arm_i = x_CP-x_CoM(i);

    Kpass_i = q_dyn(i)*A_ref*CN_alpha*arm_i;

    Vi = max(V(i),1);

    Kq_i = q_dyn(i)*A_ref ...
        *CN_alpha*arm_i^2/Vi;

    Gr_i = 1/(I_total(i)*s^2+Kq_i*s+Kpass_i);

    L_i = C*Gs*Ka_i*(Gr_i+Gf)*N_cont;

    T_i = feedback(L_i,1);

    maxPole = max(real(pole(T_i)));

    fprintf('t=%.1f s | Max pole=%+.4f | Stable=%d\n', ...
        tc,maxPole,isstable(T_i));
end

%% 9. PROJECTED AREA

A_projected = A_ref*abs(cos(alpha_AoA)) ...
    +L_body*D_body*abs(sin(alpha_AoA));

A_fins_projected = ...
    n_fins*A_fin_frontal*abs(cos(alpha_AoA)) ...
    +2*A_fin*abs(sin(alpha_AoA));

A_projected_total = A_projected+A_fins_projected;

figure
plot(t,A_projected,t,A_projected_total)
legend('Body only','Body + four fins')
xlabel('Time (s)')
ylabel('Projected Area (m^2)')
title('Effect of Fins on Projected Area')
xlim([0 5])
grid on

%% 10. STATIC STABILITY

figure
plot(t,static_margin)
yline(0,'--')
xlabel('Time (s)')
ylabel('Static Margin (calibers)')
title('Rocket Static Stability')
xlim([0 15])
grid on

%% 11. PASSIVE AERODYNAMIC MOMENT

figure
plot(t,torque_passive)
xlabel('Time (s)')
ylabel('Pitch Moment (N m)')
title('Passive Aerodynamic Restoring Moment')
xlim([0 5])
grid on

%% 12. GEOMETRY AND CONTROL SUMMARY

fprintf('\nROCKET GEOMETRY\n')
fprintf('Body reference area: %.4f m^2\n',A_ref)
fprintf('One fin area: %.4f m^2\n',A_fin)
fprintf('Four-fin area: %.4f m^2\n',A_fins_total)
fprintf('Pitch-control area: %.4f m^2\n',S_fin)
fprintf('Mean aerodynamic chord: %.4f m\n',c_mean)
fprintf('Pivot chord position: %.4f m\n',x_pivot_chord)
fprintf('Approximate frontal area: %.4f m^2\n',A_front_total)

fprintf('\nGEOMETRY-BASED AERODYNAMICS\n')
fprintf('Nose CP: %.4f m\n',x_CP_nose)
fprintf('Fin aerodynamic center: %.4f m\n',x_fin_CP)
fprintf('Rocket CP: %.4f m\n',x_CP)
fprintf('Total CN_alpha: %.2f per rad\n',CN_total)
fprintf('Initial static margin: %.2f calibers\n', ...
    static_margin(1))
fprintf('Burnout static margin: %.2f calibers\n', ...
    static_margin(burn_idx(end)))

fprintf('\nCONTROL RESULTS\n')
fprintf('Maximum fin deflection: %.4f deg\n', ...
    max(abs(rad2deg(delta(burn_idx)))))
fprintf('Maximum structural bending: %.4f deg\n', ...
    max(abs(rad2deg(q_flex(burn_idx)))))

%% 13. OSCILLATION FREQUENCY DIAGNOSTIC

idx = t>=8 & t<=12;

signal = delta(idx);
signal = signal-mean(signal);

Nfft = length(signal);

Y = abs(fft(signal));
freq = (0:Nfft-1)/(Nfft*dt);

half = 2:floor(Nfft/2);

[~,k] = max(Y(half));
dominant_freq = freq(half(k));

fprintf('\nSTRUCTURAL DIAGNOSTIC\n')
fprintf('Structural frequency: %.2f Hz\n',f_flex)
fprintf('Notch frequency: %.2f Hz\n',f_notch)
fprintf('Dominant fin oscillation: %.2f Hz\n',dominant_freq)

%% 14. ALL-MOVING FIN HINGE TORQUE

% Force on one actuated fin:
% neutral-fin AoA contribution plus control increment

F_fin_each = q_dyn .* ( ...
    A_ref*(CN_fins/n_fins).*alpha_AoA ...
    - A_fin*CL_max.*tanh(CL_alpha*delta/CL_max));

torque_hinge = F_fin_each*hinge_arm;

torque_hinge(V<=0) = 0;

[max_hinge,idx_hinge] = max(abs(torque_hinge));

fprintf('\nALL-MOVING FIN HINGE LOADS\n')
fprintf('Hinge lever arm: %.4f m\n',hinge_arm)
fprintf('Maximum hinge torque per fin: %.3f N m\n',max_hinge)
fprintf('Time of maximum hinge torque: %.2f s\n', ...
    t(idx_hinge))

figure
plot(t,torque_hinge)
xlabel('Time (s)')
ylabel('Hinge Torque per Fin (N m)')
title('Aerodynamic Hinge Torque - All-Moving Fin')
xlim([0 15])
grid on

test_angles = deg2rad([1 3 5]);

% Neutral fin loading from body angle of attack
M_aoa = maxQ*A_ref*(CN_fins/n_fins) ...
    *abs(hinge_arm).*abs(test_angles);

% Additional loading from fin deflection
M_control = maxQ*A_fin*CL_max ...
    .*tanh(CL_alpha*abs(test_angles)/CL_max) ...
    *abs(hinge_arm);

% Conservative case: contributions act in same direction
M_worst = M_aoa+M_control;

fprintf('\nMAX Q HINGE TORQUE STRESS TEST\n')

for k = 1:length(test_angles)
    fprintf(['Angle: %.0f deg | AoA: %.2f | Control: %.2f' ...
        ' | Combined: %.2f N m | Margin: %+.2f N m\n'], ...
        rad2deg(test_angles(k)),M_aoa(k), ...
        M_control(k),M_worst(k), ...
        M_servo_available-M_worst(k));
end

%% 15. FIN SPAN SENSITIVITY STUDY

test_spans = [0.06 0.08 0.10 0.12 0.15];

fprintf('\nFIN SPAN SENSITIVITY STUDY\n')
fprintf('Span | Launch SM | Burnout SM | Kpass/Ka\n')

for span = test_spans
    % Fin area
    A_test = 0.5*(fin_root+fin_tip)*span;

    % Barrowman fin contribution
    l_test = hypot(span,mid_sweep);
    K_test = 1+R_body/(R_body+span);

    CN_test = K_test ...
        *(4*n_fins*(span/D_body)^2) ...
        /(1+sqrt(1+(2*l_test/(fin_root+fin_tip))^2));

    % Center of pressure
    CN_rocket = CN_nose+CN_test;

    CP_test = (CN_nose*x_CP_nose ...
        +CN_test*x_fin_CP)/CN_rocket;

    % Static margins
    SM_launch = (CP_test-x_CoM(1))/D_body;
    SM_burn = (CP_test-x_CoM(burn_idx(end)))/D_body;

    Kpass_test = maxQ*A_ref*CN_rocket ...
        *(CP_test-x_CoM(idx_maxQ));

    Ka_test = maxQ*n_pitch_fins*A_test ...
        *CL_alpha*(x_fin_CP-x_CoM(idx_maxQ));

    fprintf('%.2f | %.2f | %.2f | %.3f\n', ...
        span,SM_launch,SM_burn,Kpass_test/Ka_test);
end

%% 16. FIN POSITION SENSITIVITY STUDY

shifts = [0 0.05 0.10 0.15];

fprintf('\nFIN POSITION SENSITIVITY STUDY\n')
fprintf('Forward shift | Launch SM | Burnout SM | Kpass/Ka\n')

for shift = shifts
    x_fin_test = x_fin_CP-shift;

    x_CP_test = (CN_nose*x_CP_nose ...
        +CN_fins*x_fin_test)/CN_total;

    % Static margins
    SM_launch = (x_CP_test-x_CoM(1))/D_body;
    SM_burn = (x_CP_test-x_CoM(burn_idx(end)))/D_body;

    % Passive and active pitch stiffness at Max Q
    Kpass_test = maxQ*A_ref*CN_total ...
        *(x_CP_test-x_CoM(idx_maxQ));

    Ka_test = maxQ*S_fin*CL_alpha ...
        *(x_fin_test-x_CoM(idx_maxQ));

    fprintf('%.2f m | %.2f | %.2f | %.3f\n', ...
        shift,SM_launch,SM_burn,Kpass_test/Ka_test);
end

%% 17. MACH NUMBER DIAGNOSTIC

T_air = 288.15 - 0.0065*min(h,11000);

a_sound = sqrt(1.4*287.05*T_air);

speed_series = hypot(vx-wind_history,V);

Mach = speed_series./a_sound;

% Maximum Mach during ascent
[M_max,idx_M] = max(Mach(1:idx_apogee));

fprintf('\nMACH NUMBER RESULTS\n')
fprintf('Maximum Mach: %.3f\n',M_max)
fprintf('Time of maximum Mach: %.2f s\n',t(idx_M))
fprintf('Mach at burnout: %.3f\n',Mach(burn_idx(end)))
fprintf('Mach at apogee: %.3f\n',Mach(idx_apogee))

figure
plot(t,Mach,'LineWidth',1.3)
yline(1,'--','Mach 1')
xlabel('Time (s)')
ylabel('Mach number')
title('Rocket Mach Number')
xlim([0 t(idx_apogee)])
grid on

%% 18. MACH-DEPENDENT DRAG COEFFICIENT

M_plot = linspace(0,2,500);
Cd_plot = Cd0*interp1(M_table,Cd_ratio,M_plot,'pchip');

figure
plot(M_plot,Cd_plot,'LineWidth',1.5)
xline(1,'--','Mach 1')
xlabel('Mach number')
ylabel('Drag coefficient Cd')
title('Mach-Dependent Rocket Drag Coefficient')
grid on

Cd_maxMach = Cd0*interp1( ...
    M_table,Cd_ratio,M_max,'pchip');

fprintf('\nMACH-DEPENDENT DRAG\n')
fprintf('Subsonic baseline Cd: %.3f\n',Cd0)
fprintf('Cd at maximum Mach: %.3f\n',Cd_maxMach)
fprintf('Increase over baseline: %.1f%%\n', ...
    100*(Cd_maxMach/Cd0-1))



%% 19. INDIVIDUAL FIN DRAG DIAGNOSTIC

% Body axes: X forward, Y in pitch plane, Z lateral
% Initial illustrative fin profile drag coefficient
Cd_fin_test = 0.02;

% Spanwise center of fin area
r_fin_cp = R_body + fin_span*(fin_root+2*fin_tip) ...
    /(3*(fin_root+fin_tip));

% Four fins spaced 90 degrees apart
fin_azimuth = (0:n_fins-1)'*2*pi/n_fins;

% Evaluate at maximum dynamic pressure
j = idx_maxQ;

% Fin force application points relative to CoM
r_fin_body = [ ...
    repmat(x_CoM(j)-x_fin_CP,n_fins,1), ...
    r_fin_cp*cos(fin_azimuth), ...
    r_fin_cp*sin(fin_azimuth)];

% Rocket velocity in body coordinates
u_body = vx(j)*sin(theta(j)) + V(j)*cos(theta(j));
v_body = vx(j)*cos(theta(j)) - V(j)*sin(theta(j));

% Drag magnitude per fin
D_fin = q_dyn(j)*A_fin*Cd_fin_test;

% Individual drag vectors, opposite local airflow
F_drag_one = -D_fin*[u_body v_body 0] ...
    /max(hypot(u_body,v_body),1e-6);

F_drag_fins = repmat(F_drag_one,n_fins,1);

% Moments about center of mass
M_drag_fins = cross(r_fin_body,F_drag_fins,2);

fprintf('\nINDIVIDUAL FIN DRAG AT MAX Q\n')
fprintf('Fin | Drag (N) | Mx | My | Mz (N m)\n')

for k = 1:n_fins
    fprintf('%d | %.2f | %+.3f | %+.3f | %+.3f\n', ...
        k,norm(F_drag_fins(k,:)),M_drag_fins(k,:));
end

fprintf('Total fin drag: %.2f N\n', ...
    norm(sum(F_drag_fins,1)))

fprintf('Total drag moment [Mx My Mz]: ')
fprintf('%.3f %.3f %.3f N m\n',sum(M_drag_fins,1))



%% 19. REYNOLDS NUMBER DIAGNOSTIC

j = idx_maxQ;

% Atmospheric temperature at Max Q
T_j = max(216.65,288.15-0.0065*h(j));

% Dynamic viscosity (Sutherland's law)
mu_j = 1.716e-5*(T_j/273.15)^1.5 ...
    *(273.15+110.4)/(T_j+110.4);

% Air-relative speed
U_j = hypot(vx(j),V(j));

% Body Reynolds number
Re_body = rho(j)*U_j*L_body/mu_j;

% Fin Reynolds number (mean aerodynamic chord)
Re_fin = rho(j)*U_j*c_mean/mu_j;

fprintf('\nREYNOLDS NUMBERS AT MAX Q\n')
fprintf('Air temperature: %.2f K\n',T_j)
fprintf('Dynamic viscosity: %.3e Pa s\n',mu_j)
fprintf('Body Reynolds number: %.3e\n',Re_body)
fprintf('Fin Reynolds number: %.3e\n',Re_fin)


%% SKIN FRICTION DRAG AT MAX Q

% Local Mach number
M_j = Mach(idx_maxQ);

% Cylindrical body geometry
L_cyl = L_body-L_nose;
S_wet_cyl = pi*D_body*L_cyl;

% Cylinder Reynolds number
Re_cyl = rho(j)*U_j*L_cyl/mu_j;

% Smooth, fully turbulent skin friction
Cf_cyl = 0.455 / ...
    ((log10(Re_cyl))^2.58*(1+0.144*M_j^2)^0.65);

Cf_fin = 0.455 / ...
    ((log10(Re_fin))^2.58*(1+0.144*M_j^2)^0.65);

% Both surfaces of all four fins
S_wet_fins = 2*n_fins*A_fin;

% Skin friction forces
D_cyl_friction = maxQ*Cf_cyl*S_wet_cyl;
D_fins_friction = maxQ*Cf_fin*S_wet_fins;

fprintf('\nSKIN FRICTION AT MAX Q\n')
fprintf('Cylinder wetted area: %.4f m^2\n',S_wet_cyl)
fprintf('Fin wetted area: %.4f m^2\n',S_wet_fins)
fprintf('Cylinder Cf: %.5f\n',Cf_cyl)
fprintf('Fin Cf: %.5f\n',Cf_fin)
fprintf('Cylinder friction drag: %.2f N\n',D_cyl_friction)
fprintf('Four-fin friction drag: %.2f N\n',D_fins_friction)



%% NOSE CONE SKIN FRICTION AT MAX Q

% Provisional straight conical nose
L_slant = hypot(L_nose,R_body);

% Conical wetted area
S_wet_nose = pi*R_body*L_slant;

% Nose Reynolds number
Re_nose = rho(j)*U_j*L_nose/mu_j;

% Compressible turbulent skin-friction coefficient
Cf_nose = 0.455 / ...
    ((log10(Re_nose))^2.58*(1+0.144*M_j^2)^0.65);

% Nose skin-friction drag
D_nose_friction = maxQ*Cf_nose*S_wet_nose;

% Total estimated skin-friction drag
D_friction_total = D_cyl_friction ...
    +D_fins_friction+D_nose_friction;

fprintf('\nNOSE AND TOTAL SKIN FRICTION\n')
fprintf('Nose wetted area: %.4f m^2\n',S_wet_nose)
fprintf('Nose Cf: %.5f\n',Cf_nose)
fprintf('Nose friction drag: %.2f N\n',D_nose_friction)
fprintf('Total skin-friction drag: %.2f N\n',D_friction_total)



%% SURFACE ROUGHNESS SENSITIVITY

% Representative surface roughness heights (micrometers)
roughness_um = [0.5 5 20 100 200];
roughness = roughness_um*1e-6;

% Critical Reynolds number for roughness effects
Recrit_body = 51*(roughness/L_body).^(-1.039);
Recrit_fin = 51*(roughness/c_mean).^(-1.039);

fprintf('\nSURFACE ROUGHNESS SENSITIVITY\n')
fprintf('Roughness | Body affected | Fins affected\n')

for k = 1:length(roughness)
    fprintf('%6.1f um | %13d | %13d\n', ...
        roughness_um(k), ...
        Re_body > Recrit_body(k), ...
        Re_fin > Recrit_fin(k));
end



%% ROUGHNESS-CORRECTED SKIN FRICTION

% Characteristic lengths, Reynolds numbers, wetted areas
L_parts = [L_cyl L_nose c_mean];
Re_parts = [Re_cyl Re_nose Re_fin];
S_parts = [S_wet_cyl S_wet_nose S_wet_fins];

% Smooth turbulent coefficients
Cf_smooth = 1./(1.5*log(Re_parts)-5.6).^2;

% Compressibility correction
if M_j < 1
    Cf_smooth = Cf_smooth*(1-0.1*M_j^2);
    rough_factor = 1-0.1*M_j^2;
else
    Cf_smooth = Cf_smooth/(1+0.15*M_j^2)^0.58;
    rough_factor = 1/(1+0.18*M_j^2);
end

fprintf('\nROUGHNESS-CORRECTED FRICTION AT MAX Q\n')
fprintf('Roughness | Cylinder | Nose | Fins | Total (N)\n')

for k = 1:length(roughness_um)

    % Roughness-limited skin friction
    Cf_rough = 0.032*(roughness(k)./L_parts).^0.2;
    Cf_rough = Cf_rough*rough_factor;

    % Use the larger smooth or roughness-limited value
    Cf_effective = max(Cf_smooth,Cf_rough);

    
    % Empirical body and fin geometry corrections
    K_body = 1 + D_body/(2*L_body);
    K_fin = 1 + 2*fin_thickness/c_mean;
    
    K_parts = [K_body K_body K_fin];
    
    D_parts = maxQ*Cf_effective.*S_parts.*K_parts;


    fprintf('%6.1f um | %8.2f | %5.2f | %5.2f | %8.2f\n', ...
        roughness_um(k),D_parts(1),D_parts(2), ...
        D_parts(3),sum(D_parts));
end



%% NOSE PRESSURE AND WAVE DRAG

% Conical nose geometry
epsilon = atan(R_body/L_nose);
s_eps = sin(epsilon);

% Low-speed and Mach 1 pressure coefficients
Cd_nose0 = 0.8*s_eps^2;
Cd_nose1 = s_eps;

% Slope of pressure coefficient at Mach 1
dCd1 = 4/(1.4+1)*(1-0.5*Cd_nose1);

% Supersonic boundary at Mach 1.3
M2 = 1.3;
Cd_nose2 = 2.1*s_eps^2 ...
    +0.5*s_eps/sqrt(M2^2-1);

dCd2 = -0.5*s_eps*M2/(M2^2-1)^1.5;

% Evaluate at Max Q and across a Mach range
M_eval = [Mach(idx_maxQ) linspace(0,2,401)];
Cd_nose_pressure = zeros(size(M_eval));

% Subsonic empirical interpolation
sub = M_eval <= 1;
b = dCd1/(Cd_nose1-Cd_nose0);

Cd_nose_pressure(sub) = Cd_nose0 ...
    +(Cd_nose1-Cd_nose0)*M_eval(sub).^b;

% Transonic cubic Hermite interpolation
trans = M_eval > 1 & M_eval < M2;
u = (M_eval(trans)-1)/(M2-1);
dM = M2-1;

Cd_nose_pressure(trans) = ...
    (2*u.^3-3*u.^2+1)*Cd_nose1 ...
    +(u.^3-2*u.^2+u)*dM*dCd1 ...
    +(-2*u.^3+3*u.^2)*Cd_nose2 ...
    +(u.^3-u.^2)*dM*dCd2;

% Supersonic pressure-drag correlation
sup = M_eval >= M2;

Cd_nose_pressure(sup) = 2.1*s_eps^2 ...
    +0.5*s_eps./sqrt(M_eval(sup).^2-1);

% Nose pressure drag at Max Q
Cd_nose_maxQ = Cd_nose_pressure(1);
D_nose_pressure = maxQ*A_ref*Cd_nose_maxQ;

fprintf('\nNOSE PRESSURE DRAG AT MAX Q\n')
fprintf('Mach: %.3f\n',Mach(idx_maxQ))
fprintf('Nose pressure Cd: %.4f\n',Cd_nose_maxQ)
fprintf('Nose pressure drag: %.2f N\n',D_nose_pressure)

figure
plot(M_eval(2:end),Cd_nose_pressure(2:end), ...
    'LineWidth',1.5)
xline(1,'--','Mach 1')
xlabel('Mach number')
ylabel('Nose pressure drag coefficient')
title('Conical Nose Pressure Drag vs Mach')
grid on



%% BASE DRAG DIAGNOSTIC

% Full rocket base area
A_base = pi*R_body^2;

% Empirical base drag coefficient
Cd_base = zeros(size(Mach));

sub = Mach <= 1;
Cd_base(sub) = 0.12 + 0.13*Mach(sub).^2;
Cd_base(~sub) = 0.25./Mach(~sub);

% Full-area base drag (power-off reference)
D_base_full = q_dyn.*Cd_base*A_base;


% Powered and unpowered base drag
A_motor = pi*D_motor^2/4;

D_base_actual = D_base_full;

powered = m_fuel > 0;

D_base_actual(powered) = q_dyn(powered) ...
    .*Cd_base(powered)*max(A_base-A_motor,0);

i_on = find(t < burn_time,1,'last');
i_off = find(t >= burn_time,1,'first');

fprintf('\nBASE DRAG WITH MOTOR GEOMETRY\n')
fprintf('Motor diameter: %.0f mm\n',D_motor*1000)
fprintf('Before burnout: %.2f N\n',D_base_actual(i_on))
fprintf('After burnout: %.2f N\n',D_base_actual(i_off))


% Evaluate at Max Q and during coast
j_coast = find(t >= burn_time+2,1);

fprintf('\nBASE DRAG DIAGNOSTIC\n')
fprintf('Base area: %.5f m^2\n',A_base)
fprintf('Full-area Cd at Max Q: %.4f\n',Cd_base(idx_maxQ))
fprintf('Full-area drag at Max Q: %.2f N\n', ...
    D_base_full(idx_maxQ))
fprintf('Coast Mach at 12 s: %.3f\n',Mach(j_coast))
fprintf('Coast base drag at 12 s: %.2f N\n', ...
    D_base_full(j_coast))



%% POWER-ON BASE DRAG SENSITIVITY

% Samples immediately before and after burnout
i_on = find(t < burn_time,1,'last');
i_off = find(t >= burn_time,1,'first');

% Motor-to-body diameter ratios (sensitivity study)
motor_ratio = [0 0.25 0.50 0.75 1.00];

% Effective exposed base area and drag during thrust
A_base_on = A_base*(1-motor_ratio.^2);
D_base_on = q_dyn(i_on)*Cd_base(i_on)*A_base_on;

fprintf('\nBASE DRAG ACROSS BURNOUT\n')
fprintf('Last powered sample: %.5f s\n',t(i_on))
fprintf('First unpowered sample: %.5f s\n',t(i_off))
fprintf('Full-area drag after burnout: %.2f N\n', ...
    D_base_full(i_off))

fprintf('\nMotor/body diameter | Power-on base drag (N)\n')

for k = 1:length(motor_ratio)
    fprintf('%18.2f | %22.2f\n', ...
        motor_ratio(k),D_base_on(k));
end



%% FIN LEADING-EDGE PRESSURE DRAG

j = idx_maxQ;
M = Mach(j);

% Leading-edge sweep angle
sweep_LE = atan(fin_tip_LE_offset/fin_span);

% Rounded leading-edge pressure coefficient
if M < 0.9
    Cd_LE_perp = (1-M^2)^(-0.417)-1;
elseif M < 1
    Cd_LE_perp = 1-1.785*(M-0.9);
else
    Cd_LE_perp = 1.214-0.502/M^2+0.1095/M^4;
end

% Sweep correction
Cd_LE = Cd_LE_perp*cos(sweep_LE)^2;

% Tapered trailing-edge assumption
Cd_TE = 0;

% Pressure drag coefficient
Cd_fin_pressure = Cd_LE+Cd_TE;

% Frontal area of one fin
A_fin_edge = fin_span*fin_thickness;

% Drag per fin and total drag
D_fin_pressure = q_dyn(j)*Cd_fin_pressure*A_fin_edge;
D_fins_pressure = n_fins*D_fin_pressure;

fprintf('\nFIN PRESSURE DRAG AT MAX Q\n')
fprintf('Mach: %.3f\n',M)
fprintf('Leading-edge sweep: %.2f deg\n',rad2deg(sweep_LE))
fprintf('Fin pressure Cd: %.4f\n',Cd_fin_pressure)
fprintf('Pressure drag per fin: %.2f N\n',D_fin_pressure)
fprintf('Pressure drag, all fins: %.2f N\n',D_fins_pressure)



%% TIME-DEPENDENT DRAG BUILD-UP

% Atmospheric viscosity
T_drag = max(216.65,288.15-0.0065*h);
mu_drag = 1.716e-5*(T_drag/273.15).^1.5 ...
    .*(273.15+110.4)./(T_drag+110.4);

% Component geometry
L_parts = [L_cyl L_nose c_mean];
S_parts = [S_wet_cyl S_wet_nose S_wet_fins];

K_parts = [K_body K_body K_fin];

D_friction_series = zeros(size(t));

for p = 1:3

    Re = max(rho.*speed_series*L_parts(p) ...
        ./mu_drag,1e4);

    Cf_smooth = 1./(1.5*log(Re)-5.6).^2;
    Cf_rough = 0.032*(roughness_selected/L_parts(p))^0.2;

    smooth_factor = zeros(size(t));
    rough_factor = zeros(size(t));

    sub = Mach < 1;

    smooth_factor(sub) = 1-0.1*Mach(sub).^2;
    rough_factor(sub) = smooth_factor(sub);

    smooth_factor(~sub) = ...
        1./(1+0.15*Mach(~sub).^2).^0.58;

    rough_factor(~sub) = ...
        1./(1+0.18*Mach(~sub).^2);

    Cf = max(Cf_smooth.*smooth_factor, ...
        Cf_rough.*rough_factor);

    D_friction_series = D_friction_series ...
        +q_dyn.*Cf*S_parts(p)*K_parts(p);
end

% Nose pressure drag
Cd_nose_series = interp1( ...
    M_eval(2:end),Cd_nose_pressure(2:end), ...
    min(max(Mach,0),2),'pchip');

D_nose_series = q_dyn.*A_ref.*Cd_nose_series;

% Four-fin pressure drag
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

% Power-on and power-off base drag
D_base_series = q_dyn.*Cd_base ...
    .*(A_base-A_motor*(m_fuel>0));

% Combined preliminary drag
D_model = D_friction_series+D_nose_series ...
    +D_fin_series+D_base_series;

Cd_model = D_model./max(q_dyn*A_ref,1e-9);

% Plot powered flight and early coast
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

fprintf('\nPRELIMINARY DRAG MODEL\n')
fprintf('Before burnout: %.2f N\n',D_model(i_on))
fprintf('After burnout: %.2f N\n',D_model(i_off))
fprintf('Powered Cd: %.4f\n',Cd_model(i_on))



%% BODY ANGLE-OF-ATTACK CROSSFLOW

% Body side-view area: conical nose + cylinder
A_side_nose = 0.5*D_body*L_nose;
A_side_cyl = D_body*(L_body-L_nose);
A_side_body = A_side_nose+A_side_cyl;

% Crossflow center of pressure from nose tip
x_side_nose = 2*L_nose/3;
x_side_cyl = (L_nose+L_body)/2;

x_CP_cross = (A_side_nose*x_side_nose ...
    +A_side_cyl*x_side_cyl)/A_side_body;

% Angle-dependent body crossflow
sa = sin(alpha_AoA);

% Signed force perpendicular to body axis
F_body_cross = q_dyn*Cd_body_cross*A_side_body ...
    .*sa.*abs(sa);

% Component opposing flight velocity
D_body_AoA = abs(F_body_cross.*sa);

% Evaluate actual flight
j = idx_maxQ;

idx_power = find(t<=burn_time & V>20);
[D_peak,k] = max(D_body_AoA(idx_power));
i_peak = idx_power(k);

fprintf('\nBODY ANGLE-OF-ATTACK DRAG\n')
fprintf('Body side area: %.4f m^2\n',A_side_body)
fprintf('Crossflow CP: %.4f m\n',x_CP_cross)
fprintf('AoA at Max Q: %.3f deg\n',rad2deg(alpha_AoA(j)))
fprintf('Crossflow normal force at Max Q: %.3f N\n',F_body_cross(j))
fprintf('Additional drag at Max Q: %.4f N\n',D_body_AoA(j))
fprintf('Maximum powered AoA drag: %.4f N\n',D_peak)
fprintf('Time of maximum: %.2f s\n',t(i_peak))

% Sensitivity to larger angles of attack at Max Q
test_alpha = [0 1 3 5 10];

fprintf('\nBODY AoA SENSITIVITY AT MAX Q\n')
fprintf('AoA | Normal force (N) | Added drag (N)\n')

for a = test_alpha
    s = sind(a);
    F_test = maxQ*Cd_body_cross*A_side_body*s^2;
    D_test = F_test*s;

    fprintf('%3.0f | %16.2f | %14.2f\n', ...
        a,F_test,D_test)
end




%% INDIVIDUAL FIN LIFT AND INDUCED DRAG

j = idx_maxQ;

% Fin geometry and positions
phi = (0:n_fins-1)'*2*pi/n_fins;

s_fin = [zeros(n_fins,1),cos(phi),sin(phi)];
n_fin = [zeros(n_fins,1),-sin(phi),cos(phi)];

r_span = R_body + fin_span*(fin_root+2*fin_tip) ...
    /(3*(fin_root+fin_tip));

r_fin = [repmat(x_CoM(j)-x_fin_CP,n_fins,1), ...
    r_span*cos(phi),r_span*sin(phi)];

% Local airflow, including pitch rotation

% Air-relative velocity at Max Q
v_air_j = vx(j)-wind_history(j);

v_body = [ ...
    v_air_j*sin(theta(j))+V(j)*cos(theta(j)), ...
    v_air_j*cos(theta(j))-V(j)*sin(theta(j)), ...
    0];


w_body = [0 0 omega(j)];

v_local = repmat(v_body,n_fins,1) ...
    +cross(repmat(w_body,n_fins,1),r_fin,2);

U_fin = vecnorm(v_local,2,2);
v_hat = v_local./max(U_fin,1e-9);
q_fin = 0.5*rho(j)*U_fin.^2;

% Individual control deflections
delta_fin = delta(j)*n_fin(:,2);

% Rotated fin chord and surface normal
c_def = [cos(delta_fin), ...
    -sin(delta_fin).*n_fin(:,2), ...
    -sin(delta_fin).*n_fin(:,3)];

n_def = [sin(delta_fin), ...
    cos(delta_fin).*n_fin(:,2), ...
    cos(delta_fin).*n_fin(:,3)];

% Effective incidence of each fin
alpha_fin = atan2( ...
    sum(v_local.*n_def,2), ...
    sum(v_local.*c_def,2));


% Passive fin loading from local airflow
alpha_flow = atan2( ...
    sum(v_local.*n_fin,2),v_local(:,1));

CL_passive_alpha = A_ref*CN_fins/S_fin;
CL_passive = CL_passive_alpha*alpha_flow;

% Incremental lift from individual fin deflection
CL_control = CL_max*tanh( ...
    CL_alpha*delta_fin/CL_max);

% Combined fin lift coefficient
CL_fin = CL_passive+CL_control;


% Lift: perpendicular to local flow and fin span
lift_dir = cross(v_hat,s_fin,2);
lift_dir = lift_dir./max(vecnorm(lift_dir,2,2),1e-9);

F_lift = -(q_fin*A_fin.*CL_fin).*lift_dir;

% Approximate induced drag (body-mirrored fin)
AR_fin = 2*fin_span^2/A_fin;
Cd_ind = CL_fin.^2/(pi*e_fin*AR_fin);

D_ind = q_fin*A_fin.*Cd_ind;
F_ind = -D_ind.*v_hat;

% Combined forces and moments
F_fin = F_lift+F_ind;
M_fin = cross(r_fin,F_fin,2);

fprintf('\nINDIVIDUAL FIN FORCES AT MAX Q\n')
fprintf('Fin | Incidence | Lift (N) | Induced drag (N)\n')

for k = 1:n_fins
    fprintf('%3d | %+9.3f | %8.3f | %16.4f\n', ...
        k,rad2deg(alpha_fin(k)), ...
        norm(F_lift(k,:)),D_ind(k));
end

fprintf('\nFOUR-FIN TOTALS\n')
fprintf('Effective aspect ratio: %.3f\n',AR_fin)
fprintf('Total induced drag: %.4f N\n',sum(D_ind))
fprintf('Fx: %+.3f N\n',sum(F_fin(:,1)))
fprintf('Fy: %+.3f N\n',sum(F_fin(:,2)))
fprintf('Fz: %+.3f N\n',sum(F_fin(:,3)))
fprintf('Roll moment:  %+.3f N m\n',sum(M_fin(:,1)))
fprintf('Yaw moment:   %+.3f N m\n',sum(M_fin(:,2)))
fprintf('Pitch moment: %+.3f N m\n',sum(M_fin(:,3)))
fprintf('Lift-flow projection: %.6f N\n', ...
    sum(sum(F_lift.*v_hat,2)))



%% INDIVIDUAL FIN PROFILE DRAG

% Local Mach and Reynolds numbers
M_fin_local = U_fin/a_sound(j);

Re_fin_local = max( ...
    rho(j)*U_fin*c_mean/mu_drag(j),1e4);

% Smooth and rough skin-friction coefficients
Cf_smooth_fin = 1./ ...
    (1.5*log(Re_fin_local)-5.6).^2;

Cf_rough_fin = 0.032* ...
    (roughness_selected/c_mean)^0.2 ...
    *ones(n_fins,1);

sub = M_fin_local < 1;

Cf_smooth_fin(sub) = Cf_smooth_fin(sub) ...
    .*(1-0.1*M_fin_local(sub).^2);

Cf_smooth_fin(~sub) = Cf_smooth_fin(~sub) ...
    ./(1+0.15*M_fin_local(~sub).^2).^0.58;

Cf_rough_fin(sub) = Cf_rough_fin(sub) ...
    .*(1-0.1*M_fin_local(sub).^2);

Cf_rough_fin(~sub) = Cf_rough_fin(~sub) ...
    ./(1+0.18*M_fin_local(~sub).^2);

Cf_fin_local = max(Cf_smooth_fin,Cf_rough_fin);

% Fin skin-friction drag
K_fin = 1+2*fin_thickness/c_mean;

D_skin_fin = q_fin.*(2*A_fin) ...
    .*Cf_fin_local*K_fin;

% Mach-dependent leading-edge pressure drag
Cd_LE_local = zeros(n_fins,1);

a = M_fin_local < 0.9;
b = M_fin_local >= 0.9 & M_fin_local < 1;
c = M_fin_local >= 1;

Cd_LE_local(a) = ...
    (1-M_fin_local(a).^2).^(-0.417)-1;

Cd_LE_local(b) = ...
    1-1.785*(M_fin_local(b)-0.9);

Cd_LE_local(c) = ...
    1.214-0.502./M_fin_local(c).^2 ...
    +0.1095./M_fin_local(c).^4;

Cd_LE_local = Cd_LE_local*cos(sweep_LE)^2;

D_pressure_fin = q_fin*A_fin_edge.*Cd_LE_local;

% Profile drag vectors
F_profile_fin = -(D_skin_fin+D_pressure_fin).*v_hat;

% Total fin forces: lift + induced drag + profile drag
F_fin_total = F_lift+F_ind+F_profile_fin;

% Moments about center of mass
M_fin_total = cross(r_fin,F_fin_total,2);

fprintf('\nINDIVIDUAL FIN DRAG AT MAX Q\n')
fprintf('Fin | Skin | Pressure | Induced | Total (N)\n')

for k = 1:n_fins
    fprintf('%3d | %5.2f | %8.2f | %7.3f | %7.2f\n', ...
        k,D_skin_fin(k),D_pressure_fin(k), ...
        D_ind(k),D_skin_fin(k)+D_pressure_fin(k)+D_ind(k));
end

fprintf('\nCOMBINED FIN AERODYNAMICS\n')
fprintf('Total fin drag: %.3f N\n', ...
    sum(D_skin_fin+D_pressure_fin+D_ind))

fprintf('Total Fx: %+.3f N\n',sum(F_fin_total(:,1)))
fprintf('Total Fy: %+.3f N\n',sum(F_fin_total(:,2)))
fprintf('Total Fz: %+.3f N\n',sum(F_fin_total(:,3)))

fprintf('Roll moment:  %+.3f N m\n',sum(M_fin_total(:,1)))
fprintf('Yaw moment:   %+.3f N m\n',sum(M_fin_total(:,2)))
fprintf('Pitch moment: %+.3f N m\n',sum(M_fin_total(:,3)))



%% FIN-BODY INTERFERENCE AND PROTRUSION DRAG

j = idx_maxQ;

% Fin-body interference (empirical approximation)
A_root_equiv = R_body*fin_root;

D_interference = q_fin.*2.*Cf_fin_local ...
    .*K_fin*A_root_equiv;

F_interference = -D_interference.*v_hat;

M_interference = cross(r_fin,F_interference,2);

% Rail-button frontal area
A_button = pi*D_rail_button^2/4;

% Stagnation pressure ratio
M = Mach(j);

if M < 1
    stag_ratio = 1+M^2/4+M^4/40;
else
    stag_ratio = 1.84-0.76/M^2 ...
        +0.166/M^4+0.035/M^6;
end

% Blunt circular protrusion coefficient
Cd_button = 0.85*stag_ratio;

% Rail-button drag
D_buttons = n_rail_buttons*maxQ ...
    *A_button*Cd_button;

% Combined additional drag
D_additional = sum(D_interference)+D_buttons;

fprintf('\nINTERFERENCE AND PROTRUSION DRAG\n')
fprintf('Fin-body interference: %.3f N\n', ...
    sum(D_interference))
fprintf('Rail-button drag: %.3f N\n',D_buttons)
fprintf('Additional drag: %.3f N\n',D_additional)

fprintf('\nINTERFERENCE MOMENTS\n')
fprintf('Roll:  %+.4f N m\n',sum(M_interference(:,1)))
fprintf('Yaw:   %+.4f N m\n',sum(M_interference(:,2)))
fprintf('Pitch: %+.4f N m\n',sum(M_interference(:,3)))

% Updated preliminary drag budget
D_total_updated = D_model(i_on) ...
    +sum(D_ind)+D_additional+D_body_AoA(i_on);

fprintf('\nUPDATED DRAG BUDGET AT MAX Q\n')
fprintf('Previous powered baseline: %.2f N\n', ...
    D_model(i_on))
fprintf('Induced fin drag: %.3f N\n',sum(D_ind))
fprintf('Interference + protrusions: %.3f N\n', ...
    D_additional)
fprintf('Body AoA drag: %.4f N\n',D_body_AoA(i_on))
fprintf('Updated estimated total: %.2f N\n', ...
    D_total_updated)
fprintf('Updated Cd: %.4f\n', ...
    D_total_updated/(q_dyn(i_on)*A_ref))



%% 4D AERODYNAMIC FORCE AND MOMENT ASSEMBLY

% Evaluate at Max Q (engine off at this exact timestep)
j = idx_maxQ;

% Body-axis velocity direction
v_body_hat = v_body/max(norm(v_body),1e-9);

% Zero-AoA drag excluding fin skin friction
% Fin friction is already in F_fin_total
D_body_axis = D_friction_series(j)-sum(D_skin_fin) ...
    +D_nose_series(j)+D_base_series(j)+D_buttons;

F_body_axis = -D_body_axis*v_body_hat;

% Nose normal force, perpendicular to airflow
lift_dir_body = [-v_body_hat(2),v_body_hat(1),0];

F_nose_normal = q_dyn(j)*A_ref*CN_nose ...
    *alpha_AoA(j)*lift_dir_body;

r_nose = [x_CoM(j)-x_CP_nose,0,0];

M_nose = cross(r_nose,F_nose_normal);

% Additional nonlinear body crossflow force
F_cross = [0,F_body_cross(j),0];

r_cross = [x_CoM(j)-x_CP_cross,0,0];

M_cross = cross(r_cross,F_cross);

% Four independent fin forces and moments
F_fins = sum(F_fin_total,1);
M_fins = sum(M_fin_total,1);

% Fin-body interference
F_int = sum(F_interference,1);
M_int = sum(M_interference,1);

% Complete aerodynamic force and moment
F_aero_total = F_body_axis+F_nose_normal ...
    +F_cross+F_fins+F_int;

M_aero_total = M_nose+M_cross+M_fins+M_int;

% Effective drag opposing rocket velocity
D_aero_vector = -dot(F_aero_total,v_body_hat);

% Independent scalar drag-budget check
D_aero_budget = D_model(j)+sum(D_ind) ...
    +sum(D_interference)+D_buttons+D_body_AoA(j);

% Compare new normal-force model with existing 2D model
F_normal_old = q_dyn(j)*A_ref*CN_alpha ...
    *alpha_AoA(j) ...
    -q_dyn(j)*S_fin*CL_max ...
    *tanh(CL_alpha*delta(j)/CL_max);

F_normal_new = F_nose_normal(2)+sum(F_lift(:,2));

fprintf('\nAERODYNAMIC FORCE ASSEMBLY\n')
fprintf('Fx: %+.3f N\n',F_aero_total(1))
fprintf('Fy: %+.3f N\n',F_aero_total(2))
fprintf('Fz: %+.3f N\n',F_aero_total(3))
fprintf('Vector-derived drag: %.3f N\n',D_aero_vector)
fprintf('Component drag budget: %.3f N\n',D_aero_budget)

fprintf('\nAERODYNAMIC MOMENT ASSEMBLY\n')
fprintf('Roll Mx:  %+.3f N m\n',M_aero_total(1))
fprintf('Yaw My:   %+.3f N m\n',M_aero_total(2))
fprintf('Pitch Mz: %+.3f N m\n',M_aero_total(3))

fprintf('\nOLD VS NEW AERODYNAMIC MODEL\n')
fprintf('Old normal force: %+.3f N\n',F_normal_old)
fprintf('New normal force: %+.3f N\n',F_normal_new)
fprintf('Old pitch moment: %+.3f N m\n',torque_total(j))
fprintf('New pitch moment: %+.3f N m\n',M_aero_total(3))



%% NORMAL FORCE AND PITCH MOMENT CHECK

% Old model: separate nose, passive fins and control
F_old_nose = q_dyn(j)*A_ref*CN_nose*alpha_AoA(j);
F_old_fin_passive = q_dyn(j)*A_ref*CN_fins*alpha_AoA(j);

F_old_fin_control = -q_dyn(j)*S_fin*CL_max ...
    *tanh(CL_alpha*delta(j)/CL_max);

% Old moment decomposition
M_old_nose = F_old_nose*(x_CoM(j)-x_CP_nose);
M_old_fin = (F_old_fin_passive+F_old_fin_control) ...
    *(x_CoM(j)-x_fin_CP);

M_old_damping = torque_total(j)-M_old_nose-M_old_fin;

% New moment decomposition
M_new_fin_lift = sum(cross(r_fin,F_lift,2),1);
M_new_fin_drag = sum(cross( ...
    r_fin,F_ind+F_profile_fin,2),1);

fprintf('\nFIN FORCE CONSISTENCY\n')
fprintf('Old equivalent fin CL slope: %.3f\n', ...
    A_ref*CN_fins/S_fin)
fprintf('New fin CL slope: %.3f\n',CL_alpha)
fprintf('Old passive fin force: %+.3f N\n',F_old_fin_passive)
fprintf('Old control fin force: %+.3f N\n',F_old_fin_control)
fprintf('New combined fin lift: %+.3f N\n', ...
    sum(F_lift(:,2)))

fprintf('\nPITCH MOMENT BREAKDOWN\n')
fprintf('Old nose: %+.3f N m\n',M_old_nose)
fprintf('Old fins: %+.3f N m\n',M_old_fin)
fprintf('Old damping: %+.3f N m\n',M_old_damping)
fprintf('New nose: %+.3f N m\n',M_nose(3))
fprintf('New fin lift: %+.3f N m\n',M_new_fin_lift(3))
fprintf('New fin drag: %+.3f N m\n',M_new_fin_drag(3))
fprintf('New crossflow: %+.3f N m\n',M_cross(3))
fprintf('New interference: %+.3f N m\n',M_int(3))


%% WIND GUST RESPONSE

gust_idx = t >= gust_start & ...
    t <= gust_start+gust_duration;

fprintf('\nWIND GUST RESPONSE\n')
fprintf('Peak gust: %.2f m/s\n',wind_peak)
fprintf('Maximum pitch during gust: %.3f deg\n', ...
    max(abs(rad2deg(theta(gust_idx)))))
fprintf('Maximum fin deflection during gust: %.3f deg\n', ...
    max(abs(rad2deg(delta(gust_idx)))))
fprintf('Maximum bending during gust: %.3f deg\n', ...
    max(abs(rad2deg(q_flex(gust_idx)))))

figure
plot(t,rad2deg(theta),'LineWidth',1.3)
xlabel('Time (s)')
ylabel('Pitch relative to vertical (deg)')
title('Rocket Attitude Under Wind Gust')
xlim([5 15])
grid on



%% STRUCTURAL SENSOR ERROR DIAGNOSTIC

% Raw sensor reading including structural bending
theta_sensor_raw = theta+q_flex;

% Burnout sample
j = burn_idx(end);

% Gust interval
gust_idx = t>=gust_start & ...
    t<=gust_start+gust_duration;

fprintf('\nSTRUCTURAL SENSOR ERROR\n')
fprintf('True pitch at burnout: %.3f deg\n', ...
    rad2deg(theta(j)))
fprintf('Bending at burnout: %.3f deg\n', ...
    rad2deg(q_flex(j)))
fprintf('Raw sensor pitch at burnout: %.3f deg\n', ...
    rad2deg(theta_sensor_raw(j)))
fprintf('Maximum bending during gust: %.3f deg\n', ...
    max(abs(rad2deg(q_flex(gust_idx)))))

figure
plot(t,rad2deg(theta),'LineWidth',1.3)
hold on
plot(t,rad2deg(theta_sensor_raw),'--','LineWidth',1.3)
legend('True rocket pitch','Raw sensor pitch')
xlabel('Time (s)')
ylabel('Pitch angle (deg)')
title('Structural Bending and Attitude Measurement Error')
xlim([5 15])
grid on



%% CONTROL TORQUE AND FLEXIBILITY COMPARISON

idx_gust = t >= gust_start & ...
    t <= gust_start+gust_duration;

fprintf('\nCONTROL AND BENDING DIAGNOSTIC\n')
fprintf('Sensor mode: %d\n',sensor_mode)
fprintf('Peak active fin torque: %.2f N m\n', ...
    max(abs(torque_fin(idx_gust))))
fprintf('RMS active fin torque: %.2f N m\n', ...
    sqrt(mean(torque_fin(idx_gust).^2)))
fprintf('Maximum bending: %.3f deg\n', ...
    max(abs(rad2deg(q_flex(idx_gust)))))
fprintf('Maximum true pitch: %.3f deg\n', ...
    max(abs(rad2deg(theta(idx_gust)))))
fprintf('Pitch at burnout: %.3f deg\n', ...
    rad2deg(theta(burn_idx(end))))



%% BENDING ESTIMATOR VALIDATION

idx_gust = t>=gust_start & ...
    t<=gust_start+gust_duration;

est_error = q_flex(idx_gust)-q_flex_est(idx_gust);

fprintf('\nBENDING ESTIMATOR RESULTS\n')
fprintf('True peak bending: %.3f deg\n', ...
    max(abs(rad2deg(q_flex(idx_gust)))))
fprintf('Estimated peak bending: %.3f deg\n', ...
    max(abs(rad2deg(q_flex_est(idx_gust)))))
fprintf('RMS estimation error: %.3f deg\n', ...
    rad2deg(sqrt(mean(est_error.^2))))
fprintf('Maximum estimation error: %.3f deg\n', ...
    max(abs(rad2deg(est_error))))



%% BENDING COMPENSATION PERFORMANCE

idx = t>=gust_start & t<=gust_start+gust_duration;
j = burn_idx(end);

fprintf('\nBENDING COMPENSATION RESULTS\n')
fprintf('Sensor mode: %d\n',sensor_mode)
fprintf('Maximum pitch during gust: %.3f deg\n', ...
    max(abs(rad2deg(theta(idx)))))
fprintf('Pitch at burnout: %.3f deg\n', ...
    rad2deg(theta(j)))
fprintf('Peak fin torque: %.2f N m\n', ...
    max(abs(torque_fin(idx))))
fprintf('RMS fin torque: %.2f N m\n', ...
    sqrt(mean(torque_fin(idx).^2)))
fprintf('Maximum bending: %.3f deg\n', ...
    max(abs(rad2deg(q_flex(idx)))))
fprintf('Maximum fin deflection: %.3f deg\n', ...
    max(abs(rad2deg(delta(idx)))))
fprintf('RMS bending estimator error: %.3f deg\n', ...
    rad2deg(sqrt(mean( ...
    (q_flex(idx)-q_flex_est(idx)).^2))))

%% FIN TORQUE ESTIMATOR ACCURACY

idx = t>=gust_start & t<=gust_start+gust_duration;

actual = torque_fin(idx);
estimated = torque_est_history(idx);

fprintf('\nFIN TORQUE ESTIMATOR ACCURACY\n')
fprintf('Actual RMS torque: %.2f N m\n', ...
    sqrt(mean(actual.^2)))
fprintf('Estimated RMS torque: %.2f N m\n', ...
    sqrt(mean(estimated.^2)))
fprintf('Torque RMS error: %.2f N m\n', ...
    sqrt(mean((actual-estimated).^2)))
fprintf('Torque relative RMS error: %.1f%%\n', ...
    100*sqrt(mean((actual-estimated).^2)) ...
    /max(sqrt(mean(actual.^2)),1e-9))

figure
plot(t(idx),actual,'LineWidth',1.3)
hold on
plot(t(idx),estimated,'--','LineWidth',1.3)
legend('Actual active fin torque','Estimated torque')
xlabel('Time (s)')
ylabel('Torque (N m)')
title('Fin Torque Estimation During Gust')
grid on
