% ROCKET_FLIGHT_STEP5A  Planar rocket flight with component-based axial drag.
% Run rocket_flight_step5a from this folder. SI units, angles in radians.
% Set show_plots = true before running to display the flight and control figure.
% Requires MATLAB and Control System Toolbox for the Tustin notch filter.
% Pitch and flight-path angles are measured from vertical. V is vertical speed.

%% 1. ADJUSTABLE PARAMETERS

% Simulation
if ~exist('show_plots','var')
    show_plots = false;
end
dt = 0.00025;
t = 0:dt:120;
n = length(t);

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
roughness_selected = 5e-6;  % m, equivalent surface roughness

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

% Rail buttons (provisional)
n_rail_buttons = 2;
D_rail_button = 0.008;

%% 2. GEOMETRY AND MASS PROPERTIES

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

CN_alpha = CN_total;

% Aerodynamic moment arms
l_fin = x_fin_CP-x_CoM;

% Static stability
static_margin = (x_CP-x_CoM)/D_body;

% Axial drag geometry: cylinder, conical nose, and both faces of all fins.

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

A_base = A_ref;
A_motor = pi*D_motor^2/4;
A_button = pi*D_rail_button^2/4;

%% 3. INITIALIZE STATES

% Flight
V = zeros(1,n);
h = zeros(1,n);

% Horizontal motion
vx = zeros(1,n);
x_pos = zeros(1,n);

rho = zeros(1,n);
q_dyn = zeros(1,n);
Mach = zeros(1,n);
drag_force = zeros(1,n);

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

for i = 1:n
    % 4.1  Atmosphere and airspeed

    rho(i) = rho0*exp(-h(i)/H);

    speed_air = hypot(vx(i),V(i));

    q_dyn(i) = 0.5*rho(i)*speed_air^2;

    % 4.2  Flight-path angle and angle of attack

    gamma_i = atan2(vx(i),V(i));

    alpha_i = atan2( ...
        sin(theta(i)-gamma_i), ...
        cos(theta(i)-gamma_i));

    % 4.3  Aerodynamic normal forces

    F_passive = 0;
    F_active = 0;
    F_normal_i = 0;

    % Apply simplified normal-force model during ascent
    if speed_air > 1 && V(i) > 0
        F_passive = q_dyn(i)*A_ref ...
            *CN_alpha*alpha_i;

        CL_i = CL_max*tanh( ...
            CL_alpha*delta(i)/CL_max);

        F_active = -q_dyn(i)*S_fin*CL_i;

        F_normal_i = F_passive+F_active;
    end

    % 4.4  Propulsion and drag

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
    Mach(i) = Mach_i;
    if i == n
        break  % Final state has air-data histories but no subsequent update.
    end

    % Preserve existing flight-equation interface
    drag_factor = D_drag_i/max(speed_air,1e-9);

    % 4.5  Aerodynamic force directions

    if speed_air > 1
        nx = V(i)/speed_air;
        nz = -vx(i)/speed_air;
    else
        nx = 0;
        nz = 0;
    end

    % 4.6  Two-dimensional flight dynamics

    ax = (thrust*sin(theta(i)) ...
        -drag_factor*vx(i) ...
        +F_normal_i*nx)/m_total(i);

    az = (thrust*cos(theta(i)) ...
        -drag_factor*V(i) ...
        +F_normal_i*nz)/m_total(i)-g;

    vx(i+1) = vx(i)+ax*dt;
    V(i+1) = V(i)+az*dt;

    x_pos(i+1) = x_pos(i)+vx(i+1)*dt;
    h(i+1) = max(0,h(i)+V(i+1)*dt);

    if h(i+1)==0 && V(i+1)<0
        V(i+1) = 0;
        vx(i+1) = 0;
    end

    % 4.7  Sensor measurements

    theta_raw = theta(i) + q_flex(i);
    omega_raw = omega(i) + qdot_flex(i);

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

    % 4.8  PD controller targeting vertical attitude

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

    % 4.9  Fin deflection limit and servo

    cmd = max(-delta_max,min(delta_max,cmd));

    delta(i+1) = delta(i) ...
        +(cmd-delta(i))*dt/tau_s;

    % 4.10  Aerodynamic pitch moments

    torque_fin(i) = -F_active*l_fin(i);

    torque_passive(i) = ...
        -F_passive*(x_CP-x_CoM(i));

    if speed_air > 1 && V(i) > 0
        arm = x_CP-x_CoM(i);

        Kq_i = q_dyn(i)*A_ref ...
            *CN_alpha*arm^2/speed_air;

        torque_passive(i) = ...
            torque_passive(i)-Kq_i*omega(i);
    end

    torque_total(i) = ...
        torque_fin(i)+torque_passive(i);

    % 4.11  Rigid-body rotation

    pitch_accel = torque_total(i)/I_total(i);

    omega(i+1) = omega(i)+pitch_accel*dt;
    theta(i+1) = theta(i)+omega(i+1)*dt;

    % 4.12  Structural flexibility

    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;
end

%% 5. FLIGHT AND CONTROL RESULTS

[apogee,idx_apogee] = max(h);
burn_time = m_fuel0/mdot;
idx_burn = find(t <= burn_time,1,'last');
% Match the original Max Q window, including the burnout sample.
[maxQ,idx_maxQ] = max(q_dyn(1:idx_burn));
[M_max,idx_M] = max(Mach(1:idx_apogee));
max_bend_deg = rad2deg(max(abs(q_flex)));
Cd_at_maxQ = drag_force(idx_maxQ)/max(maxQ*A_ref,1e-9);

% Retain the original per-fin hinge-load estimate for servo sizing.
gamma = atan2(vx,V);
alpha_AoA = atan2(sin(theta-gamma),cos(theta-gamma));
F_fin_each = q_dyn.*(A_ref*(CN_fins/n_fins).*alpha_AoA ...
    -A_fin*CL_max.*tanh(CL_alpha*delta/CL_max));
torque_hinge = F_fin_each*hinge_arm;
torque_hinge(V <= 0) = 0;
max_hinge = max(abs(torque_hinge));

fprintf('\nFLIGHT RESULTS\nApogee: %.2f m at %.3f s\nPowered-flight Max Q: %.3f kPa at %.3f s\n', ...
    apogee,t(idx_apogee),maxQ/1000,t(idx_maxQ));
fprintf(['\nCONTROL RESULTS\nPitch at burnout: %+.4f deg\nPeak absolute bending: %.4f deg\n' ...
    'Peak hinge torque per fin: %.3f N m (available %.2f N m)\n'], ...
    rad2deg(theta(idx_burn)),max_bend_deg,max_hinge,M_servo_available);
fprintf('\nMACH NUMBER RESULTS\nPeak ascent Mach: %.5f at %.3f s\n',M_max,t(idx_M));
fprintf('\nINTEGRATED DRAG MODEL\nDrag at powered-flight Max Q: %.3f N\nEffective Cd: %.5f\n', ...
    drag_force(idx_maxQ),Cd_at_maxQ);

%% 6. OPTIONAL FLIGHT AND CONTROL PLOTS

if show_plots
    figure('Name','Step 5A flight and control');
    tiledlayout(2,2);
    nexttile
    plot(t,h)
    ylabel('Altitude (m)')
    xlabel('Time (s)')
    grid on

    nexttile
    plot(t,q_dyn/1000)
    ylabel('Dynamic pressure (kPa)')
    xlabel('Time (s)')
    grid on

    ascent = 1:idx_apogee;
    nexttile
    plot(t(ascent),rad2deg(theta(ascent)),t(ascent),rad2deg(delta(ascent)))
    ylabel('Angle (deg)')
    xlabel('Ascent time (s)')
    legend('Pitch','Fin deflection')
    grid on

    nexttile
    plot(t(ascent),rad2deg(q_flex(ascent)))
    ylabel('Structural bending (deg)')
    xlabel('Ascent time (s)')
    grid on
end
