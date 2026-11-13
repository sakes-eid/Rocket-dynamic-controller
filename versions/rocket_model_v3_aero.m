
clear; clc; close all

%% 1. ADJUSTABLE PARAMETERS

% Simulation
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

% Four all-moving fins
n_fins = 4;
n_pitch_fins = 2;

fin_span = 0.15;
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


%% 2. GEOMETRY AND MASS PROPERTIES

% Body geometry
A_ref = pi*D_body^2/4;
R_body = D_body/2;

% Fin geometry
A_fin = 0.5*(fin_root+fin_tip)*fin_span;
A_fins_total = n_fins*A_fin;
A_moving_total = A_fins_total;
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


%% 3. INITIALIZE STATES

% Flight
V = zeros(1,n);
h = zeros(1,n);
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

    % Atmospheric density
    rho(i) = rho0*exp(-h(i)/H);

    % Dynamic pressure
    q_dyn(i) = 0.5*rho(i)*V(i)^2;

    % Propulsion
    thrust = Isp*g*mdot*(m_fuel(i)>0);

    % Aerodynamic drag
    drag = 0.5*rho(i)*Cd0*A_ref*V(i)*abs(V(i));

    % Vertical flight dynamics
    a = (thrust-drag)/m_total(i)-g;

    V(i+1) = V(i)+a*dt;
    h(i+1) = max(0,h(i)+V(i+1)*dt);

    if h(i+1)==0 && V(i+1)<0
        V(i+1)=0;
    end

    % Sensor measurements including bending
    theta_raw = theta(i)+q_flex(i);
    omega_raw = omega(i)+qdot_flex(i);

    % Notch-filtered pitch
    theta_m = bN(1)*theta_raw ...
        +bN(2)*angle_x(1)+bN(3)*angle_x(2) ...
        -aN(2)*angle_y(1)-aN(3)*angle_y(2);

    % Notch-filtered pitch rate
    omega_m = bN(1)*omega_raw ...
        +bN(2)*rate_x(1)+bN(3)*rate_x(2) ...
        -aN(2)*rate_y(1)-aN(3)*rate_y(2);

    % Update filter histories
    angle_x = [theta_raw angle_x(1)];
    angle_y = [theta_m angle_y(1)];

    rate_x = [omega_raw rate_x(1)];
    rate_y = [omega_m rate_y(1)];

    % PD controller
    cmd = -Kp_fin*theta_m-Kd_fin*omega_m;

    cmd = max(-delta_max,min(delta_max,cmd));

    % Preliminary torque-limited fin command
    M_scale = q_dyn(i)*A_fin ...
        *abs(hinge_arm)*CL_max;

    if M_scale>M_servo_available

        CL_limit = M_servo_available/M_scale;

        alpha_limit = (CL_max/CL_alpha) ...
            *atanh(CL_limit);

        cmd = max(-theta(i)-alpha_limit, ...
            min(-theta(i)+alpha_limit,cmd));

    end

    % Mechanical deflection limit
    cmd = max(-delta_max,min(delta_max,cmd));

    % Servo dynamics
    delta(i+1) = delta(i) ...
        +(cmd-delta(i))*dt/tau_s;

    % Nonlinear all-moving fin coefficient
    CL = CL_max*tanh(CL_alpha*delta(i)/CL_max);

    % Active fin aerodynamic force
    F_fin = q_dyn(i)*S_fin*CL;

    % Active fin pitch moment
    torque_fin(i) = F_fin*l_fin(i);

    % Passive aerodynamic restoring moment
    if V(i)>0

        % Approximate angle of attack
        alpha_i = theta(i);

        F_normal = q_dyn(i)*A_ref ...
            *CN_alpha*alpha_i;

        torque_passive(i) = ...
            -F_normal*(x_CP-x_CoM(i));

        % Aerodynamic pitch-rate damping
        arm = x_CP-x_CoM(i);

        Kq_i = q_dyn(i)*A_ref ...
            *CN_alpha*arm^2/max(V(i),1);

        torque_passive(i) = ...
            torque_passive(i)-Kq_i*omega(i);

    end

    % Total pitch moment
    torque_total(i) = ...
        torque_fin(i)+torque_passive(i);

    % Rigid-body rotation
    pitch_accel = torque_total(i)/I_total(i);

    omega(i+1) = omega(i)+pitch_accel*dt;
    theta(i+1) = theta(i)+omega(i+1)*dt;

    % Flexible-body dynamics
    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;

end

rho(end) = rho0*exp(-h(end)/H);
q_dyn(end) = 0.5*rho(end)*V(end)^2;


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


%% 7. FROZEN-TIME BODE ANALYSIS AT MAX Q

s = tf('s');

I0 = I_total(idx_maxQ);

% Active control effectiveness
Ka = maxQ*S_fin*CL_alpha*l_fin(idx_maxQ);

% Passive restoring stiffness
Kpass = maxQ*A_ref*CN_alpha ...
    *(x_CP-x_CoM(idx_maxQ));

% Controller
C = Kp_fin+Kd_fin*s;

% Servo
Gs = 1/(tau_s*s+1);

% Aerodynamic pitch-rate damping
V0 = max(V(idx_maxQ),1);
arm0 = x_CP-x_CoM(idx_maxQ);

Kq0 = maxQ*A_ref*CN_alpha*arm0^2/V0;

% Rigid-body transfer function
Gr = 1/(I0*s^2+Kq0*s+Kpass);

% Flexible sensor transfer function
Gf = K_flex*wn_flex^2 ...
    /(s^2+2*zeta_flex*wn_flex*s+wn_flex^2);

N = N_cont;

L_rigid = C*Gs*Ka*Gr*N;
L_flexible = C*Gs*Ka*(Gr+Gf)*N;

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

    L_i = C*Gs*Ka_i*(Gr_i+Gf)*N;

    T_i = feedback(L_i,1);

    maxPole = max(real(pole(T_i)));

    fprintf('t=%.1f s | Max pole=%+.4f | Stable=%d\n', ...
        tc,maxPole,isstable(T_i));

end


%% 9. PROJECTED AREA

% Vertical-ascent AoA approximation
alpha_AoA = theta;

% Body projection
A_projected = A_ref*abs(cos(alpha_AoA)) ...
    +L_body*D_body*abs(sin(alpha_AoA));

% Four-fin projection
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

% Simplified fin incidence
alpha_fin = theta+delta;

% Nonlinear fin coefficient
CL_hinge = CL_max ...
    *tanh(CL_alpha*alpha_fin/CL_max);

% Aerodynamic hinge moment per fin
torque_hinge = q_dyn.*A_fin ...
    .*CL_hinge*hinge_arm;

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

% Max Q hinge torque stress test
test_angles = deg2rad([1 3 5]);

CL_test = CL_max ...
    *tanh(CL_alpha*test_angles/CL_max);

hinge_torque_test = maxQ*A_fin ...
    *abs(hinge_arm).*abs(CL_test);

torque_margin = M_servo_available-hinge_torque_test;

fprintf('\nMAX Q ACTUATOR TORQUE MARGIN\n')
fprintf('Available torque: %.2f N m\n',M_servo_available)

for k = 1:length(test_angles)

    fprintf(['Incidence: %.0f deg | Required: %.3f N m' ...
        ' | Margin: %+.3f N m\n'], ...
        rad2deg(test_angles(k)), ...
        hinge_torque_test(k),torque_margin(k));

end
