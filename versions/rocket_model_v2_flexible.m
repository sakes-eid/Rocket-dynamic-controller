clear; clc; close all

% Adjustable parameters
m_dry = 12;       % Dry mass (kg)
m_fuel0 = 5;      % Initial propellant mass (kg)
mdot = 0.5;       % Propellant consumption (kg/s)

% Mass evolution
t = 0:0.01:120;
m_fuel = max(0, m_fuel0 - mdot*t);
m_total = m_dry + m_fuel;

% Plot
plot(t,m_total,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Rocket mass (kg)')
title('Rocket Mass During Propellant Burn')
grid on



% Adjustable positions measured from rocket nose
x_dry = 1.2;       % Dry mass center (m)
x_fuel = 2.0;      % Propellant center (m)

% Center of mass
x_CoM = (m_dry*x_dry + m_fuel*x_fuel)./m_total;

figure
plot(t,x_CoM,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Center of mass from nose (m)')
title('Rocket Center of Mass During Burn')
grid on



% Adjustable dry rocket inertia
I_dry = 8;       % kg*m^2 (illustrative)

% Total pitch moment of inertia
I_total = I_dry ...
    + m_dry*(x_dry - x_CoM).^2 ...
    + m_fuel.*(x_fuel - x_CoM).^2;

figure
plot(t,I_total,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Moment of inertia (kg m^2)')
title('Rocket Moment of Inertia During Burn')
grid on



% Adjustable applied torque
torque = 1;           % N*m

% Angular acceleration from rest
alpha = torque./I_total;

figure
plot(t,alpha,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Angular acceleration (rad/s^2)')
title('Angular Acceleration vs Propellant Burn')
grid on



% Angular velocity and pitch angle
omega = cumtrapz(t,alpha);
theta = cumtrapz(t,omega);

figure
plot(t,rad2deg(theta),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Rocket Pitch Angle vs Time')
grid on



% Adjustable PD controller
Kp = 10;
Kd = 8;

theta_PD = zeros(size(t));
omega_PD = zeros(size(t));
theta_PD(1) = deg2rad(5);

dt = t(2)-t(1);

for i = 1:length(t)-1
    torque_PD = -Kp*theta_PD(i) - Kd*omega_PD(i);
    alpha_PD = torque_PD/I_total(i);

    omega_PD(i+1) = omega_PD(i) + alpha_PD*dt;
    theta_PD(i+1) = theta_PD(i) + omega_PD(i)*dt;
end

figure
plot(t,rad2deg(theta_PD),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Rocket Pitch Correction with PD Control')
grid on



% Damping ratio during propellant burn
zeta = Kd./(2*sqrt(I_total*Kp));

figure
plot(t,zeta,'LineWidth',1.5)
yline(1,'--')
xlabel('Time (s)')
ylabel('Damping ratio')
title('Pitch Controller Damping Ratio')
grid on



% Adjustable actuator response time
tau_s = 0.04;

theta_servo = zeros(size(t));
omega_servo = zeros(size(t));
torque_actual = zeros(size(t));

theta_servo(1) = deg2rad(5);

for i = 1:length(t)-1
    torque_cmd = -Kp*theta_servo(i) - Kd*omega_servo(i);

    torque_actual(i+1) = torque_actual(i) ...
        + (torque_cmd - torque_actual(i))*dt/tau_s;

    alpha_servo = torque_actual(i)/I_total(i);

    omega_servo(i+1) = omega_servo(i) + alpha_servo*dt;
    theta_servo(i+1) = theta_servo(i) + omega_servo(i)*dt;
end

figure
plot(t,rad2deg(theta_PD),t,rad2deg(theta_servo),'LineWidth',1.5)
legend('Ideal actuator','Actuator with delay')
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Effect of Actuator Dynamics')
grid on



% Adjustable structural flexibility
f_flex = 4;          % Natural frequency (Hz)
zeta_flex = 0.05;    % Damping ratio (assumed)
K_flex = 1;          % Mode gain (assumed)

% Flexible mode transfer function
s = tf('s');
wn_flex = 2*pi*f_flex;

G_flex = K_flex*wn_flex^2 / ...
    (s^2 + 2*zeta_flex*wn_flex*s + wn_flex^2);

figure
bode(G_flex)
grid on
title('Rocket Flexible Mode - Bode Diagram')



% Adjustable flexible-mode disturbance
disturbance_amplitude = deg2rad(0.1);

% Flexible vibration excited at t = 1 second
u_flex = disturbance_amplitude * double(t >= 0);

theta_flex = lsim(G_flex,u_flex,t);

% Sensor measures rigid-body pitch plus bending
theta_measured = theta_servo(:) + theta_flex(:);

figure
plot(t,rad2deg(theta_servo), ...
    t,rad2deg(theta_measured),'LineWidth',1.5)
legend('Rigid-body pitch','Sensor measurement')
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Effect of Structural Flexibility on Pitch Measurement')
grid on



% Adjustable flight parameters
Isp = 180;       % Specific impulse (s)
CdA = 0.015;     % Drag coefficient * reference area (m^2)
rho0 = 1.225;    % Sea-level air density (kg/m^3)
H = 8500;        % Atmospheric scale height (m)
g = 9.81;        % Gravity (m/s^2)

% Vertical flight
V = zeros(size(t));
h = zeros(size(t));

for i = 1:length(t)-1
    rho = rho0*exp(-h(i)/H);
    thrust = Isp*g*mdot*(m_fuel(i)>0);
    drag = 0.5*rho*CdA*V(i)*abs(V(i));

    a = (thrust-drag)/m_total(i)-g;

    V(i+1) = V(i)+a*dt;
h(i+1) = max(0,h(i)+V(i+1)*dt);

if h(i+1) == 0 && V(i+1) < 0
    V(i+1) = 0;
end
end

figure
plot(t,V,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Vertical velocity (m/s)')
title('Rocket Velocity During Flight')
grid on

figure
plot(t,h,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('Rocket Altitude During Flight')
grid on



% Dynamic pressure
rho = rho0*exp(-h/H);
q_dyn = 0.5*rho.*V.^2;

figure
plot(t,q_dyn/1000,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Dynamic pressure (kPa)')
title('Rocket Dynamic Pressure During Flight')
grid on



% Adjustable fin aerodynamic parameters
S_fin = 0.01;            % Effective fin area (m^2)
CL_alpha = 3;            % Lift slope (per radian)
CL_max = 1.0;            % Assumed coefficient limit
delta = deg2rad(5);      % Test fin deflection (rad)

% Nonlinear aerodynamic coefficient
CL = CL_max*tanh(CL_alpha*delta/CL_max);

% Aerodynamic fin force
F_fin = q_dyn*S_fin*CL;

figure
plot(t,F_fin,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Fin aerodynamic force (N)')
title('Fin Force at 5 Degree Deflection')
grid on



% Adjustable fin aerodynamic-center location
x_fin = 2.5;        % Distance from nose (m)

% Pitch moment arm
l_fin = x_fin - x_CoM;

% Aerodynamic pitch torque
torque_fin = F_fin.*l_fin;

figure
plot(t,torque_fin,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Pitch torque (N m)')
title('Aerodynamic Pitch Torque from Fin Deflection')
grid on



% Adjustable fin controller parameters
Kp_fin = 0.5;
Kd_fin = 0.2;
delta_max = deg2rad(5);

% Fin angle command using previous pitch response
delta_cmd = -Kp_fin*theta_servo - Kd_fin*omega_servo;

% Limit commanded fin deflection
delta_cmd = max(-delta_max,min(delta_max,delta_cmd));

figure
plot(t,rad2deg(delta_cmd),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Commanded fin angle (degrees)')
title('PD Controller Fin Deflection Command')
grid on



% Simulated actual fin angle
delta_actual = zeros(size(t));

for i = 1:length(t)-1
    delta_actual(i+1) = delta_actual(i) ...
        + (delta_cmd(i)-delta_actual(i))*dt/tau_s;
end

figure
plot(t,rad2deg(delta_cmd),t,rad2deg(delta_actual),'LineWidth',1.5)
legend('Commanded angle','Actual servo angle')
xlabel('Time (s)')
ylabel('Fin angle (degrees)')
title('Servo Response to Fin Commands')
xlim([0 15])
grid on



% Aerodynamic force from actual servo deflection
CL_actual = CL_max*tanh(CL_alpha*delta_actual/CL_max);

F_fin_actual = q_dyn.*S_fin.*CL_actual;

% Pitch torque
torque_fin_actual = F_fin_actual.*l_fin;

figure
plot(t,torque_fin_actual,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Pitch torque (N m)')
title('Pitch Torque from Actual Fin Movement')
xlim([0 15])
grid on


delta_history = zeros(size(t));
torque_history = zeros(size(t));


% Closed-loop aerodynamic pitch control
theta_aero = zeros(size(t));
omega_aero = zeros(size(t));
delta_aero = zeros(size(t));

theta_aero(1) = deg2rad(5);

for i = 1:length(t)-1
    % PD controller
    cmd = -Kp_fin*theta_aero(i) - Kd_fin*omega_aero(i);
    cmd = max(-delta_max,min(delta_max,cmd));

    % Servo dynamics
    delta_aero(i+1) = delta_aero(i) ...
        + (cmd-delta_aero(i))*dt/tau_s;

    % Nonlinear fin aerodynamics
    CL_i = CL_max*tanh(CL_alpha*delta_aero(i)/CL_max);
    F_i = q_dyn(i)*S_fin*CL_i;
    delta_history(i) = delta_aero(i);
    torque_history(i) = F_i*l_fin(i);

    % Pitch dynamics
    alpha_i = F_i*l_fin(i)/I_total(i);
    omega_aero(i+1) = omega_aero(i) + alpha_i*dt;
    theta_aero(i+1) = theta_aero(i) + omega_aero(i)*dt;
end


figure
plot(t,rad2deg(delta_history))
xlim([0 5])
xlabel('Time (s)')
ylabel('Actual fin angle (degrees)')
grid on

figure
plot(t,torque_history)
xlim([0 5])
xlabel('Time (s)')
ylabel('Pitch torque (N m)')
grid on



% Adjustable fin-to-structure coupling
K_flex = 0.001;     % rad/(N*m), assumed

G_flex = K_flex*wn_flex^2 / ...
    (s^2 + 2*zeta_flex*wn_flex*s + wn_flex^2);

% Structural response to actual fin torque
theta_flex = lsim(G_flex,torque_history,t);

figure
plot(t,rad2deg(theta_flex),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Flexible angular contribution (degrees)')
title('Structural Bending Excited by Fin Torque')
xlim([0 5])
grid on



% Separate dynamic bending from quasi-static bending
theta_static = K_flex*torque_history(:);
theta_dynamic = theta_flex(:) - theta_static;

figure
plot(t,rad2deg(theta_dynamic),'LineWidth',1.5)
xlim([0 4])
xlabel('Time (s)')
ylabel('Dynamic bending component (degrees)')
title('Flexible Mode - Dynamic Response')
grid on



% Sensor response including structural bending
theta_sensor = theta_aero(:) + theta_flex(:);

figure
plot(t,rad2deg(theta_aero), ...
    t,rad2deg(theta_sensor),'LineWidth',1.5)
legend('Rigid-body pitch','Measured pitch')
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Pitch Sensor with Structural Flexibility')
xlim([0 5])
grid on



% Coupled rigid-body and flexible-body feedback
theta_c = zeros(size(t));
omega_c = zeros(size(t));
delta_c = zeros(size(t));
q_c = zeros(size(t));
qdot_c = zeros(size(t));
torque_c = zeros(size(t));

theta_c(1) = deg2rad(5);

for i = 1:length(t)-1

    % Sensor measurements
    theta_m = theta_c(i) + q_c(i);
    omega_m = omega_c(i) + qdot_c(i);

    % PD controller
    cmd = -Kp_fin*theta_m - Kd_fin*omega_m;
    cmd = max(-delta_max,min(delta_max,cmd));

    % Servo
    delta_c(i+1) = delta_c(i) ...
        + (cmd-delta_c(i))*dt/tau_s;

    % Aerodynamic fin force and torque
    CL_i = CL_max*tanh(CL_alpha*delta_c(i)/CL_max);
    F_i = q_dyn(i)*S_fin*CL_i;
    torque_c(i) = F_i*l_fin(i);

    % Rigid-body pitch dynamics
    alpha_i = torque_c(i)/I_total(i);
    omega_c(i+1) = omega_c(i) + alpha_i*dt;
    theta_c(i+1) = theta_c(i) + omega_c(i+1)*dt;

    % Flexible bending dynamics
    qdd = wn_flex^2*(K_flex*torque_c(i)-q_c(i)) ...
        - 2*zeta_flex*wn_flex*qdot_c(i);

    qdot_c(i+1) = qdot_c(i) + qdd*dt;
    q_c(i+1) = q_c(i) + qdot_c(i+1)*dt;
end

figure
plot(t,rad2deg(theta_aero), ...
    t,rad2deg(theta_c),'LineWidth',1.5)
legend('Rigid-only feedback','Flexible feedback')
xlabel('Time (s)')
ylabel('Pitch angle (degrees)')
title('Effect of Structural Feedback on Pitch Control')
xlim([0 5])
grid on


figure
plot(t,rad2deg(q_c),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Flexible sensor angle (degrees)')
title('Structural Bending with Closed-Loop Feedback')
xlim([0 5])
grid on



figure
plot(t,rad2deg(theta_flex), ...
    t,rad2deg(q_c),'LineWidth',1.5)
legend('Without structural feedback','With structural feedback')
xlabel('Time (s)')
ylabel('Flexible sensor angle (degrees)')
title('Effect of Feedback on Structural Bending')
xlim([0 5])
grid on



% Frozen-time linear model at maximum dynamic pressure
[~,idx] = max(q_dyn);

I0 = I_total(idx);
Ka = q_dyn(idx)*S_fin*CL_alpha*l_fin(idx);

% Controller and servo
C = Kp_fin + Kd_fin*s;
Gs = 1/(tau_s*s+1);

% Rigid-body and flexible sensor dynamics
Gr = 1/(I0*s^2);
Gf = K_flex*wn_flex^2 / ...
    (s^2 + 2*zeta_flex*wn_flex*s + wn_flex^2);

% Complete open-loop transfer function
L = C*Gs*Ka*(Gr+Gf);

figure
margin(L)
grid on
title('Open-Loop Stability Margins at Max Q')



L_rigid = C*Gs*Ka*Gr;

[~,PM_rigid] = margin(L_rigid)
[~,PM_flexible] = margin(L)

figure
bode(L_rigid,L)
legend('Rigid rocket','Flexible rocket')
grid on


T_rigid = feedback(L_rigid,1);
T_flexible = feedback(L,1);

pole(T_rigid)
pole(T_flexible)

isstable(T_rigid)
isstable(T_flexible)