clear; clc; close all

% Adjustable parameters
m_dry = 12;       % Dry mass (kg)
m_fuel0 = 5;      % Initial propellant mass (kg)
mdot = 0.5;       % Propellant consumption (kg/s)

% Mass evolution
t = 0:0.01:15;
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


