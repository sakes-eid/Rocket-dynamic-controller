
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

% Vertical trajectory guidance
K_gamma = 0.2;
alpha_guidance_max = deg2rad(0.75);

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


% Experimental dynamic-pressure guidance schedule
alpha_lowQ_max = deg2rad(1.5);
q_schedule = 60000;   % Pa


% Mach-dependent drag coefficient
M_table = [0 0.5 0.7 0.8 0.9 1.0 ...
    1.1 1.2 1.4 1.6 2.0];

Cd_ratio = [1 1 1.03 1.10 1.25 1.55 ...
    1.70 1.60 1.35 1.25 1.20];

% Precompute smooth interpolation
Cd_curve = pchip(M_table,Cd_ratio);


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

    % Atmospheric density
    rho(i) = rho0*exp(-h(i)/H);

    % Total air-relative speed
    speed_air = hypot(vx(i),V(i));

    % Dynamic pressure
    q_dyn(i) = 0.5*rho(i)*speed_air^2;


    %% 4.2 FLIGHT-PATH ANGLE AND ANGLE OF ATTACK

    % Flight-path angle measured from vertical
    gamma_i = atan2(vx(i),V(i));

    % Actual angle of attack
    alpha_i = atan2( ...
        sin(theta(i)-gamma_i), ...
        cos(theta(i)-gamma_i));


    %% 4.3 AERODYNAMIC NORMAL FORCES

    F_passive = 0;
    F_active = 0;
    F_normal_i = 0;

    % Apply simplified normal-force model during ascent
    if speed_air > 1 && V(i) > 0

        % Nose and neutral fin contributions
        F_passive = q_dyn(i)*A_ref ...
            *CN_alpha*alpha_i;

        % Incremental force from moving the fins
        CL_i = CL_max*tanh( ...
            CL_alpha*delta(i)/CL_max);

        F_active = -q_dyn(i)*S_fin*CL_i;

        % Total aerodynamic normal force
        F_normal_i = F_passive+F_active;

    end


    %% 4.4 PROPULSION AND DRAG

    % Engine thrust
    thrust = Isp*g*mdot*(m_fuel(i)>0);

    
    % Atmospheric temperature
    T_i = max(216.65,288.15-0.0065*h(i));
    
    % Speed of sound
    a_i = sqrt(1.4*287.05*T_i);
    
    % Current Mach number
    Mach_i = speed_air/a_i;
    
    % Mach-dependent drag coefficient
    M_lookup = max(M_table(1),min(M_table(end),Mach_i));
    Cd_i = Cd0*ppval(Cd_curve,M_lookup);
    
    % Aerodynamic drag factor
    drag_factor = 0.5*rho(i)*Cd_i*A_ref*speed_air;


    %% 4.5 AERODYNAMIC FORCE DIRECTIONS

    % Unit vector perpendicular to flight velocity
    if speed_air > 1

        nx = V(i)/speed_air;
        nz = -vx(i)/speed_air;

    else

        nx = 0;
        nz = 0;

    end


    %% 4.6 TWO-DIMENSIONAL FLIGHT DYNAMICS

    % Horizontal acceleration
    ax = (thrust*sin(theta(i)) ...
        -drag_factor*vx(i) ...
        +F_normal_i*nx)/m_total(i);

    % Vertical acceleration
    az = (thrust*cos(theta(i)) ...
        -drag_factor*V(i) ...
        +F_normal_i*nz)/m_total(i)-g;

    % Update velocities
    vx(i+1) = vx(i)+ax*dt;
    V(i+1) = V(i)+az*dt;

    % Update positions
    x_pos(i+1) = x_pos(i)+vx(i+1)*dt;
    h(i+1) = max(0,h(i)+V(i+1)*dt);

    % Ground contact
    if h(i+1)==0 && V(i+1)<0

        V(i+1) = 0;
        vx(i+1) = 0;

    end


    %% 4.7 SENSOR MEASUREMENTS
    
    % Sensor measurements including structural bending
    theta_raw = theta(i) + q_flex(i);
    omega_raw = omega(i) + qdot_flex(i);


    % Notch-filtered pitch angle
    theta_m = bN(1)*theta_raw ...
        +bN(2)*angle_x(1)+bN(3)*angle_x(2) ...
        -aN(2)*angle_y(1)-aN(3)*angle_y(2);

    % Notch-filtered pitch rate
    omega_m = bN(1)*omega_raw ...
        +bN(2)*rate_x(1)+bN(3)*rate_x(2) ...
        -aN(2)*rate_y(1)-aN(3)*rate_y(2);

    % Update notch filter histories
    angle_x = [theta_raw angle_x(1)];
    angle_y = [theta_m angle_y(1)];

    rate_x = [omega_raw rate_x(1)];
    rate_y = [omega_m rate_y(1)];


    %% 4.8 PD CONTROLLER - TRUE VERTICAL ATTITUDE

    % Control during powered ascent and coast
    if V(i) > 1 && speed_air > 1

        % Desired rocket attitude relative to true vertical
        theta_ref = 0;

        % Attitude tracking error
        theta_error = atan2( ...
            sin(theta_m-theta_ref), ...
            cos(theta_m-theta_ref));

        % PD attitude controller
        cmd = -Kp_fin*theta_error-Kd_fin*omega_m;

    else

        % Neutral fins outside the controlled ascent regime
        cmd = 0;

    end


    %% 4.9 PRELIMINARY SERVO TORQUE LIMIT

    % Approximate aerodynamic hinge load scale
%    M_scale = q_dyn(i)*A_fin ...
%        *abs(hinge_arm)*CL_max;

 %   if M_scale > M_servo_available

   %     CL_limit = M_servo_available/M_scale;

%        alpha_limit = (CL_max/CL_alpha) ...
 %           *atanh(CL_limit);

        % Existing simplified incidence constraint
  %      cmd = max(-theta(i)-alpha_limit, ...
  %          min(-theta(i)+alpha_limit,cmd));

    %end

    % Enforce mechanical limit
    cmd = max(-delta_max,min(delta_max,cmd));

    % Servo response
    delta(i+1) = delta(i) ...
        +(cmd-delta(i))*dt/tau_s;


    %% 4.10 AERODYNAMIC PITCH MOMENTS

    % Active fin moment
    torque_fin(i) = -F_active*l_fin(i);

    % Passive restoring moment from actual AoA
    torque_passive(i) = ...
        -F_passive*(x_CP-x_CoM(i));

    % Aerodynamic pitch-rate damping
    if speed_air > 1 && V(i) > 0

        arm = x_CP-x_CoM(i);

        Kq_i = q_dyn(i)*A_ref ...
            *CN_alpha*arm^2/speed_air;

        torque_passive(i) = ...
            torque_passive(i)-Kq_i*omega(i);

    end

    % Total pitch moment
    torque_total(i) = ...
        torque_fin(i)+torque_passive(i);


    %% 4.11 RIGID-BODY ROTATION

    % Pitch angular acceleration
    pitch_accel = torque_total(i)/I_total(i);

    % Integrate pitch dynamics
    omega(i+1) = omega(i)+pitch_accel*dt;
    theta(i+1) = theta(i)+omega(i+1)*dt;


    %% 4.12 STRUCTURAL FLEXIBILITY

    % Flexible mode excited by active fin moment
    qdd = wn_flex^2 ...
        *(K_flex*torque_fin(i)-q_flex(i)) ...
        -2*zeta_flex*wn_flex*qdot_flex(i);

    % Integrate structural dynamics
    qdot_flex(i+1) = qdot_flex(i)+qdd*dt;
    q_flex(i+1) = q_flex(i)+qdot_flex(i+1)*dt;

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
gamma = atan2(vx,V);

% Angle of attack
alpha_AoA = atan2(sin(theta-gamma),cos(theta-gamma));

% Powered flight, excluding very low speeds
speed_air = hypot(vx,V);
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

figure
plot(t(burn_idx),rad2deg(theta(burn_idx)), ...
    t(burn_idx),rad2deg(gamma(burn_idx)), ...
    'LineWidth',1.2)

legend('Rocket pitch','Flight-path angle')
xlabel('Time (s)')
ylabel('Angle from vertical (deg)')
title('Powered Ascent: Attitude vs Trajectory')
grid on


%% GUIDANCE TRACKING DIAGNOSTIC

% Reconstruct requested AoA
alpha_requested = max(-alpha_guidance_max, ...
    min(alpha_guidance_max,-K_gamma*gamma));

% Compare requested and actual AoA during powered flight
figure
plot(t(valid),rad2deg(alpha_requested(valid)), ...
    '--','LineWidth',1.3)
hold on
plot(t(valid),rad2deg(alpha_AoA(valid)), ...
    'LineWidth',1.3)

legend('Requested AoA','Actual AoA')
xlabel('Time (s)')
ylabel('Angle of Attack (deg)')
title('Guidance Command Tracking')
xlim([0 burn_time])
grid on

% Burnout comparison
i = burn_idx(end);

fprintf('\nGUIDANCE TRACKING AT BURNOUT\n')
fprintf('Flight-path angle: %.3f deg\n',rad2deg(gamma(i)))
fprintf('Requested AoA: %.3f deg\n',rad2deg(alpha_requested(i)))
fprintf('Actual AoA: %.3f deg\n',rad2deg(alpha_AoA(i)))
fprintf('Horizontal velocity: %.3f m/s\n',vx(i))

fprintf('\nROCKET ATTITUDE RELATIVE TO VERTICAL\n')
fprintf('Initial pitch: %.3f deg\n',rad2deg(theta(1)))
fprintf('Pitch at burnout: %.3f deg\n', ...
    rad2deg(theta(burn_idx(end))))
fprintf('Pitch at apogee: %.3f deg\n', ...
    rad2deg(theta(idx_apogee)))
fprintf('Maximum powered-flight pitch: %.3f deg\n', ...
    max(abs(rad2deg(theta(burn_idx)))))


%% SENSOR TRACKING DIAGNOSTIC AT BURNOUT

i = burn_idx(end);

alpha_ref_i = max(-alpha_guidance_max, ...
    min(alpha_guidance_max,-K_gamma*gamma(i)));

alpha_sensed = atan2( ...
    sin(theta(i)+q_flex(i)-gamma(i)), ...
    cos(theta(i)+q_flex(i)-gamma(i)));

fprintf('\nBURNOUT SENSOR DIAGNOSTIC\n')
fprintf('Requested AoA: %.3f deg\n',rad2deg(alpha_ref_i))
fprintf('Actual AoA: %.3f deg\n',rad2deg(alpha_AoA(i)))
fprintf('Sensor bending: %.3f deg\n',rad2deg(q_flex(i)))
fprintf('Sensed AoA (unfiltered): %.3f deg\n',rad2deg(alpha_sensed))
fprintf('Active fin torque: %.2f N m\n',torque_fin(i))


%% POWERED FLIGHT AND COAST GUIDANCE DIAGNOSTIC

fprintf('\nASCENT GUIDANCE PERFORMANCE\n')
fprintf('Time | Pitch | Path | AoA | Vx | Q(kPa) | Downrange\n')

for tc = [5 10 15 20 25 30 33]

    [~,i] = min(abs(t-tc));

    fprintf('%4.0f | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f | %8.1f\n', ...
        t(i),rad2deg(theta(i)),rad2deg(gamma(i)), ...
        rad2deg(alpha_AoA(i)),vx(i), ...
        q_dyn(i)/1000,x_pos(i));

end

% Horizontal displacement accumulated after burnout
i_burn = burn_idx(end);

fprintf('\nDownrange at burnout: %.1f m\n',x_pos(i_burn))
fprintf('Downrange at apogee: %.1f m\n',x_pos(idx_apogee))
fprintf('Coast downrange increase: %.1f m\n', ...
    x_pos(idx_apogee)-x_pos(i_burn))


%% HORIZONTAL FORCE BALANCE DIAGNOSTIC

fprintf('\nHORIZONTAL ACCELERATION BALANCE\n')
fprintf('Time | Thrust | Drag | Aero | Net (m/s^2)\n')

for tc = [3 5 8 9.5 12 15 20]

    [~,i] = min(abs(t-tc));

    speed_i = hypot(vx(i),V(i));

    thrust_i = Isp*g*mdot*(m_fuel(i)>0);

    a_thrust = thrust_i*sin(theta(i))/m_total(i);

    a_drag = -0.5*rho(i)*Cd0*A_ref ...
        *speed_i*vx(i)/m_total(i);

    a_net = (vx(i+1)-vx(i))/dt;

    a_aero = a_net-a_thrust-a_drag;

    fprintf('%4.1f | %+6.2f | %+6.2f | %+6.2f | %+6.2f\n', ...
        t(i),a_thrust,a_drag,a_aero,a_net);

end


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


% Force on one actuated fin:
% neutral-fin AoA contribution plus control increment

F_fin_each = q_dyn .* ( ...
    A_ref*(CN_fins/n_fins).*alpha_AoA ...
    - A_fin*CL_max.*tanh(CL_alpha*delta/CL_max));

% Aerodynamic hinge torque per fin
torque_hinge = F_fin_each*hinge_arm;

% Restrict diagnostic to forward ascent
torque_hinge(V<=0) = 0;

idx_power = find(t<=burn_time & V>20);

[max_hinge,k] = max(abs(torque_hinge(idx_power)));
idx_hinge = idx_power(k);

fprintf('Peak powered-flight hinge torque: %.3f N m\n',max_hinge)
fprintf('Time: %.2f s\n',t(idx_hinge))

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



%% DESCENT ORIENTATION DIAGNOSTIC

speed = hypot(vx,V);

% +1 = nose-first, -1 = tail-first, 0 = sideways
nose_alignment = ...
    (sin(theta).*vx + cos(theta).*V) ...
    ./max(speed,1e-6);

% Exclude low-speed conditions and ground contact
descending = V < -20 & h > 0;

figure
plot(t(descending),nose_alignment(descending))
yline(1,'--','Nose-first')
yline(-1,'--','Tail-first')
yline(0,':','Sideways')
xlabel('Time (s)')
ylabel('Nose alignment')
title('Rocket Orientation During Descent')
ylim([-1.1 1.1])
grid on



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

    % Passive pitch stiffness at Max Q
    Kpass_test = maxQ*A_ref*CN_rocket ...
        *(CP_test-x_CoM(idx_maxQ));

    % Active fin control effectiveness at Max Q
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

    % Move fin aerodynamic center forward
    x_fin_test = x_fin_CP-shift;

    % Recalculate rocket center of pressure
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

% Standard atmosphere approximation
T_air = 288.15 - 0.0065*min(h,11000);

% Local speed of sound
a_sound = sqrt(1.4*287.05*T_air);

% Air-relative speed
speed_series = hypot(vx,V);

% Mach number
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

% Smooth interpolation
M_plot = linspace(0,2,500);
Cd_plot = Cd0*interp1(M_table,Cd_ratio,M_plot,'pchip');

figure
plot(M_plot,Cd_plot,'LineWidth',1.5)
xline(1,'--','Mach 1')
xlabel('Mach number')
ylabel('Drag coefficient Cd')
title('Mach-Dependent Rocket Drag Coefficient')
grid on

% Drag coefficient at current maximum Mach
Cd_maxMach = Cd0*interp1( ...
    M_table,Cd_ratio,M_max,'pchip');

fprintf('\nMACH-DEPENDENT DRAG\n')
fprintf('Subsonic baseline Cd: %.3f\n',Cd0)
fprintf('Cd at maximum Mach: %.3f\n',Cd_maxMach)
fprintf('Increase over baseline: %.1f%%\n', ...
    100*(Cd_maxMach/Cd0-1))
