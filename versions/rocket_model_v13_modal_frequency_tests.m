%% V13 CHECKPOINT: EI-based free vibration and 4 Hz forcing / frequency sweep
% Standalone follow-up to rocket_model_v12_beam_4hz_test.m.
% This is a STRUCTURAL TEST ONLY: it does not rerun or alter V11/V12 flight.
% It reproduces V12's first free-free Euler-Bernoulli mode approximation.
% No two-way aeroelastic feedback, no aerodynamic damping, no torsion.
% The older 4 Hz sensor/mount oscillator is deliberately NOT included here.
clear; clc; close all

%% USER-ADJUSTABLE TEST PARAMETERS
EI_body = 40000;          % N m^2, assumed uniform bending rigidity
zeta_body = 0.03;         % first-mode damping ratio (assumed)
F_test_N = 15;            % N, peak sinusoidal force applied at fin station
f_excitation = 4;        % Hz, imposed excitation (NOT natural frequency)
q_initial_mm = 0.10;     % mm, initial maximum-mode deformation for free decay

dt = 0.00025;            % s; resolves ~46 Hz mode at 87 steps per cycle
L_body = 2.7;            % m
m_dry = 12;              % kg, uniformly distributed approximation
m_fuel0 = 5;             % kg, point mass at x_fuel
x_fuel = 2.0;            % m aft of nose
x_sensor = 0.55;         % m aft of nose (V12 assumption)

% Same four-fin geometry and reference fin aerodynamic-center station as V12.
fin_root = 0.22;
fin_tip = 0.10;
fin_span = 0.08;
fin_tip_LE_offset = 0.06;
x_fin_root_LE = L_body-fin_root;
c_mean = (2/3)*(fin_root^2+fin_root*fin_tip+fin_tip^2) ...
    /(fin_root+fin_tip);
x_MAC_LE = x_fin_root_LE+(fin_tip_LE_offset/3) ...
    *(fin_root+2*fin_tip)/(fin_root+fin_tip);
x_fin = x_MAC_LE+0.25*c_mean;

assert(EI_body>0 && zeta_body>0 && zeta_body<1 && F_test_N>0, ...
    'EI, damping and test force must be positive; require 0 < zeta < 1');
assert(f_excitation>0 && dt>0 && q_initial_mm>=0, ...
    'Excitation frequency, timestep, initial amplitude must be valid');

%% EXACT SAME TRIAL SHAPE AS V12 (not a re-fitted mode)
B = 4.730040744862704; % first non-rigid free-free root
x = linspace(0,L_body,801);
b = B/L_body;
eta = (cosh(B)-cos(B))/(sinh(B)-sin(B));
phi_raw = cosh(b*x)+cos(b*x) ...
    -eta*(sinh(b*x)+sin(b*x));
normalizer = max(abs(phi_raw));
phi = phi_raw/normalizer;
dphi = b*(sinh(b*x)-sin(b*x) ...
    -eta*(cosh(b*x)+cos(b*x)))/normalizer;
ddphi = b^2*(cosh(b*x)-cos(b*x) ...
    -eta*(sinh(b*x)-sin(b*x)))/normalizer;

% Generalized mode q has units of meters because phi is dimensionless.
K = EI_body*trapz(x,ddphi.^2); % N/m
M_dry = (m_dry/L_body)*trapz(x,phi.^2); % kg
phi_fuel = interp1(x,phi,x_fuel,'pchip');
phi_fin = interp1(x,phi,x_fin,'pchip');
phi_tail = phi(end);
slope_sensor = interp1(x,dphi,x_sensor,'pchip');

M = [M_dry+m_fuel0*phi_fuel^2, M_dry]; % liftoff / burnout
fn = sqrt(K./M)/(2*pi);
wd = (2*pi*fn)*sqrt(1-zeta_body^2);
C = 2*zeta_body*sqrt(K.*M);
scenario_names = {'Liftoff (17 kg)', 'Burnout (12 kg)'};

%% TEST 1: EXACT FREE-DECAY SOLUTION WITH ZERO EXTERNAL FORCE
% Same INITIAL MODAL DISPLACEMENT, zero initial modal velocity, both masses.
t_free = 0:dt:0.70;
q0 = q_initial_mm/1000;
q_free = zeros(2,numel(t_free));
qd_free = zeros(2,numel(t_free));
E_free = zeros(2,numel(t_free));
free_freq_est = zeros(1,2);
for s = 1:2
    wn = 2*pi*fn(s);
    a = zeta_body*wn;
    w = wd(s);
    q_free(s,:) = q0*exp(-a*t_free).* ...
        (cos(w*t_free)+(a/w)*sin(w*t_free));
    qd_free(s,:) = -q0 .* exp(-a*t_free) .* ((wn^2/w)*sin(w*t_free));
    E_free(s,:) = 0.5*K*q_free(s,:).^2+0.5*M(s)*qd_free(s,:).^2;
    % Frequency from zero crossings of the free response (independent check).
    crossings = find(q_free(s,1:end-1).*q_free(s,2:end)<0);
    crossing_times = t_free(crossings) ...
        +dt*abs(q_free(s,crossings)) ...
        ./(abs(q_free(s,crossings))+abs(q_free(s,crossings+1)));
    free_freq_est(s) = 1/(2*mean(diff(crossing_times)));
end

%% TEST 2: HARMONIC FREQUENCY RESPONSE FOR THE SAME FIN FORCE
% M qdd + C qd + K q = phi_fin * F_test_N * sin(2*pi*f*t)
% Steady-state response: |q(f)| = |phi_fin F| / |K-M*w^2+i*C*w|.
f_sweep = 0.2:0.025:80;
w_sweep = 2*pi*f_sweep;
Q0 = abs(phi_fin*F_test_N); % generalized force amplitude, N
amp_q = zeros(2,numel(f_sweep));
amp_sensor_deg = zeros(2,numel(f_sweep));
f_peak = zeros(1,2);
q_amp_4 = zeros(1,2);
q_amp_peak = zeros(1,2);
rot_amp_4 = zeros(1,2);
for s = 1:2
    amp_q(s,:) = Q0./sqrt((K-M(s)*w_sweep.^2).^2 ...
        +(C(s)*w_sweep).^2);
    amp_sensor_deg(s,:) = rad2deg(abs(slope_sensor)*amp_q(s,:));
    q_amp_4(s) = Q0/sqrt((K-M(s)*(2*pi*f_excitation)^2)^2 ...
        +(C(s)*2*pi*f_excitation)^2);
    [q_amp_peak(s),peak_idx] = max(amp_q(s,:));
    f_peak(s) = f_sweep(peak_idx);
    rot_amp_4(s) = rad2deg(abs(slope_sensor)*q_amp_4(s));
end

%% TEST 3: FINITE 15 N FORCING, 4 Hz VERSUS BODY RESONANCE
% Both use the same liftoff mass and same force envelope. Observe the
% difference between forced vibration and the ringdown after shutoff.
t_force = 0:dt:3.60;
N_force = numel(t_force);
force_env = zeros(1,N_force);
turn_on = 0.35;
turn_off_start = 2.65;
turn_off = 3.00;
mask_up = t_force<=turn_on;
mask_flat = t_force>turn_on & t_force<turn_off_start;
mask_down = t_force>=turn_off_start & t_force<=turn_off;
force_env(mask_up) = sin(0.5*pi*t_force(mask_up)/turn_on).^2;
force_env(mask_flat) = 1;
force_env(mask_down) = sin(0.5*pi*(turn_off-t_force(mask_down)) ...
    /(turn_off-turn_off_start)).^2;

f_compare = [f_excitation, fn(1)];
q_forced = zeros(2,N_force);
qd_forced = zeros(2,N_force);
force_inputs = zeros(2,N_force);
for s = 1:2
    force_inputs(s,:) = F_test_N*force_env.*sin(2*pi*f_compare(s)*t_force);
    for j = 1:N_force-1
        qdd = (phi_fin*force_inputs(s,j)-C(1)*qd_forced(s,j) ...
            -K*q_forced(s,j))/M(1);
        qd_forced(s,j+1) = qd_forced(s,j)+qdd*dt;
        q_forced(s,j+1) = q_forced(s,j)+qd_forced(s,j+1)*dt;
    end
end
steady_interval = t_force>=1.4 & t_force<=2.5;
measured_amps_mm = max(abs(q_forced(:,steady_interval)),[],2)'*1000;
predicted_amps_mm = zeros(1,2);
for s = 1:2
    omega_f = 2*pi*f_compare(s);
    predicted_amps_mm(s) = 1000*Q0/sqrt( ...
        (K-M(1)*omega_f^2)^2+(C(1)*omega_f)^2);
end
relative_amp_err = abs(measured_amps_mm-predicted_amps_mm) ...
    ./max(predicted_amps_mm,eps);

%% NUMERICAL / PHYSICAL SANITY CHECKS
assert(all(M>0) && K>0 && all(isfinite(amp_q(:))), ...
    'Invalid modal properties or frequency response');
assert(all(abs(free_freq_est-fn*sqrt(1-zeta_body^2))./fn<0.002), ...
    'The free-decay zero-crossing frequency disagrees with the model');
assert(all(diff(E_free(1,:))<=1e-11) && all(diff(E_free(2,:))<=1e-11), ...
    'Free-decay mechanical energy should not grow for positive damping');
% The plateau comparison checks the numerical integrator against the
% closed-form forced solution, without requiring MATLAB toolboxes.
assert(all(relative_amp_err<0.10), ...
    'Forced-response amplitude disagrees with theoretical steady state >10%%');

%% PRINT CHECKPOINT
fprintf('\nV13 - BODY MODE FREQUENCY / FREE DECAY / 4 Hz SWEEP\n');
fprintf('EI = %.0f N m^2 | damping = %.1f%% | fin force = %.1f N\n', ...
    EI_body,100*zeta_body,F_test_N);
fprintf('Modal K = %.1f N/m | phi(fin)=%.4f | phi(tail)=%.4f\n', ...
    K,phi_fin,phi_tail);
fprintf('Sensor station = %.2f m | phi slope = %.5f 1/m\n', ...
    x_sensor,slope_sensor);
for s = 1:2
    fprintf('\n%s\n',scenario_names{s});
    fprintf('  Modal mass: %.4f kg\n',M(s));
    fprintf('  Natural frequency: %.3f Hz\n',fn(s));
    fprintf('  Measured free-decay frequency: %.3f Hz\n',free_freq_est(s));
    fprintf('  4 Hz tail deflection amplitude: %.5f mm\n', ...
        1000*abs(phi_tail)*q_amp_4(s));
    fprintf('  4 Hz sensor slope amplitude: %.6f deg\n',rot_amp_4(s));
    fprintf('  Response peak: %.3f Hz; tail amplitude: %.5f mm\n', ...
        f_peak(s),1000*abs(phi_tail)*q_amp_peak(s));
    fprintf('  Peak / 4 Hz amplitude ratio: %.2f x\n', ...
        q_amp_peak(s)/q_amp_4(s));
end
fprintf('\nTIME-DOMAIN CHECK (liftoff mass, same 15 N force)\n');
fprintf('  4 Hz: measured %.5f mm / predicted %.5f mm\n', ...
    measured_amps_mm(1),predicted_amps_mm(1));
fprintf('  %.2f Hz: measured %.5f mm / predicted %.5f mm\n', ...
    f_compare(2),measured_amps_mm(2),predicted_amps_mm(2));
fprintf('  Numerical amplitude agreement: %.2f%% / %.2f%%\n', ...
    100*relative_amp_err(1),100*relative_amp_err(2));
fprintf('All built-in checks passed. NO flight/aerodynamic feedback in this test.\n');

%% PLOTS (ALL DEDICATED TO THE NEW FREQUENCY CHECKPOINT)
figure('Name','V13 - Mode Shape and Frozen Natural Frequencies');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
plot(x,phi,'LineWidth',1.3); hold on
xline(x_sensor,':','Sensor');xline(x_fin,':','Fin')
xlabel('Position aft of nose (m)');ylabel('Normalized modal shape');
title('First free-free elastic mode (same assumed shape as V12)');grid on
nexttile
bar([1 2],fn);hold on
yline(f_excitation,'--','4 Hz excitation');
set(gca,'XTick',[1 2],'XTickLabel',scenario_names);
ylabel('Frequency (Hz)');title('Body frequencies from EI and modal mass');grid on

figure('Name','V13 - Free Decay at Body Natural Frequencies');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
plot(t_free,1000*q_free(1,:),t_free,1000*q_free(2,:), ...
    'LineWidth',1.1); hold on
plot(t_free, q_initial_mm*exp(-zeta_body*2*pi*fn(1)*t_free), ...
    '--','Color',[.5 .5 .5]);
plot(t_free,-q_initial_mm*exp(-zeta_body*2*pi*fn(1)*t_free), ...
    '--','Color',[.5 .5 .5]);
xlabel('Time (s)');ylabel('Modal displacement (mm)');
legend(scenario_names{1},scenario_names{2},'Liftoff envelopes', ...
    'Location','northeast');
title('Free vibration: same initial deformation, no applied force');
xlim([0 0.35]);grid on
nexttile
plot(t_free,1e6*E_free(1,:),t_free,1e6*E_free(2,:),'LineWidth',1.15)
xlabel('Time (s)');ylabel('Elastic + kinetic energy (microjoules)');
legend(scenario_names,'Location','northeast');
title('Mechanical energy decreases under positive structural damping');
grid on

figure('Name','V13 - Harmonic Frequency Sweep');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
semilogy(f_sweep,1000*abs(phi_tail)*amp_q(1,:), ...
         f_sweep,1000*abs(phi_tail)*amp_q(2,:), ...
         'LineWidth',1.25);hold on
xline(f_excitation,':','4 Hz force');
xline(fn(1),'--','Liftoff mode');xline(fn(2),'--','Burnout mode');
xlabel('Forcing frequency (Hz)');ylabel('Steady tail amplitude (mm)');
legend(scenario_names,'Location','best');
title('Same 15 N fin forcing: 4 Hz vs natural bending resonance');grid on
nexttile
semilogy(f_sweep,amp_sensor_deg(1,:), ...
         f_sweep,amp_sensor_deg(2,:),'LineWidth',1.25);hold on
xline(f_excitation,':','4 Hz force');
xline(fn(1),'--','Liftoff mode');xline(fn(2),'--','Burnout mode');
xlabel('Forcing frequency (Hz)');ylabel('Sensor-local elastic slope (deg)');
legend(scenario_names,'Location','best');
title('Physical sensor rotation from the same structural mode');grid on

figure('Name','V13 - Equal Force: 4 Hz vs Body Resonance');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
plot(t_force,1000*q_forced(1,:),t_force,1000*q_forced(2,:), ...
    'LineWidth',1.0);hold on
xline(3,':','Force off');
legend(sprintf('%.1f Hz',f_compare(1)), ...
       sprintf('%.2f Hz (body mode)',f_compare(2)),'Location','best');
ylabel('Modal displacement (mm)');
title('Finite-duration sinusoidal input: same 15 N peak force');grid on
nexttile
plot(t_force,1000*q_forced(1,:),t_force,1000*q_forced(2,:), ...
    'LineWidth',1.0);
xlabel('Time (s)');ylabel('Modal displacement (mm)');
title('Zoom to compare instantaneous oscillations');
xlim([1.4 1.75]);grid on

figure('Name','V13 - Bending Stiffness Sensitivity');
EI_values = linspace(5000,80000,160);
f_lift_EI = fn(1)*sqrt(EI_values/EI_body);
f_burn_EI = fn(2)*sqrt(EI_values/EI_body);
plot(EI_values,f_lift_EI,EI_values,f_burn_EI,'LineWidth',1.25);hold on
xline(EI_body,':','Selected EI');yline(f_excitation,'--','4 Hz input');
legend(scenario_names,'Location','best');
xlabel('Assumed bending rigidity EI (N m^2)');
ylabel('First-mode frequency (Hz)');
title('Stiffness and natural frequency are NOT independent for a given mode');
grid on
