clc; close all;
clear all;
%% Simulation Parameters
N = 3600;           % 1-hour simulation
dt = 1;             % 1-second sampling
C_rated = 2.3;      % Ah

% True battery states
SoC_true = zeros(1,N); SoH_true = zeros(1,N); T_true = zeros(1,N);
V_true = zeros(1,N); 
V1 = 0;  % RC voltage

% Battery parameters
R0_ref = 0.01; R1 = 0.015; C1 = 2000;

%% Generate Variable Current Profile (simulating EV driving)
I = zeros(1,N);
for k = 1:N
    if mod(k,600) < 200          % Acceleration
        I(k) = 5 + 0.5*randn;   % 5A +/- noise
    elseif mod(k,600) < 400      % Cruising
        I(k) = 2 + 0.2*randn;   % 2A +/- noise
    else                         % Braking / regen
        I(k) = -3 + 0.3*randn;  % -3A +/- noise (charging)
    end
end

T_amb = 25;  % Ambient temperature
SoC_true(1) = 1; SoH_true(1) = 1; T_true(1) = T_amb;

%% Simulate True Battery Dynamics with Variable Load
for k = 2:N
    SoH_true(k) = SoH_true(k-1) - 5e-7;    % slow fade
    C_current = C_rated * SoH_true(k);
    R0 = R0_ref*(1 + 0.005*(T_true(k-1)-25));
    
    % SoC dynamics
    dSoC = -I(k-1)/(C_current*3600);
    SoC_true(k) = SoC_true(k-1) + dSoC*dt;
    
    % RC dynamics
    dV1 = dt*(-V1/(R1*C1) + I(k-1)/C1);
    V1 = V1 + dV1;
    
    % Terminal voltage
    V_true(k) = OCV_nonlinear(SoC_true(k)) - V1 - I(k)*R0;
    
    % Temperature dynamics
    dT = dt*(I(k-1)^2*R0/50 + (T_amb - T_true(k-1))/100);
    T_true(k) = T_true(k-1) + dT;
end

%% Add Measurement Noise
V_meas = V_true + 0.02*randn(1,N);  % 20 mV noise
I_meas = I + 0.05*randn(1,N);       % 50 mA noise
T_meas = T_true + 0.5*randn(1,N);   % 0.5°C noise

%% Run EKF Estimator
[SoC_est, SoH_est, T_est] = EKF_Battery_Noisy(V_meas, I_meas, T_meas, dt, C_rated);

%% Plot Results
figure;
subplot(4,1,1);
plot(1:N, I, 'k'); xlabel('Time (s)'); ylabel('Current (A)'); title('Battery Current Profile'); grid on;

subplot(4,1,2);
plot(1:N, SoC_true, 'b', 1:N, SoC_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('SoC'); legend('True','Estimated'); grid on;

subplot(4,1,3);
plot(1:N, SoH_true, 'b', 1:N, SoH_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('SoH'); legend('True','Estimated'); grid on;

subplot(4,1,4);
plot(1:N, T_true, 'b', 1:N, T_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Temperature (°C)'); legend('True','Estimated'); grid on;












%% Simulation Parameters
N = 3600;       % 1-hour simulation
dt = 1;         % 1-second sampling
C_rated = 2.3;  % Ah

% Initialize true battery states
SoC_true = zeros(1,N); SoH_true = zeros(1,N); T_true = zeros(1,N);
V_true = zeros(1,N); V1 = 0; V2 = 0;

% Battery parameters
R0_ref = 0.01; R1 = 0.015; C1 = 2000; R2 = 0.02; C2 = 1500;
T_amb = 25;

%% Generate Variable EV Current Profile
I = zeros(1,N);
for k = 1:N
    if mod(k,600) < 200          % Acceleration
        I(k) = 5 + 0.5*randn;   
    elseif mod(k,600) < 400      % Cruising
        I(k) = 2 + 0.2*randn;   
    else                         % Braking / regenerative
        I(k) = -3 + 0.3*randn;  
    end
end

SoC_true(1) = 1; SoH_true(1) = 1; T_true(1) = T_amb;

%% Simulate True Battery Dynamics (2 RC Branches)
for k = 2:N
    SoH_true(k) = SoH_true(k-1) - 5e-7;  % slow fade
    C_current = C_rated*SoH_true(k);
    R0 = R0_ref*(1 + 0.005*(T_true(k-1)-25));
    
    dSoC = -I(k-1)/(C_current*3600);
    SoC_true(k) = SoC_true(k-1) + dSoC*dt;
    
    dV1 = dt*(-V1/(R1*C1) + I(k-1)/C1);
    V1 = V1 + dV1;
    dV2 = dt*(-V2/(R2*C2) + I(k-1)/C2);
    V2 = V2 + dV2;
    
    V_true(k) = OCV_nonlinear(SoC_true(k)) - V1 - V2 - I(k)*R0;
    
    % Temperature dynamics
    dT = dt*(I(k-1)^2*R0/50 + (T_amb - T_true(k-1))/100);
    T_true(k) = T_true(k-1) + dT;
end

%% Add Measurement Noise
V_meas = V_true + 0.02*randn(1,N);  % voltage noise
I_meas = I + 0.05*randn(1,N);       % current noise
T_meas = T_true + 0.5*randn(1,N);   % temperature noise

%% Run 2-RC EKF
[SoC_est, SoH_est, T_est, V1_est, V2_est] = EKF_Battery_2RC(V_meas, I_meas, T_meas, dt, C_rated);

%% Plot Results
figure;
subplot(5,1,1);
plot(1:N, I, 'k'); xlabel('Time (s)'); ylabel('Current (A)'); title('Battery Current Profile'); grid on;

subplot(5,1,2);
plot(1:N, SoC_true, 'b', 1:N, SoC_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('SoC'); legend('True','Estimated'); grid on;

subplot(5,1,3);
plot(1:N, SoH_true, 'b', 1:N, SoH_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('SoH'); legend('True','Estimated'); grid on;

subplot(5,1,4);
plot(1:N, T_true, 'b', 1:N, T_est, 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Temperature (°C)'); legend('True','Estimated'); grid on;

subplot(5,1,5);
plot(1:N, V1_est, 'r', 1:N, V2_est, 'g', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('RC Voltages (V)'); legend('V1','V2'); grid on;

%% Nonlinear OCV function
function v = OCV_nonlinear(soc)
    v = 3 + 0.5*soc + 0.7*soc.^2;
end

% %% Nonlinear OCV function
% function v = OCV_nonlinear(soc)
%     v = 3 + 0.5*soc + 0.7*soc.^2;
% end