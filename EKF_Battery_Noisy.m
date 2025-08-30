function [SoC_est, SoH_est, T_est] = EKF_Battery_Noisy(V_meas, I_meas, T_meas, dt, C_rated)
% EKF_Battery_Noisy: Estimate SoC, SoH, and battery temperature from noisy measurements
%
% Inputs:
%   V_meas   : measured terminal voltage [V] (vector)
%   I_meas   : measured battery current [A] (vector, positive = discharge)
%   T_meas   : measured battery temperature [°C] (vector)
%   dt       : sampling time [s]
%   C_rated  : rated battery capacity [Ah]
%
% Outputs:
%   SoC_est  : estimated state of charge [0-1]
%   SoH_est  : estimated state of health [0-1]
%   T_est    : estimated battery temperature [°C]

    N = length(V_meas);
    
    % Battery parameters
    R0_ref = 0.01;  % Ohm at 25°C
    R1 = 0.015;
    C1 = 2000;
    
    % EKF initialization
    x = [0.9; 0; 0.95; T_meas(1)];  % [SoC; V1; SoH; T]
    P = diag([0.01 0.01 0.001 0.5]);
    Q = diag([1e-6 1e-5 1e-8 0.01]); % process noise
    Rv = 0.05;  % measurement noise voltage
    Rt = 0.5;   % measurement noise temperature
    
    SoC_est = zeros(1,N);
    SoH_est = zeros(1,N);
    T_est = zeros(1,N);
    
    for k = 2:N
        I_k = I_meas(k-1);
        
        % Predicted battery parameters
        C_current = C_rated * x(3);
        R0 = R0_ref*(1 + 0.005*(x(4)-25));
        
        % --- Prediction ---
        x_pred = [x(1) - dt*I_k/(C_current*3600);          % SoC
                  x(2) + dt*(-x(2)/(R1*C1) + I_k/C1);    % V1
                  x(3);                                  % SoH
                  x(4) + dt*(I_k^2*R0/50 + (T_meas(1) - x(4))/100)]; % T
              
        F = [1, 0, dt*I_k/(3600*C_rated*x(3)^2), 0;
             0, 1-dt/(R1*C1), 0, 0;
             0, 0, 1, 0;
             0, 0, 0, 1-dt/100];
        P_pred = F*P*F' + Q;
        
        % --- Measurement Update ---
        % Measurement vector: voltage + temperature
        z = [V_meas(k); T_meas(k)];
        
        % Predicted measurement
        h = [OCV_nonlinear(x_pred(1)) - x_pred(2) - I_k*R0;
             x_pred(4)];
        
        % Jacobian H
        H = [OCV_nonlinear_dOCVdx(x_pred(1)), -1, 0, -0.005*I_k;
             0, 0, 0, 1];
        
        % Measurement noise covariance
        R = diag([Rv^2, Rt^2]);
        
        % Kalman Gain
        y = z - h;
        S = H*P_pred*H' + R;
        K = P_pred*H'/S;
        
        % Update state
        x = x_pred + K*y;
        P = (eye(4) - K*H)*P_pred;
        
        % Save estimates
        SoC_est(k) = x(1);
        SoH_est(k) = x(3);
        T_est(k) = x(4);
    end
end

%% Nonlinear OCV function
function v = OCV_nonlinear(soc)
    v = 3 + 0.5*soc + 0.7*soc.^2;
end

%% Derivative of nonlinear OCV
function dv = OCV_nonlinear_dOCVdx(soc)
    dv = 0.5 + 1.4*soc;
end
