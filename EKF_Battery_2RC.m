function [SoC_est, SoH_est, T_est, V1_est, V2_est] = EKF_Battery_2RC(V_meas, I_meas, T_meas, dt, C_rated)
% EKF_Battery_2RC: EKF estimation of SoC, SoH, temperature, and 2 RC voltages
%
% Inputs:
%   V_meas, I_meas, T_meas : measured voltage, current, temperature
%   dt                      : sampling time [s]
%   C_rated                 : rated battery capacity [Ah]
%
% Outputs:
%   SoC_est, SoH_est, T_est : estimated SoC, SoH, temperature
%   V1_est, V2_est          : estimated RC voltages

    N = length(V_meas);
    
    % Battery parameters
    R0_ref = 0.01; 
    R1 = 0.015; C1 = 2000; 
    R2 = 0.02;  C2 = 1500;
    
    % EKF initialization
    x = [0.9; 0; 0; 0; T_meas(1)]; % [SoC; V1; V2; SoH; T]
    P = diag([0.01 0.01 0.01 0.001 0.5]);
    Q = diag([1e-6 1e-5 1e-5 1e-8 0.01]);
    Rv = 0.05; Rt = 0.5;
    
    SoC_est = zeros(1,N); SoH_est = zeros(1,N); T_est = zeros(1,N);
    V1_est = zeros(1,N); V2_est = zeros(1,N);
    
    for k = 2:N
        I_k = I_meas(k-1);
        
        % Temperature-dependent R0
        R0 = R0_ref*(1 + 0.005*(x(5)-25));
        C_current = C_rated * x(4);
        
        % --- Prediction ---
        x_pred = zeros(5,1);
        x_pred(1) = x(1) - dt*I_k/(C_current*3600);          % SoC
        x_pred(2) = x(2) + dt*(-x(2)/(R1*C1) + I_k/C1);     % V1
        x_pred(3) = x(3) + dt*(-x(3)/(R2*C2) + I_k/C2);     % V2
        x_pred(4) = x(4);                                    % SoH
        x_pred(5) = x(5) + dt*(I_k^2*R0/50 + (25 - x(5))/100); % T
        
        % Jacobian F
        F = eye(5);
        F(1,4) = dt*I_k/(3600*C_rated*x(4)^2);
        F(2,2) = 1 - dt/(R1*C1);
        F(3,3) = 1 - dt/(R2*C2);
        F(5,5) = 1 - dt/100;
        
        P_pred = F*P*F' + Q;
        
        % --- Measurement Update ---
        z = [V_meas(k); T_meas(k)];
        h = [OCV_nonlinear(x_pred(1)) - x_pred(2) - x_pred(3) - I_k*R0;
             x_pred(5)];
        
        H = zeros(2,5);
        H(1,1) = OCV_nonlinear_dOCVdx(x_pred(1));
        H(1,2) = -1; H(1,3) = -1; H(1,4) = -0.005*I_k;
        H(2,5) = 1;
        
        R = diag([Rv^2, Rt^2]);
        y = z - h;
        S = H*P_pred*H' + R;
        K = P_pred*H'/S;
        x = x_pred + K*y;
        P = (eye(5) - K*H)*P_pred;
        
        % Save estimates
        SoC_est(k) = x(1);
        V1_est(k) = x(2);
        V2_est(k) = x(3);
        SoH_est(k) = x(4);
        T_est(k) = x(5);
    end
end

%% Nonlinear OCV
function v = OCV_nonlinear(soc)
    v = 3 + 0.5*soc + 0.7*soc.^2;
end

%% Derivative of nonlinear OCV
function dv = OCV_nonlinear_dOCVdx(soc)
    dv = 0.5 + 1.4*soc;
end
