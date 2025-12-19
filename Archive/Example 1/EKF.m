% EKF.m
% Extended Kalman Filter (EKF)

clear; clc; close all;

% Read data file

load data1

Ts          = dt;
x_k_1k_1    = stdx_0;           % x(0|0)=E{x_0}
P_k_1k_1    = P_0;              % P(0|0)=P(0)

% Extended Kalman Filter (EKF)

ti          = 0;
tf          = Ts;          
n           = length(x_k_1k_1); % n: state dimension

for k=1:N

    % Prediction
   [t,x]            = ode45(@func1,[ti tf],x_k_1k_1);   % One stage ahead prediction using integration
   x_kk_1           = x(length(t),:)'; 	                % x(k+1|k) (prediction)
   x_pred(:,k)      = x_kk_1;
   z_kk_1           = func2(x_kk_1);                    % z(k+1|k) (predicted output)
   z_pred(:,k)      = z_kk_1;

   [Phi,Gamma]      = func3(x_k_1k_1,t,n,Ts);          	% Phi(k+1,k), Gamma(k+1,k)
   P_kk_1           = Phi*P_k_1k_1*Phi'+Gamma*Q*Gamma'; % P(k+1|k) (prediction covariance matrix)
   P_pred(:,k)      = diag(P_kk_1);
   stdx_pred(:,k)   = sqrt(diag(P_kk_1));

   % Correction
   H                = func4(x_kk_1,n);                  % Jacobian H
   Ve               = (H*P_kk_1*H'+R);                  % Pz(k+1|k) (covariance matrix of innovation)
   Vek(:,k)         = sqrt(diag(Ve));
   K                = P_kk_1*H'*inv(Ve);                % K(k+1) (gain)

   disp(size(K))
   
   x_kk             = x_kk_1+K*(z_k(:,k)-z_kk_1);       % x(k|k) (correction)
   x_cor(:,k)       = x_kk;
   
   P_kk             = (eye(n)-K*H)*P_kk_1*(eye(n)-K*H)'+K*R*K';	% P(k|k) (correction)
   P_cor(:,k)       = diag(P_kk);                          	
   stdx_cor(:,k)    = sqrt(diag(P_kk));                 

   % Next step
   x_k_1k_1         = x_kk;
   P_k_1k_1         = P_kk; 
   ti               = tf;
   tf               = tf+Ts;
end

% Plotting
figure
plot(x_k)
title('Real state')
legend('Real state')
xlabel('Samples');

figure
plot(z_k)
title('Measurement')
legend('Measurement')
xlabel('Samples');

figure
plot(x_k,'b')
hold on
plot(x_pred,'r')
title('Real and estimated states')
legend({'Real state','Estimated state'});
xlabel('Samples');

figure
plot(x_pred-x_k)
title('Estimation error')
legend('Error between real and estimated state')
xlabel('Samples');

figure
plot(x_pred-x_k)
hold on
plot(stdx_pred,'r')
hold on
plot(-stdx_pred,'g')
title('Estimation error and STD')
legend({'Error between real and estimated state','Predicted STD bounds','Predicted STD bounds'})
xlabel('Samples');

