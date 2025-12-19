clc
clear all
close all

% ------------------------------------------------
% EULER ANGLE PD CONTROLLER - SIMULATION
% ------------------------------------------------

function M = M_matrix(theta_1, theta_2, theta_3)
    % Inverse kinematics omega -> theta_dot
    M = [cos(theta_2), sin(theta_1)*sin(theta_2), cos(theta_1)*sin(theta_2);
        0, cos(theta_1)*cos(theta_2), -sin(theta_1)*cos(theta_2);
        0, sin(theta_1), cos(theta_1)];
end 



% ------------------------------------------------

% MMOI of the spacecraft
J11 = 124.531; % kg m^2
J22 = 124.586; % kg m^2
J33 = 0.704; % kg m^2

% MMOI matrix
J = [124.531 0 0;
    0 124.586 0;
    0 0 0.704];

% Disturbance torque
Td = [0.0001;
      0.0001;
      0.0001]; % Nm

% mu, radius of the orbit, ...
% orbital parameter
mu  = 398600000000000;      % m^3/s^2
R = 700*10^3 + 6378*10^3;   % m 
n   = sqrt(mu/(R * R * R)); % rad/s

% Initial attitude and rate
theta_1 = 5*pi/180;
theta_2 = 5*pi/180;
theta_3 = 5*pi/180;

theta_1_dot = 0;
theta_2_dot = 0;
theta_3_dot = 0;

P = [1 0 -sin(theta_2);
     0 cos(theta_1) sin(theta_1)*cos(theta_2);
     0 -sin(theta_1) cos(theta_1)*cos(theta_2)];

theta_d =[theta_1_dot;
          theta_2_dot;
          theta_3_dot];

Q = [cos(theta_2)*sin(theta_3);
    sin(theta_1)*sin(theta_2)*sin(theta_3)+cos(theta_1)*cos(theta_3);
    cos(theta_1)*sin(theta_2)*sin(theta_3)-sin(theta_1)*cos(theta_3)];

% Starting omegas. 
omega = P*theta_d-n*Q;

% Coefficients for dynamics
k1 = (J22-J33)/J11;
k2 = (J33-J11)/J22;
k3 = (J11-J22)/J33;

% Simulation parameters. 
t       = 0;
dt      = 0.1;
N       = 15000;
th      = zeros(3,N);
om      = zeros(3,N);

thetaNoisy = zeros(3,N);
omegaNoisy = zeros(3,N);

stateFiltered = zeros(9,N);
stateFiltered(:,1) = [theta_1; theta_2; theta_3; omega; 0; 0; 0];

% For plotting.
theta_desired = zeros(1, N);

% For noise.
theta1Noise = 0.1*pi/180 * randn(1,N);
theta2Noise = 0.1*pi/180 * randn(1,N);
theta3Noise = 0.1*pi/180 * randn(1,N);
omega1Bias = 0.1 * pi/180; 
omega2Bias = -0.1* pi/180; 
omega3Bias = 0.15* pi/180;

% For the KF. 
syms t1 t2 t3 w1 w2 w3 w1b w2b w3b To1 To2 To3
xDotVector = [1/cos(t2) * ((w1 - w1b)*cos(t2) + (w2 - w2b)*(sin(t1)*sin(t2)) + (w3 - w3b)*(cos(t1)*sin(t1)));
              1/cos(t2) * ((w2 - w2b)*(cos(t1)*cos(t2)) - (w3 - w3b)*sin(t1)*cos(t2));
              1/cos(t2) * ((w2 - w2b)*sin(t1) + (w3 - w3b)*cos(t1));
              (To1 + Td(1,1)/J11 + k1*(w2 - w2b)*(w3 - w3b));
              (To2 + Td(2,1)/J22 + k2*(w1 - w1b)*(w3 - w3b));
              (To3 + Td(3,1)/J33 + k3*(w2 - w2b)*(w1 - w1b));
              0;0;0];
xDotFunc = matlabFunction(xDotVector, 'Vars', {t1,t2,t3,w1,w2,w3,w1b,w2b,w3b,To1,To2,To3});
FMatx = jacobian(xDotVector, [t1,t2,t3,w1,w2,w3,w1b,w2b,w3b]);
FMatxFunc = matlabFunction(FMatx, 'Vars', {t1,t2,t3,w1,w2,w3,w1b,w2b,w3b});

hVector = [t1; t2; t3; w1+w1b; w2+w2b; w3+w3b];
HMatx = jacobian(hVector, [t1,t2,t3,w1,w2,w3,w1b,w2b,w3b]);
HMatxFunc = matlabFunction(HMatx, 'Vars', {t1,t2,t3,w1,w2,w3,w1b,w2b,w3b});

Pk1k1 = eye(9);
RMatx = diag([(0.1*pi/180)^2, (0.1*pi/180)^2, (0.1*pi/180)^2, 1e-6^2, 1e-6^2, 1e-6^2]);

% ------------------------------------------------
% SIMULATION
% ------------------------------------------------
for i=1:N
    % Reference control angles
    theta_des = 0;
    theta_desired(i) = theta_des;

    % ROLL CONTROL
    Kp1 = 0.4;
    Kd1 = 10;

    % PITCH CONTROL
    Kp2 = 0.225;
    Kd2 = 9;

    % YAW CONTROL
    Kp3 = 0.15;
    Kd3 = 3.5;

    % Control torque T_c
    Tc = [Kp1*(theta_des-theta_1) - Kd1*theta_d(1, 1);
          Kp2*(theta_des-theta_2) - Kd2*theta_d(2, 1);
          Kp3*(theta_des-theta_3) - Kd3*theta_d(3, 1)];

    % Dynamics
    omega_1_dot = (Tc(1, 1)+Td(1,1))/J11 + k1*omega(2, 1)*omega(3, 1);
    omega_2_dot = (Tc(2, 1)+Td(2,1))/J22 + k2*omega(1, 1)*omega(3, 1);
    omega_3_dot = (Tc(3, 1)+Td(3,1))/J33 + k3*omega(2, 1)*omega(1, 1);
    
    % omega dot
    omega_d = [omega_1_dot;
        omega_2_dot;
        omega_3_dot];

    % Integrating to get omega
    omega = omega + omega_d*dt;

    % Inverse kinematics omega -> theta_dot
    M = M_matrix(theta_1, theta_2, theta_3);

    % theta_dot update
    theta_dot = (1/cos(theta_2))*M*omega; 
    
    % Integrating to get theta
    theta_1 = theta_1 + theta_dot(1, 1)*dt;
    theta_2 = theta_2 + theta_dot(2, 1)*dt;
    theta_3 = theta_3 + theta_dot(3, 1)*dt;

    % Storing updated omega, theta & theta_dot
    om(1, i) = omega(1);
    om(2, i) = omega(2);
    om(3, i) = omega(3);

    th(1,i) = theta_1;
    th(2,i) = theta_2;
    th(3,i) = theta_3;

    theta_d(1, 1) = theta_dot(1, 1);
    theta_d(2, 1) = theta_dot(2, 1);
    theta_d(3, 1) = theta_dot(3, 1);

    % Storing noisy values. 
    omegaNoisy(1,i) = omega(1) + omega1Bias;
    omegaNoisy(2,i) = omega(2) + omega2Bias;
    omegaNoisy(3,i) = omega(3) + omega3Bias;

    thetaNoisy(1,i) = theta_1 + theta1Noise(i);
    thetaNoisy(2,i) = theta_2 + theta2Noise(i);
    thetaNoisy(3,i) = theta_3 + theta3Noise(i);

    zMeasured = [thetaNoisy(:,i); omegaNoisy(:,i)];

    %% State Prediction 
    % Dynamics

    omega_1_dot_pred = (Tc(1, 1)+Td(1,1))/J11 + ...
        k1*(stateFiltered(5,i))*(stateFiltered(6, i));
    omega_2_dot_pred = (Tc(2, 1)+Td(2,1))/J22 + ...
        k2*(stateFiltered(4,i))*(stateFiltered(6, i));
    omega_3_dot_pred = (Tc(3, 1)+Td(3,1))/J33 + ...
        k3*(stateFiltered(5,i))*(stateFiltered(4, i));
    % omega dot
    omega_d_predicted = [omega_1_dot_pred;
        omega_2_dot_pred;
        omega_3_dot_pred];

    % Integrating to get omega
    predictedState(4:6,1) = stateFiltered(4:6,i) + omega_d_predicted*dt;

    % Inverse kinematics omega -> theta_dot
    M = M_matrix(stateFiltered(1,i), stateFiltered(2,i), stateFiltered(3,i));

    % theta_dot update
    theta_dot_predict = (1/cos(stateFiltered(2,i)))*M*predictedState(4:6,1); 

    % Integrating to get theta
    predictedState(1,1) = stateFiltered(1,i) + theta_dot_predict(1, 1)*dt;
    predictedState(2,1) = stateFiltered(2,i) + theta_dot_predict(2, 1)*dt;
    predictedState(3,1) = stateFiltered(3,i) + theta_dot_predict(3, 1)*dt;

    % Applying bias. 
    predictedState(7:9,1) = stateFiltered(7:9,i);

    %% Filtering
    % Numerical Fx matrix using previous best state estimate. 
    t1 = stateFiltered(1,i); t2 = stateFiltered(2,i); t3 = stateFiltered(3,i); 
    w1 = stateFiltered(4,i); w2 = stateFiltered(5,i); w3 = stateFiltered(6,i);
    w1b = stateFiltered(7,i); w2b = stateFiltered(8,i); w3b = stateFiltered(9,i);

    FMatxNum = FMatxFunc(t1,t2,t3,w1+w1b,w2+w2b,w3+w3b,w1b,w2b,w3b);

    zPredicted = [t1; t2; t3; w1 + w1b; w2 + w2b; w3 + w3b];

    % Numerical Hx matrix and measurement prediction using predicted state.
    t1 = zMeasured(1); t2 = zMeasured(2); t3 = zMeasured(3); 
    w1 = zMeasured(4); w2 = zMeasured(5); w3 = zMeasured(6);
    w1b = predictedState(7,1); w2b = predictedState(8,1); w3b = predictedState(9,1);

    HMatxNum = HMatxFunc(t1,t2,t3,w1,w2,w3,w1b,w2b,w3b);
    
    % Discretization of Fx Matrix.
    phiMatx = eye(9) + FMatxNum * dt;

    % State prediction covariance. 
    Pk1k = phiMatx*Pk1k1*phiMatx.';

    % Kalman gain matrix.
    B = HMatxNum*Pk1k*HMatxNum.';
    C = RMatx;
    K = Pk1k*HMatxNum.' / (B + C);

    % Filtered best estimate state. 
    stateFiltered(:,i+1) = (predictedState + K*(zMeasured - zPredicted));

    % Updates state covariance matrix. 
    Pk1k1 = (eye(9)-K*HMatxNum) * Pk1k * (eye(9)-K*HMatxNum).' + ...
            K*RMatx*K.';
    
end

figure
% Attitude plots
subplot(2,3,1)
plot(theta_desired, 'Color', [1.0, 0.1, 0.6], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold on;
plot(th(1,:),'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'Actual');
plot(thetaNoisy(1,:), "b--", "displayName", "Noisy");
plot(stateFiltered(1,:), "r--", "displayName", "Predicted");
legend();
ylabel('\theta_1 [rads]')
xlabel('Steps')
title('Roll Angle \theta_1 Response')
grid on;

subplot(2,3,2)
plot(theta_desired, 'Color', [1.0, 0.1, 0.6], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold on;
plot(th(2,:),'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'Actual');
plot(thetaNoisy(2,:), "b--", "displayName", "Noisy");
plot(stateFiltered(2,:), "r--", "displayName", "Predicted");
legend();
ylabel('\theta_2 [rads]')
xlabel('Steps')
title('Pitch Angle \theta_2 Response')
grid on;

subplot(2,3,3)
plot(theta_desired, 'Color', [1.0, 0.1, 0.6], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold on;
plot(th(3,:),'Color', [0 0 0], 'LineWidth', 1.5, 'DisplayName', 'Actual');
plot(thetaNoisy(3,:), "b--", "displayName", "Noisy");
plot(stateFiltered(3,:), "r--", "displayName", "Predicted");
legend();
xlabel('Steps')
ylabel('\theta_3 [rads]')
title('Yaw Angle \theta_3 Response')
grid on;

% omega plots
subplot(2,3,4)
hold on
plot(om(1,:), 'Color', [0 0 0], 'LineWidth', 1.5);
plot(omegaNoisy(1,:), "b--", "displayName", "Noisy");
plot(stateFiltered(4,:), "r--", "displayName", "Predicted");
ylim([-0.25 0.25]);
xlabel('Steps')
ylabel('\omega_1 [rads/s]');
title('\omega_1 Response')
grid on;

subplot(2,3,5)
hold on
plot(om(2,:), 'Color', [0 0 0], 'LineWidth', 1.5);
plot(omegaNoisy(2,:), "b--", "displayName", "Noisy");
plot(stateFiltered(5,:), "r--", "displayName", "Predicted");
ylim([-0.25 0.25]);
xlabel('Steps')
ylabel('\omega_2 [rads/s]')
title('\omega_2 Response')
grid on;

subplot(2,3,6)
hold on
plot(om(3,:), 'Color', [0 0 0], 'LineWidth', 1.5);
plot(omegaNoisy(3,:), "b--", "displayName", "Noisy");
plot(stateFiltered(6,:), "r--", "displayName", "Predicted");
ylim([-0.25 0.25]);
xlabel('Steps')
ylabel('\omega_3 [rads/s]')
title('\omega_3 Response')
grid on;
