%% Initializes
clc
clear;
close all;

%% Environment
mu = 3.986004418e5;     % km^3/s^-2
semiMajor = 700 + 6378; % km

timeStep = 0.01;         % Seconds
timeFinal = 10;       % Seconds

meanMot = mean_mot(mu, semiMajor);
J = [124.531, 0, 0;
     0, 124.586, 0;
     0, 0, 0.704];

times = 0:timeStep:timeFinal;
numRows = size(times);

%% System Descriptions
syms xReal [7 1] real
syms trqs [3 1] real

% Real state change function.
% Used for propagating system dynamics (Outside of control) 
function xDotReal = fRealFunc(~, xReal, trqs, J, meanMot)
    xDotReal = [0.5 * (xReal(2)*(xReal(7)) - xReal(3)*(xReal(6)) + xReal(4)*(xReal(5)));
         0.5 * (-xReal(1)*(xReal(7) ) + xReal(3)*(xReal(5)) + xReal(4)*(xReal(6) ));
         0.5 * (xReal(1)*(xReal(6) ) - xReal(2)*(xReal(5)) + xReal(4)*(xReal(7) ));
         0.5 * (-xReal(1)*(xReal(5)) - xReal(2)*(xReal(6) ) - xReal(3)*(xReal(7) ));
         1/J(1,1) * (3*meanMot^2 * (-J(2,2)*2*(xReal(2)*xReal(3) + xReal(1)*xReal(4))*(1 - 2*(xReal(1)^2 + xReal(2)^2)) + J(3,3)*(1 - 2*(xReal(1)^2 + xReal(2)^2))*2*(xReal(2)*xReal(3) + xReal(1)*xReal(4))) + trqs(1) - (-J(2,2)*(xReal(5))*(xReal(7) ) + J(3,3)*(xReal(7) )*(xReal(6) )));
         1/J(2,2) * (3*meanMot^2 * (J(1,1)*2*(xReal(1)*xReal(3) + xReal(2)*xReal(4))*(1 - 2*(xReal(1)^2 + xReal(2)^2))) + J(3,3)*(2*(xReal(2)*xReal(3) + xReal(1)*xReal(4))*2*(xReal(1)*xReal(3) + xReal(2)*xReal(4))) + trqs(2) - (J(1,1)*(xReal(6) )*(xReal(7) ) - J(3,3)*(xReal(7) )*(xReal(5) )));
         1/J(3,3) * (3*meanMot^2 * (-J(1,1)*2*(xReal(1)*xReal(3) + xReal(2)*xReal(4))*2*(xReal(2)*xReal(3) + xReal(1)*xReal(4)) + J(2,2)*2*(xReal(2)*xReal(3) + xReal(1)*xReal(4))*2*(xReal(1)*xReal(3) + xReal(2)*xReal(4))) + trqs(3) - (-J(1,1)*(xReal(5) )*(xReal(6) ) + J(2,2)*(xReal(6) )*(xReal(5) )));];
end

% EKF state change function.
% Used in control system and filtering. 
syms t xEKF [7 1] real 
syms wMeas [3 1] real
fEKF = [0.5 * (xEKF(2)*(wMeas(3)-xEKF(7)) - xEKF(3)*(wMeas(2)-xEKF(6)) + xEKF(4)*(wMeas(1)-xEKF(5)));...
        0.5 * (-xEKF(1)*(wMeas(3)-xEKF(7)) + xEKF(3)*(wMeas(1)-xEKF(5)) + xEKF(4)*(wMeas(2)-xEKF(6)));...
        0.5 * (xEKF(1)*(wMeas(2)-xEKF(6)) - xEKF(2)*(wMeas(1)-xEKF(5)) + xEKF(4)*(wMeas(3)-xEKF(7)));...
        0.5 * (-xEKF(1)*(wMeas(1)-xEKF(5)) - xEKF(2)*(wMeas(2)-xEKF(6)) - xEKF(3)*(wMeas(3)-xEKF(7)));...
        0; 0; 0];
% Creates a jacobian matrix function. 
FMatxEKF = matlabFunction(jacobian(fEKF, xEKF), 'Vars', {xEKF, wMeas});

% Redifining the previous function in order to be able to use matlab ode45.
function xDotEKF = fEKFFunc(~, xEKF, wMeas)
    xDotEKF = [0.5 * (xEKF(2)*(wMeas(3)-xEKF(7)) - xEKF(3)*(wMeas(2)-xEKF(6)) + xEKF(4)*(wMeas(1)-xEKF(5)));...
        0.5 * (-xEKF(1)*(wMeas(3)-xEKF(7)) + xEKF(3)*(wMeas(1)-xEKF(5)) + xEKF(4)*(wMeas(2)-xEKF(6)));...
        0.5 * (xEKF(1)*(wMeas(2)-xEKF(6)) - xEKF(2)*(wMeas(1)-xEKF(5)) + xEKF(4)*(wMeas(3)-xEKF(7)));...
        0.5 * (-xEKF(1)*(wMeas(1)-xEKF(5)) - xEKF(2)*(wMeas(2)-xEKF(6)) - xEKF(3)*(wMeas(3)-xEKF(7)));...
        0; 0; 0];
end

% Observation vector. 
hEKF =  [atan2( 2*(xEKF(4)*xEKF(1)+xEKF(2)*xEKF(3)), 1-2*(xEKF(1)^2+xEKF(2)^2) );
         asin(  2*(xEKF(2)*xEKF(4)-xEKF(1)*xEKF(3)) );
         atan2( 2*(xEKF(1)*xEKF(2)+xEKF(3)*xEKF(4)), 1-2*(xEKF(2)^2+xEKF(3)^2) )];
% Observation vector function. 
function zVec = hEKFFunc (xEKF)
    zVec =  [atan2( 2*(xEKF(4)*xEKF(1)+xEKF(2)*xEKF(3)), 1-2*(xEKF(1)^2+xEKF(2)^2) );
         asin(  2*(xEKF(2)*xEKF(4)-xEKF(1)*xEKF(3)) );
         atan2( 2*(xEKF(1)*xEKF(2)+xEKF(3)*xEKF(4)), 1-2*(xEKF(2)^2+xEKF(3)^2) )];
end
% Jacobian matrix of observation vector. 
HMatxEKF = matlabFunction(jacobian(hEKF, xEKF), 'Vars', {xEKF});
disp(jacobian(hEKF, xEKF))

%% Sensor characteristics. 
% Euler angle std deviations.  
thetaNoiseSTDev = deg2rad(0.1) * ones([1,3]);

% Quaternion std deviations for initialization of P matrix. 
quatNoiseSTDev = eul_to_quat(thetaNoiseSTDev(1), thetaNoiseSTDev(2), thetaNoiseSTDev(3));

theta1Noises = thetaNoiseSTDev(1) * randn(numRows(2),1);
theta2Noises = thetaNoiseSTDev(2) * randn(numRows(2),1);
theta3Noises = thetaNoiseSTDev(3) * randn(numRows(2),1);
thetaNoise = [theta1Noises, theta2Noises, theta3Noises];

% Angular velocity bias. 
wBias = deg2rad([0.1, -0.1, 0.15]);

% Initializes EKF R matrix. 
RMatx = diag(thetaNoiseSTDev);

% Initializes EKF P matrix. (P00)
P_k_k = diag([quatNoiseSTDev, 0, 0, 0]);

%% Initial Conditions %% 
% Initial euler angles. 
theta1 = 5*pi/180; theta2 = 5*pi/180; theta3 = 5*pi/180;

% Creates angular velocity vector of LVLH frame wrt inertial (J2000).
% Defined within Body frame. 
wOrbital = - meanMot * [...
    cos(theta2)*sin(theta3),...
    sin(theta1)*sin(theta2)*sin(theta3) + cos(theta1)*cos(theta3), ... 
    cos(theta1)*sin(theta2)*sin(theta3) - sin(theta1)*cos(theta3), ...
    0];

% Assumes initial qDot=0. 
qDot = [0, 0, 0, 0];
% Converts initial euler attitude angles to quaternion.
qInit = eul_to_quat(theta1, theta2, theta3);

% Obtains initial omega vector. 
Q = [qInit(4), -qInit(3), qInit(2), qInit(1);
     qInit(3), qInit(4), -qInit(1), qInit(2);
     -qInit(2), -qInit(1), qInit(4), qInit(3);
     -qInit(1), -qInit(2), -qInit(3), qInit(4)];

w = ((2 *Q.' *qDot.') + wOrbital.').';

wDot = w_dot(J, meanMot, qInit, w, [0,0,0], [0,0,0]);

%% Initial Measurements. 
% Creates initial error quaternion and applies. 
qInitError = eul_to_quat(theta1Noises(1), theta2Noises(1), ...
                         theta3Noises(1));
qInitMeasure = quat_mul(qInitError, qInit);

%% State array initialization. 
stateReal       = zeros(7, numRows(2));                   % Real states. 
zEKF            = zeros(3, numRows(2));                    % Measured states. 
xEKFBest        = zeros(7, numRows(2));                % Estimated states.

xEKFPredicted = zeros(7, numRows(2));            % Predicted states. 

% Initial values
stateReal(1:4, 1) = qInit.';
stateReal(5:7, 1) = w(1:3).';

% Initial estimate for optimal EKFiltered state. 
startQuaternion = [1,2,3,4] / norm([1,2,3,4]);
xEKFBest(:,1) = [qInit,0,0,0];

zEKF(1:3, 1) = quat_to_eul(qInitMeasure);

% Initial control torque. 
ctrlTorque = zeros(1, 3);

%% Input array for ode45.
timespan = [0 timeStep];

for i = 1:numRows(2)-1
    disp(i)
    %% Propagation of real state.
    % Uses torques from previous step to propagate to next.  ---//---
    [~, xRealOdeOutput] =  ode45(@(t,xReal) fRealFunc(t, xReal, ctrlTorque, J, meanMot), ...
        timespan, [stateReal(:,i)]);
    stateReal(:, i+1) = xRealOdeOutput(end,:);

    % Apply bias and noise to measurements.   ---//---
    % wMeasured not saved since only required for input. 
    wMeasured = stateReal(5:7, i) + wBias.';
    % Noise added to Euler angles. Requires use of observation func. 
    zEKF(:,i+1) = hEKFFunc(stateReal(1:4,i+1)).' + thetaNoise(i+1,:);
    disp(zEKF(:,i+1))

    %% EKF Filtering of next state. 
 
    % Prediction of next state. ---//---
    [~, xEKFOdeOutput] =  ode45(@(t,xEKF) fEKFFunc(t,xEKF, wMeasured), timespan, [xEKFBest(:,i)]);
    


    xEKFPredicted(:,i+1) = xEKFOdeOutput(end,:);
    % Maintains quaternion magnitude one. 
    xEKFPredicted(1:4,i+1) = xEKFPredicted(1:4,i+1) / norm(xEKFPredicted(1:4,i+1));
    % Enforces predicted quaternions to one hemisphere. 
    if xEKFPredicted(4,i+1) < 0 
        xEKFPredicted(4,i+1) = - xEKFPredicted(4,i+1);
    end

    disp("Prediction of next state. ---//---")
    disp(xEKFPredicted(:,i+1))

    % Update covariance matrix. ---//---
    % Calculate jacobian of step. 
    FMatxEKF_i = FMatxEKF(xEKFBest(:, i), wMeasured);
    % Calculate phi matrix. 
    phiMatx = exp(timeStep * FMatxEKF_i);
    % Update covariance matrix.
    P_k1_k = phiMatx * P_k_k * phiMatx.';



    % Update Kalman gain matrix. ---//---
    disp("Update Kalman gain matrix. ---//---")
    disp(HMatxEKF(xEKFPredicted(1:4, i+1)))
    disp(P_k1_k)
    % H matrix for current step. 
    H_k1 = HMatxEKF(xEKFPredicted(1:4, i+1));
    % Update Kalman gain matrix. 
    K_k1 = P_k1_k*H_k1.' / (H_k1*P_k1_k*H_k1.' + RMatx);

    

    % Update state estimate. ---//---
    disp("Update state estimate. ---//---")
    disp(K_k1)
    disp(zEKF(:,i+1))
    disp(hEKFFunc(xEKFPredicted(1:4, i+1)))
    % Calculate innovation vector. 
    xEKFInnovation = K_k1*(zEKF(:,i+1) - hEKFFunc(xEKFPredicted(1:4, i+1)));
    % Applies angular rate bias directly. 
    xEKFBest(5:7, i+1) = xEKFPredicted(5:7,i+1) + xEKFInnovation(5:7); 
    % Applies quaterternion correction via multiplicative EKF method. 
    % (Via quaternion multiplication)
    qPredicted = quaternion(xEKFPredicted(4,i+1), xEKFPredicted(1,i+1), xEKFPredicted(2,i+1), xEKFPredicted(3,i+1));
    qInnovation = quaternion(xEKFInnovation(4), xEKFInnovation(1), xEKFInnovation(2), xEKFInnovation(3));
    qBest = compact(qPredicted * qInnovation);
    % Necessary since course notation doesn't match matlab... 
    xEKFBest(1:3, i+1) = qBest(2:4);
    xEKFBest(4, i+1) = qBest(1);


    % Maintains magnitude 1 for quaternion elements. 
    xEKFBest(1:4, i+1) = xEKFBest(1:4, i+1) / norm(xEKFBest(1:4, i+1));

    % Update error covariance matrix. ---//---
    P_k_k = (eye(7) - K_k1*H_k1)*P_k1_k;

end

figure
subplot(2,2,1)
hold on

plot( times, stateReal(1,:),"r")
plot( times, xEKFBest(1,:), "g.-")
legend({"Real q1.", "EKFiltered q1."})

subplot(2,2,2)
hold on
plot( times, stateReal(2,:),"r")
plot( times, xEKFBest(2,:), "g.-")
legend({"Real q1.", "EKFiltered q2."})

subplot(2,2,3)
hold on
plot( times, stateReal(3,:),"r")
plot( times, xEKFBest(3,:), "g.-")
legend({"Real q1.", "EKFiltered q3."})

subplot(2,2,4)
hold on
plot( times, stateReal(4,:),"r")
plot( times, xEKFBest(4,:), "g.-")
legend({"Real q1.", "EKFiltered q4."})