%% Initializes
clc
clear all;
close all;

%% Environment
mu = 3.986004418e5;     % km^3/s^-2
semiMajor = 700 + 6378; % km

timeStep = 0.01;         % Seconds
timeFinal = 500;       % Seconds

meanMot = mean_mot(mu, semiMajor);
J = [124.531, 0, 0;
     0, 124.586, 0;
     0, 0, 0.704];

times = 0:timeStep:timeFinal;
numRows = size(times);

%% System Description
syms q1 q2 q3 q4 w1 w2 w3 T1 T2 T3 w1bias w2bias w3bias

% State change function. 
fVec = [0.5 * (q2*(w3) - q3*(w2) + q4*(w1));
         0.5 * (-q1*(w3 ) + q3*(w1) + q4*(w2 ));
         0.5 * (q1*(w2 ) - q2*(w1) + q4*(w3 ));
         0.5 * (-q1*(w1) - q2*(w2 ) - q3*(w3 ));
         1/J(1,1) * (3*meanMot^2 * (-J(2,2)*2*(q2*q3 + q1*q4)*(1 - 2*(q1^2 + q2^2)) + J(3,3)*(1 - 2*(q1^2 + q2^2))*2*(q2*q3 + q1*q4)) + T1 - (-J(2,2)*(w1)*(w3 ) + J(3,3)*(w3 )*(w2 )));
         1/J(2,2) * (3*meanMot^2 * (J(1,1)*2*(q1*q3 + q2*q4)*(1 - 2*(q1^2 + q2^2))) + J(3,3)*(2*(q2*q3 + q1*q4)*2*(q1*q3 + q2*q4)) + T2 - (J(1,1)*(w2 )*(w3 ) - J(3,3)*(w3 )*(w1 )));
         1/J(3,3) * (3*meanMot^2 * (-J(1,1)*2*(q1*q3 + q2*q4)*2*(q2*q3 + q1*q4) + J(2,2)*2*(q2*q3 + q1*q4)*2*(q1*q3 + q2*q4)) + T3 - (-J(1,1)*(w1 )*(w2 ) + J(2,2)*(w2 )*(w1 )));];
fFunc = matlabFunction(fVec, 'Vars', {q1, q2, q3, q4, w1, w2, w3, T1, T2, T3});

% Observation vector. 
hVec = [ atan2( 2*(q4*q1+q2*q3), 1-2*(q1^2+q2^2) );
           asin(  2*(q2*q4-q1*q3) );
          atan2( 2*(q1*q2+q3*q4), 1-2*(q2^2+q3^2) );
          w1; w2; w3];
hFunc = matlabFunction(hVec, 'Vars', {q1, q2, q3, q4, w1, w2, w3, ...
    w1bias, w2bias, w3bias});


%% Error values. 
% Euler angle std deviations.  
theta1Noise = 0.1*pi/180; theta2Noise = 0.1*pi/180; 
theta3Noise = 0.1*pi/180;

theta1Noises = theta1Noise * randn(numRows(2),1);
theta2Noises = theta2Noise * randn(numRows(2),1);
theta3Noises = theta3Noise * randn(numRows(2),1);
thetaNoise = [theta1Noises, theta2Noises, theta3Noises];

% Angular velocity bias. 
wBias = [0.1, -0.1, 0.15] * pi/180;


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

% Applies omega bias. 
wInitMeasure = w(1:3) + wBias;

%% State array initialization. 
realX = zeros(numRows(2), 7);            % Real states. 
realE = zeros(numRows(2), 3);                   % Real Euler Angles. 
measuredZ = zeros(numRows(2), 6);        % Measured states. 
measuredQ = zeros(numRows(2), 4);
filteredX = zeros(numRows(2), 7);        % Estimated states.

predictedXArr = zeros(numRows(2),7);       % Predicted states. 

realX(1,1:4) = qInit;
realX(1,5:7) = w(1:3);
realE(1,:) = [theta1, theta2, theta3];

measuredZ(1,1:3) = quat_to_eul(qInitMeasure);
measuredZ(1,4:6) = wInitMeasure;
measuredQ(1,:) = qInitMeasure;

filteredX(1,:) = realX(1,:);

%% Controller Values. 
% Initializes target quaternion. 
qTarget = eul_to_quat(0, 0, 0);

% Defines controller gains. 
K0 = 1.6;
K = 5.6;
K1 = K; K2 = K; K3 = K;

%% Filter Matrices. 
FMatx = jacobian(fVec, [q1, q2, q3, q4, w1, w2, w3]);
FMatxFunc = matlabFunction(FMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3});
HMatx = jacobian(hVec, [q1, q2, q3, q4, w1, w2, w3]);
HMatxFunc = matlabFunction(HMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3});

RMatx = diag([theta1Noise^2, theta2Noise^2, theta3Noise^2,...
              1e-6^2, 1e-6^2, 1e-6^2]);         % 6x6

G = eye(7,7);
QMatx = (1e-8)^2 * eye(7);


% Initial covariance matrix P. 
Pk1k1 = 1*eye(7);

for i = 1:numRows(2)-1
    %% Control logic. 
    % Updates control torque. 
    qError = quat_error(realX(i,1:4), qTarget);
     
    t1 = -(K0 * (qError(1) ) + K1 * realX(i,5));
    t2 = -(K0 * (qError(2) ) + K2 * realX(i,6));
    t3 = -(K0 * (qError(3) ) + K3 * realX(i,7));

    controlTorque = [t1, t2, t3];

    %% Real state evolution.
    q1 = realX(i,1); q2 = realX(i,2); q3 = realX(i,3); q4 = realX(i,4);
    w1 = realX(i,5); w2 = realX(i,6); w3 = realX(i,7);

    realX(i+1,:) = realX(i,:) + fFunc(q1, q2, q3, q4, w1, w2, w3, t1, t2, t3)'*timeStep;
    realX(i+1,1:4) = realX(i+1,1:4)/norm(realX(i+1,1:4));
    realE(i+1,:) = quat_to_eul(realX(i+1,1:4));

    %% Update measured values. 
    
    measuredZ(i+1,1:3) = realE(i+1,:) + thetaNoise(i+1,:);
    measuredZ(i+1,4:6) = realX(i+1,5:7) + wBias;

    measuredQ(i+1,:) = eul_to_quat(measuredZ(i+1,1), measuredZ(i+1,2), measuredZ(i+1,3));

    %% Predicted state. 
    q1 = filteredX(i,1); q2 = filteredX(i,2); q3 = filteredX(i,3); q4 = filteredX(i,4);
    w1 = filteredX(i,5); w2 = filteredX(i,6); w3 = filteredX(i,7);
    

    predictedX = filteredX(i,:) + timeStep * fFunc(q1, q2, q3, q4, ...
        w1, w2, w3, t1, t2, t3)';
    predictedX(1:4) = predictedX(1:4) / norm(predictedX(1:4));

    % FMatx using previous best values.
    FMatxNum = FMatxFunc(q1, q2, q3, q4, w1, w2, w3);

    % Predicted Euler angles. 
    predictedE = quat_to_eul(predictedX(1:4));

    % Predicted measurement vector. 
    zPredicted = [predictedE.';predictedX(5:7).'];

    % disp(measuredZ(i+1,:))
    % disp(zPredicted.')

    %% Filter
    q1 = measuredQ(i+1,1); q2 = measuredQ(i+1,2); q3 = measuredQ(i+1,3); q4 = measuredQ(i+1,4);
    w1 = measuredZ(i+1,4); w2 = measuredZ(i+1,5); w3 = measuredZ(i+1,6);
    
    % HMatx using predicted values. 
    HMatxNum = HMatxFunc(q1, q2, q3, q4, w1, w2, w3);
    % Discretization of matrices.
    phiMatx = exp(FMatxNum* timeStep);

    % State prediction covariance. 
    Pk1k = phiMatx*Pk1k1*phiMatx.';

    % Kalman gain matrix.
    B = HMatxNum*Pk1k*HMatxNum.';
    C = RMatx;
    K = Pk1k*HMatxNum.' / (B + C);

    % Filtered best estimate state. 
    filteredX(i+1,:) = (predictedX.' + K*(measuredZ(i+1,:).' - zPredicted)).';
    filteredX(i+1, 1:4) = filteredX(i+1, 1:4) / norm(filteredX(i+1, 1:4));
    
    % Updates state covariance matrix. 
    Pk1k1 = (eye(7)-K*HMatxNum) * Pk1k * (eye(7)-K*HMatxNum).' + ...
            K*RMatx*K.';
end

figure
subplot(2,2,1)
hold all
plot( times, realX(:, 1),"b")
plot( times, measuredQ(:, 1),"r--")
plot( times, filteredX(:, 1), "g.-")
legend({"Real q1.", "Measured q1.", "Predicted q1."})
%ylim([-0.1, 0.1])
subplot(2,2,2)
hold all
plot( times, realX(:, 2), "b")
plot( times, measuredQ(:, 2),"r--")
plot( times, filteredX(:, 2), "g.-")
legend({"Real q2.", "Measured q2.", "Predicted q2."})
%ylim([-0.1, 0.1])
subplot(2,2,3)
hold all
plot( times, realX(:, 3), "b")
plot( times, measuredQ(:, 3),"r--")
plot( times, filteredX(:, 3), "g.-")
legend({"Real q3.", "Measured q3.", "Predicted q3."})
%ylim([-0.1, 0.1])
subplot(2,2,4)
hold all
plot( times, realX(:, 4), "b")
plot( times, measuredQ(:, 4),"r--")
plot( times, filteredX(:, 4), "g.-")
legend({"Real q4.", "Measured q4.", "Predicted q4."})
%ylim([0.9, 1.1])


% figure
% subplot(2,2,1)
% hold all
% plot( times, realX(:, 1),"b")
% plot( times, filteredX(:, 1), "g.-")
% legend({"Real q1.", "Predicted q1."})
% %ylim([-0.1, 0.1])
% subplot(2,2,2)
% hold all
% plot( times, realX(:, 2), "b")
% plot( times, filteredX(:, 2), "g.-")
% legend({"Real q2.",  "Predicted q2."})
% %ylim([-0.1, 0.1])
% subplot(2,2,3)
% hold all
% plot( times, realX(:, 3), "b")
% plot( times, filteredX(:, 3), "g.-")
% legend({"Real q3.", "Predicted q3."})
% %ylim([-0.1, 0.1])
% subplot(2,2,4)
% hold all
% plot( times, realX(:, 4), "b")
% plot( times, filteredX(:, 4), "g.-")
% legend({"Real q4.", "Predicted q4."})
% %ylim([0.9, 1.1])
% 
% figure
% subplot(3,1,1)
% hold all
% plot(times, filteredX(:,5), "b")
% plot(times, measuredZ(:,4), "--r")
% legend("Estimated omega1.", "Measured omega1.")
% subplot(3,1,2)
% hold all
% plot(times, filteredX(:,6), "b")
% plot(times, measuredZ(:,5), "--r")
% legend("Estimated omega2.", "Measured omega2.")
% subplot(3,1,3)
% hold all
% plot(times, filteredX(:,7), "b")
% plot(times, measuredZ(:,6), "--r")
% legend("Estimated omega3.", "Measured omega3.")


% figure
% subplot(3,1,1)
% hold all
% plot(times, realE(:,1), "b")
% plot(times, measuredZ(:,1), "--r")
% legend("Real theta1.", "Measured theta1.")
% subplot(3,1,2)
% hold all
% plot(times, realE(:,2), "b")
% plot(times, measuredZ(:,2), "--r")
% legend("Real theta2.", "Measured theta2.")
% subplot(3,1,3)
% hold all
% plot(times, realE(:,3), "b")
% plot(times, measuredZ(:,3), "--r")
% legend("Real theta3.", "Measured theta3.")