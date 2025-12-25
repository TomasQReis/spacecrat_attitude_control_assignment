%% Initializes
clear;
close all;

%% Environment
mu = 3.986004418e5;     % km^3/s^-2
semiMajor = 700 + 6378; % km

timeStep = 0.1;         % Seconds
timeFinal = 1000;       % Seconds

meanMot = mean_mot(mu, semiMajor);
J = [124.531, 0, 0;
     0, 124.586, 0;
     0, 0, 0.704];

times = 0:timeStep:timeFinal;
numRows = size(times);

%% System Description
syms q1 q2 q3 q4 w1 w2 w3 T1 T2 T3 w1bias w2bias w3bias

% State change function. 
fVec = [0.5 * (q2*(w3-w3bias ) - q3*(w2-w2bias ) + q4*(w1-w1bias));
         0.5 * (-q1*(w3-w3bias ) + q3*(w1-w1bias) + q4*(w2-w2bias ));
         0.5 * (q1*(w2-w2bias ) - q2*(w1-w1bias) + q4*(w3-w3bias ));
         0.5 * (-q1*(w1-w1bias) - q2*(w2-w2bias ) - q3*(w3-w3bias ));
         1/J(1,1) * (3*meanMot^2 * (-J(2,2)*2*(q2*q3 + q1*q4)*(1 - 2*(q1^2 + q2^2)) + J(3,3)*(1 - 2*(q1^2 + q2^2))*2*(q2*q3 + q1*q4)) + T1 - (-J(2,2)*(w1-w1bias)*(w3-w3bias ) + J(3,3)*(w3-w3bias )*(w2-w2bias )));
         1/J(2,2) * (3*meanMot^2 * (J(1,1)*2*(q1*q3 + q2*q4)*(1 - 2*(q1^2 + q2^2))) + J(3,3)*(2*(q2*q3 + q1*q4)*2*(q1*q3 + q2*q4)) + T2 - (J(1,1)*(w2-w2bias )*(w3-w3bias ) - J(3,3)*(w3-w3bias )*(w1-w1bias )));
         1/J(3,3) * (3*meanMot^2 * (-J(1,1)*2*(q1*q3 + q2*q4)*2*(q2*q3 + q1*q4) + J(2,2)*2*(q2*q3 + q1*q4)*2*(q1*q3 + q2*q4)) + T3 - (-J(1,1)*(w1-w1bias )*(w2-w2bias ) + J(2,2)*(w2-w2bias )*(w1-w1bias )));
         0; 0; 0];
fFunc = matlabFunction(fVec, 'Vars', {q1, q2, q3, q4, w1, w2, w3, T1, T2, T3, w1bias, w2bias, w3bias});

% Observation vector. 
hVec = [ atan2( 2*(q4*q1+q2*q3), 1-2*(q1^2+q2^2) );
           asin(  2*(q2*q4-q1*q3) );
          atan2( 2*(q1*q2+q3*q4), 1-2*(q2^2+q3^2) );
          w1+w1bias; w2+w2bias; w3+w3bias];
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
% Pre-allocates saved values.
realX = zeros(numRows(2), 10);            % Real states. 
realE = zeros(numRows(2), 3);                   % Real Euler Angles. 
measuredZ = zeros(numRows(2), 6);        % Measured states. 
measuredQ = zeros(numRows(2), 4);
filteredX = zeros(numRows(2), 10);        % Estimated states.

torques = zeros(numRows(2),3);

predictedXArr = zeros(numRows(2),10);       % Predicted states. 

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
FMatx = jacobian(fVec, [q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias]);
FMatxFunc = matlabFunction(FMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias});
HMatx = jacobian(hVec, [q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias]);
HMatxFunc = matlabFunction(HMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias});

RMatx = diag([theta1Noise^2, theta2Noise^2, theta3Noise^2,...
              1e-6^2, 1e-6^2, 1e-6^2]);         % 6x6

G = eye(10,10);
QMatx = 1.5e-6^2 * eye(10);

% Initial covariance matrix P. 
Pk1k1 = 1*eye(10);

for i = 1:numRows(2)-1
    %% Control logic. 
    % Change inputs between realX, measuredZ and measuredQ and filteredX 
    % for the different controller types. I had separate files but they are
    % now broken for some reason...

    % Updates control torque. 
    qError = quat_error(realX(i,1:4), qTarget);
     
    t1 = -(K0 * (qError(1) ) + K1 * realX(i,5));
    t2 = -(K0 * (qError(2) ) + K2 * realX(i,6));
    t3 = -(K0 * (qError(3) ) + K3 * realX(i,7));

    controlTorque = [t1, t2, t3];
    torques(i,:) = controlTorque;

    %% Real state evolution.
    q1 = realX(i,1); q2 = realX(i,2); q3 = realX(i,3); q4 = realX(i,4);
    w1 = realX(i,5); w2 = realX(i,6); w3 = realX(i,7);

    realX(i+1,:) = realX(i,:) + fFunc(q1, q2, q3, q4, w1, w2, w3, t1, t2, t3, 0, 0, 0)'*timeStep;
    realX(i+1,1:4) = realX(i+1,1:4)/norm(realX(i+1,1:4));
    realE(i+1,:) = quat_to_eul(realX(i+1,1:4));

    %% Update measured values. 
    measuredZ(i+1,1:3) = realE(i+1,:) + thetaNoise(i+1,:);
    measuredZ(i+1,4:6) = realX(i,5:7) + wBias;

    measuredQ(i+1,:) = eul_to_quat(measuredZ(i+1,1), measuredZ(i+1,2), measuredZ(i+1,3));

    %% Predicted state. 
    q1 = filteredX(i,1); q2 = filteredX(i,2); q3 = filteredX(i,3); q4 = filteredX(i,4);
    w1 = filteredX(i,5); w2 = filteredX(i,6); w3 = filteredX(i,7);
    w1bias = filteredX(i,8); w2bias = filteredX(i,9); w3bias = filteredX(i,10);
    
    predictedX = filteredX(i,:) + timeStep * fFunc(q1, q2, q3, q4, ...
        w1, w2, w3, t1, t2, t3, w1bias, w2bias, w3bias)';
    predictedX(1:4) = predictedX(1:4) / norm(predictedX(1:4));

    % FMatx using previous best values.
    FMatxNum = FMatxFunc(q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias);

    % Predicted measurement vector. 
    zPredicted = hFunc(q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias);
    
    %% Filter
    q1 = measuredQ(i+1,1); q2 = measuredQ(i+1,2); q3 = measuredQ(i+1,3); q4 = measuredQ(i+1,4);
    w1 = measuredZ(i+1,4); w2 = measuredZ(i+1,5); w3 = measuredZ(i+1,6);
    
    % HMatx using predicted values. 
    HMatxNum = HMatxFunc(q1, q2, q3, q4, w1, w2, w3, w1bias, w2bias, w3bias);
    % Discretization of matrices.
    [phiMatx, Gamma] = c2d(FMatxNum, G, timeStep);

    % State prediction covariance. 
    Pk1k = phiMatx*Pk1k1*phiMatx.' + Gamma*QMatx*Gamma.';

    % Kalman gain matrix.
    B = HMatxNum*Pk1k*HMatxNum.';
    C = RMatx;
    K = Pk1k*HMatxNum.' / (B + C);

    % Filtered best estimate state. 
    filteredX(i+1,:) = (predictedX.' + K*(measuredZ(i+1,:).' - zPredicted)).';
    filteredX(i+1, 1:4) = filteredX(i+1, 1:4) / norm(filteredX(i+1, 1:4));
    
    % Updates state covariance matrix. 
    Pk1k1 = (eye(10)-K*HMatxNum) * Pk1k * (eye(10)-K*HMatxNum).' + ...
            K*RMatx*K.';
end


%% PLOTTING FOR THETA AND TORQUE.
target = zeros(size(times));
figure 
fontsize(scale=1.5)
subplot(3,1,1)
hold on
plot(times, realE(:,1), "b")
plot(times, target, "--r")
legend("Real theta1.", "Target theta1.")
ylabel("Angle [rad]")
xlabel("Time from start [s]")
subplot(3,1,2)
hold on
plot(times, realE(:,2), "b")
plot(times, target, "--r")
legend("Real theta2.", "Target theta2.")
ylabel("Angle [rad]")
xlabel("Time from start [s]")
subplot(3,1,3)
hold on
plot(times, realE(:,3), "b")
plot(times, target, "--r")
legend("Real theta3.", "Target theta3.")
ylabel("Angle [rad]")
xlabel("Time from start [s]")

figure
subplot(3,1,1)
hold on
plot(times, torques(:,1), "b")
ylabel("Torque T_1 [Nm].")
xlabel("Time from start [s]")
subplot(3,1,2)
hold on
plot(times, torques(:,2), "b")
ylabel("Torque T_2 [Nm].")
xlabel("Time from start [s]")
subplot(3,1,3)
hold on
plot(times, torques(:,3), "b")
ylabel("Torque T_3 [Nm].")
xlabel("Time from start [s]")

%% PLOTTING FOR BIAS ESTIMATE
figure
subplot(3,1,1)
hold on
plot(times, filteredX(:,8), "b")
plot(times, measuredZ(:,4), "--r")
legend("Estimated omega1 bias.", "Measured omega1.")
xlabel("Time from start [s]")
ylabel("[rad/s]")

subplot(3,1,2)
hold on
plot(times, filteredX(:,9), "b")
plot(times, measuredZ(:,5), "--r")
legend("Estimated omega2 bias.", "Measured omega2.")
xlabel("Time from start [s]")
ylabel("[rad/s]")

subplot(3,1,3)
hold on
plot(times, filteredX(:,10), "b")
plot(times, measuredZ(:,6), "--r")
legend("Estimated omega3 bias.", "Measured omega3.")
xlabel("Time from start [s]")
ylabel("[rad/s]")

