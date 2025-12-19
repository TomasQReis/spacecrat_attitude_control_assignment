%% Initializes
clear all;
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
fMatx = [0.5 * (q2*(w3+w3bias) - q3*(w2+w2bias) + q4*(w1+w1bias));
         0.5 * (-q1*(w3+w3bias) + q3*(w1+w1bias) + q4*(w2+w2bias));
         0.5 * (q1*(w2+w2bias) - q2*(w1+w1bias) + q4*(w3+w3bias));
         0.5 * (-q1*(w1+w1bias) - q2*(w2+w2bias) - q3*(w3+w3bias));
         0 ; 0 ; 0];
fFunc = matlabFunction(fMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3, ...
    w1bias, w2bias, w3bias});

% Observation vector. 
hMatx = [atan2(2*(q4*q1 + q2*q3) ,1-2*(q1^2 + q2^2));
         asin(2*(q2*q4 - q1*q3));
         atan2(2*(q1*q2 + q3*q4) ,1-2*(q2^2 + q3^2));
         w1; w2; w3];
hFunc = matlabFunction(hMatx, 'Vars', {q1, q2, q3, q4, w1, w2, w3});

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
w = w(1:3);

wDot = w_dot(J, meanMot, qInit, w, [0,0,0], [0,0,0]);

%% Initial Measurements. 
% Creates initial error quaternion and applies. 
qInitError = eul_to_quat(theta1Noises(1), theta2Noises(1), ...
                         theta3Noises(1));
qInitMeasure = quat_mul(qInitError, qInit);
qInitMeasure = qInitMeasure/norm(qInitMeasure);

% Applies omega bias. 
wInitMeasure = w + wBias;

%% State array initialization. 
realX = zeros(numRows(2), 7);            % Real states. 
realE = zeros(numRows(2), 3);                   % Real Euler Angles. 
measuredZ = zeros(numRows(2), 6);        % Measured states. 
measuredQ = zeros(numRows(2), 4);
filteredX = zeros(numRows(2), 7);        % Estimated states.

torques = zeros(numRows(2),3);

realX(1,1:4) = qInit;
realE(1,:) = [theta1, theta2, theta3];

measuredZ(1,1:3) = [theta1 + thetaNoise(1,1), theta2 + thetaNoise(1,2),...
    theta3 + thetaNoise(1,3)];
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

%% Filter Initialization. 

for i = 1:numRows(2)-1
    %% Control logic. 
    % Updates control torque. 
    qError = quat_error(realX(i,1:4), qTarget);
     
    t1 = -(K0 * (qError(1) ) + K1 * w(1));
    t2 = -(K0 * (qError(2) ) + K2 * w(2));
    t3 = -(K0 * (qError(3) ) + K3 * w(3));

    controlTorque = [t1, t2, t3];
    torques(i,:) = controlTorque;

    %% Real state evolution.
    q1 = realX(i,1); q2 = realX(i,2); q3 = realX(i,3); q4 = realX(i,4);
    w1 = w(1); w2 = w(2); w3 = w(3);

    realX(i+1,:) = realX(i,:) + fFunc(q1, q2, q3, q4, w1, w2, w3, 0, 0, 0)'...
                                *timeStep;
    realX(i+1,1:4) = realX(i+1,1:4)/norm(realX(i+1,1:4));
    realE(i+1,:) = quat_to_eul(realX(i+1,1:4));

    % Updates w and wDot. 
    w = w + wDot * timeStep;
    wDot = w_dot(J, meanMot, realX(i+1,1:4), w, controlTorque, [0,0,0]);

    %% Update measured values. 
    measuredZ(i+1,1:3) = quat_to_eul(realX(i+1,1:4)) + thetaNoise(i+1,:);
    measuredZ(i+1,4:6) = w + wBias;

    measuredQ(i+1,:) = eul_to_quat(measuredZ(i+1,1), measuredZ(i+1,2), ...
                            measuredZ(i+1,3));
end





% figure
% subplot(2,2,1)
% hold all
% plot( times, realX(:, 1),"b")
% plot( times, measuredQ(:, 1),"r--")
% legend({"Real q1.", "Measured q1."})
% %ylim([-0.1, 0.1])
% subplot(2,2,2)
% hold all
% plot( times, realX(:, 2), "b")
% plot( times, measuredQ(:, 2),"r--")
% legend({"Real q2.", "Measured q2."})
% %ylim([-0.1, 0.1])
% subplot(2,2,3)
% hold all
% plot( times, realX(:, 3), "b")
% plot( times, measuredQ(:, 3),"r--")
% legend({"Real q3.", "Measured q3."})
% %ylim([-0.1, 0.1])
% subplot(2,2,4)
% hold all
% plot( times, realX(:, 4), "b")
% plot( times, measuredQ(:, 4),"r--")
% legend({"Real q4.", "Measured q4."})
% %ylim([0.9, 1.1])

target = zeros(numRows(2));
figure
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