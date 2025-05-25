%% Environment
mu = 3.986004418e5;     % km^3/s^-2
semiMajor = 700 + 6378; % km

meanMot = mean_mot(mu, semiMajor);

%% Initial Conditions %% 
% Creates angular velocity vector of LVLH frame wrt inertial (J2000).
% Defined within LVLH frame. 
orbitalW = [0 -meanMot 0];

% Converts initial euler attitude angles to quaternion.
theta1 = 0.0; theta2 = 0.0; theta3 = 0.0;
quatInit = eul_to_quat(theta1, theta2, theta3);



%% FUNCTIONS %%
% Verified. 
% Mean motion of orbit with given semiMajor around body mu. 
function meanMot = mean_mot(mu, semiMajor)
    meanMot = sqrt(mu / semiMajor^3);
end

% Verified. 
% Quaternion matrix multiplication function. 
% Follows the notation where q_4 is the real part of the quaternion. 
function quatMul = quat_mul(quat1, quat2)
    intermMatrix = [quat2(4), quat2(3), -quat2(2), quat2(1);
                    -quat2(3), quat2(4), quat2(1), quat2(2);
                    quat2(2), -quat2(1), quat2(4), quat2(3);
                    -quat2(1), -quat2(2), -quat2(3), quat2(4)];
    quatMul = (intermMatrix * quat1.').';
end

% Verified. 
% Euler to quaternion transformation. 
function quat = eul_to_quat(theta1, theta2, theta3)

    % Precompute half-angles. 
    t1 = theta1/2;
    t2 = theta2/2;
    t3 = theta3/2;

    % sines and cosines of half-angles. 
    c1 = cos(t1);  s1 = sin(t1);
    c2 = cos(t2);  s2 = sin(t2);
    c3 = cos(t3);  s3 = sin(t3);

    % Quaternion components. 
    q1 =  s1*c2*c3  -  c1*s2*s3;
    q2 =  c1*s2*c3  +  s1*c2*s3;
    q3 =  c1*c2*s3  -  s1*s2*c3;
    q4 =  c1*c2*c3  +  s1*s2*s3;

    % Pack into a row-vector. 
    quat = [q1, q2, q3, q4];
end


