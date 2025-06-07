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