% Verified
% Quaternion to Euler-angle transformation. 
function thetas = quat_to_eul(qs)
    % Initialize Euler angles array. 
    thetas = zeros(size(qs(:,1:3)));

    % separate components. 
    q1=qs(:,1); q2=qs(:,2); q3=qs(:,3); q4=qs(:,4);

    % Calculate relevant cosine matrix elements. 
    c11 = 1 - 2*(q2^2 + q3^2);
    c12 = 2 * (q1*q2 + q3*q4);
    c13 = 2 * (q1*q3 - q2*q4);
    c23 = 2 * (q2*q3 + q1*q4);
    c33 = 1 - 2*(q1^2 + q2^2);

    % Calculate Euler angles. 
    thetas(:,1) = atan2(c23, c33);
    thetas(:,2) = -asin(c13);
    thetas(:,3) = atan2(c12, c11);
end