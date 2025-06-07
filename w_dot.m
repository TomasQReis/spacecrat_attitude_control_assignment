% Verified
% Outputs angular velocity change vector. 
function wDot = w_dot(J, meanMot, q, w, controlTorque, wBias)
    % Defines cosine matrix elements via quaternion elems.
    C13 = 2*(q(1)*q(3) + q(2)*q(4));
    C23 = 2*(q(2)*q(3) + q(1)*q(4));
    C33 = 1 - 2*(q(1)^2 + q(2)^2);
   
    % Defines intermediate matrices. 
    C = [0, -C33, C23;
        C33, 0, -C13;
        -C23, C13, 0];

    w1 = w(1) + wBias(1); 
    w2 = w(2) + wBias(2);
    w3 = w(3) + wBias(3);

    W = [0, -w3, w2;
         w3, 0, -w1;
         -w2, w1, 0];

    leftSide = 3*meanMot^2*C*J*[C13;C23;33];
    rightSide = [0.001; 0.001; 0.001] - W*J*+[w1, w2, w3].';

    wDot = [inv(J)*((leftSide + rightSide + controlTorque.'))]';
end