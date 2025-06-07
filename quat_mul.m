% Verified. 
% Quaternion matrix multiplication function. 
% Follows the notation where q_4 is the real part of the quaternion. 
function qMul = quat_mul(q1, q2)
    % quat_mul.
    %   Returns the quaternion product q = q1 ∘ q2,
    %   where each input is a 4-element row [qx, qy, qz, q0], and the output
    %   is in the same format.
    intermMatrix = [q2(4), q2(3), -q2(2), q2(1);
                    -q2(3), q2(4), q2(1), q2(2);
                    q2(2), -q2(1), q2(4), q2(3);
                    -q2(1), -q2(2), -q2(3), q2(4)];
    qMul = (intermMatrix * q1.').';
end