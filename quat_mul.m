% Verified. 
% Quaternion matrix multiplication function. 
% Follows the notation where q_4 is the real part of the quaternion. 
% Done simply to skip the required shift of elements from the differing
% matlab to course notation. 
function qMul = quat_mul(q1, q2)
    % quat_mul.
    %   Returns the quaternion product q = q1 ∘ q2,
    %   where each input is a 4-element row [qx, qy, qz, q0], and the output
    %   is in the same format.
    q1Mat = quaternion(q1(4), q1(1), q1(2), q1(3));
    q2Mat = quaternion(q2(4), q2(1), q2(2), q2(3));

    qMulMat = compact(q1Mat * q2Mat);
    qMul = [qMulMat(2:4), qMulMat(1)].';
end