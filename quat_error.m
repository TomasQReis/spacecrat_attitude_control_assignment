% Verified. 
% Quaternion matrix error function. 
% Follows the notation where q_4 is the real part of the quaternion. 
function quatError = quat_error(quatCurrent, quatTarget)
    intermMatrix = [quatTarget(4), quatTarget(3), -quatTarget(2), -quatTarget(1);
                    -quatTarget(3), quatTarget(4), quatTarget(1), -quatTarget(2);
                    quatTarget(2), -quatTarget(1), quatTarget(4), -quatTarget(3);
                    quatTarget(1), quatTarget(2), quatTarget(3), quatTarget(4)];
    quatError = (intermMatrix * quatCurrent.').';
end