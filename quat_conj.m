
% Made due to course quaternion notation not matching that of matlab. 
function quatConj = quat_conj(quaternion)
    quatConj = [-quaternion(1) -quaternion(2) -quaternion(3) quaternion(4)];
end