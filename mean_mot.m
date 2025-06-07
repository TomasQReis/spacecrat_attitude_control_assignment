% Verified. 
% Mean motion of orbit with given semiMajor around body mu. 
function meanMot = mean_mot(mu, semiMajor)
    meanMot = sqrt(mu / semiMajor^3);
end