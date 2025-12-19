% FUNC1.M

function derivx = func1(t,x)
derivx(1) =-0.3*cos(x(1))^3;
derivx = derivx';