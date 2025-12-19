% FUNC3.M

function [Phi,Gamma] = func3(x,t,n,Ts)
F=0.9*cos(x)^2*sin(x);
G = [1];
[Phi,Gamma]=c2d(F,G,Ts);