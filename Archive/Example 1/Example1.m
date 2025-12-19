clc; clear;

Ex_0    = 10;           % Estimated initial condition
stdx_0  = 10;           % Estimated initial standard deviation
P_0     = stdx_0^2;     % Estimated initial covariance
x_0     = 3;            % Real initial condition

% Simulation data:
dt      = 0.01;         % Time step
N       = 500;          % Number of steps
t_end   = N*dt;         % End time

% System noise:
Ew      = 0;
stdw    = 10;
Q       = stdw^2;
w_k     = stdw*randn(1,N);

% Measurement noise:
Ev      = 0;
stdv    = 10;
R       = stdv^2;
v_k     = stdv*randn(1,N);

% Real simulated state-variable and measurements data:
x       = x_0;
x_k     = zeros(1,N);   % Initialise vector that will hold the states
z_k     = zeros(1,N);   % Initialise vector that will hold the outputs

for i=1:N
   dx       = func1(i,x);
   x        = x +(dx+w_k(i))*dt;
   x_k(i)   = x;
   z_k(i)   = func2(x)+v_k(i);
end

% Save data in datafile data1.mat
%---------------------------------
save data1 x_0 x_k z_k Ex_0 stdx_0 P_0 Ew stdw Q Ev stdv R dt N v_k w_k

% Now you can run EKF.m

