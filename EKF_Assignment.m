clc;
clear all;
close all;
randn('seed', 7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------
% Processing Cm
%----------------------------------------------
function Cm_new = Cm_processing(Cm)
    Cm_new = Cm;

    % Manually removing the jump in sensor measurement
    for i = 1001:3001
        diff1 = Cm(1001)-Cm(1000);
        Cm_new(i, 1) = Cm(i)-diff1;
    end
    
    for i = 4001:4301
        diff2 = Cm(4701)-Cm(4700);
        Cm_new(i, 1) = Cm(i)-diff2;
    end
    
    for i = 4701:5001
        diff3 = Cm(4701)-Cm(4700);
        Cm_new(i, 1) = Cm(i)-diff3;
    end
    
    for i = 5201:5401
        diff4 = Cm(5201)-Cm(5200);
        Cm_new(i, 1) = Cm(i)-diff4;
    end
    
    for i = 7401:8201
        diff5 = Cm(7401)-Cm(7400);
        Cm_new(i, 1) = Cm(i)-diff5;
    end
    
    for i = 9101:9201
        diff6 = Cm(9101)-Cm(9100);
        Cm_new(i, 1) = Cm(i)-diff6;
    end
end

%----------------------------------------------
% Runge Kutta Integration Estimator
%----------------------------------------------
function [t,w] = rk4(fn, xin, uin, t)
    a = t(1); 
    b = t(2);
    w = xin;
    N = 4;
    h = (b-a) / N;
    t = a;

    for j=1:N
        K1 = h * fn(w, uin);
        K2 = h * fn(w+K1/2, uin);
        K3 = h * fn(w+K2/2, uin);
        K4 = h * fn(w+K3, uin);

        w = w + (K1 + 2*K2 + 2*K3 + K4) / 6;
        t = a + j*h;
    end
end

%----------------------------------------------
% Function for xdot 
%----------------------------------------------
function xdot = kf_calc_f(x, u)
    % x -> (4, 1)
    n = 4;

    B = [1 0 0;
        0 1 0;
        0 0 1;
        0 0 0];
    
    xdot = zeros(n, n)*x + B*u;
end

%----------------------------------------------
% Function for matrix Z prediction
%----------------------------------------------
function zpred = kf_calc_h(t, x)
    
    u = x(1, 1);
    v = x(2, 1);
    w = x(3, 1);
    Calpha_up = x(4, 1);

    zpred = [atan(w/u)*(1+Calpha_up);
             atan(v/sqrt(u^2+w^2));
             sqrt(u^2+w^2+v^2)];
end

%----------------------------------------------
% Function for matrix Jacobian of H 
%----------------------------------------------

function hx = kf_calc_Hx(t, x)
    u = x(1, 1);
    v = x(2, 1);
    w = x(3, 1);
    Calpha_up = x(4, 1);

    C11 = (-w/(u^2+w^2))*(1+Calpha_up);
    C12 = 0;
    C13 = (u/(u^2+w^2))*(1+Calpha_up);
    C14 = atan(w/u);
    C21 = (-u*v/(u^2+w^2)^(3/2))*(1/(1+(v/sqrt(u^2+w^2))^2));
    C22 = (1/sqrt(u^2+w^2))*(1/(1+(v/sqrt(u^2+w^2))^2));
    C23 = (-w*v/(u^2+w^2)^(3/2))*(1/(1+(v/sqrt(u^2+w^2))^2));
    C31 = u/sqrt(u^2+v^2+w^2);
    C32 = v/sqrt(u^2+v^2+w^2);
    C33 = w/sqrt(u^2+v^2+w^2);

    hx = [C11 C12 C13 C14;
        C21 C22 C23 0;
        C31 C32 C33 0];
end

%----------------------------------------------
% True values of states using Z predictions and 
% reconstructed alpha
%----------------------------------------------

function x_true = true_values(Z, Z_pred)

    alpha_true = Z_pred;
    beta = Z(1);
    V = Z(2);

    denominator = sqrt(1 + (tan(beta) * sec(alpha_true))^2 + (tan(alpha_true))^2);
    
    u_true = V / denominator;
    w_true = u_true * tan(alpha_true);
    v_true = u_true * tan(beta) * sec(alpha_true);

    x_true = [u_true;
        v_true;
        w_true;
        0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
load('F16traindata_CMabV_2025.mat', 'Cm', 'Z_k', 'U_k');
load('F16validationdata_2025', 'Cm_val', 'alpha_val', 'beta_val');

% Measurements
alpha = Z_k(:, 1);
beta = Z_k(:, 2);
V = Z_k(:, 3);

% Inputs
udot = U_k(:, 1);
vdot = U_k(:, 2);
wdot = U_k(:, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial state estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ASSUMPTION: Initial measurements are close to the optimal
u0 = sqrt((V(1, 1)^2)/(1+(1+tan(alpha(1, 1))^2)*(tan(beta(1, 1)))+ tan(alpha(1, 1))^2));
w0 = u0*tan(alpha(1, 1));
v0 = sqrt(u0^2+w0^2)*tan(beta(1, 1));

Calpha_up0 = 0;

x_0 = [u0;
       v0;
       w0;
       Calpha_up0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n           = 4;        % number of states (u, v, w, Calphaup)
nm          = 3;        % number of measurements (alpha, beta, V)
m           = 3;        % number of inputs (udot, vdot, wdot)
dt          = 0.01;     % time step (s)
N           =10001;     % sample size

printfigs   = 0;
figpath     = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values for statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B           = [1 0 0;
               0 1 0
               0 0 1;
               0 0 0];          % input matrix

G           = [1 0 0 0;
               0 1 0 0;
               0 0 1 0;
               0 0 0 1];        % noise input matrix

% System noise statistics:
E_w         = 0;                           % bias of system noise
std_wu      = 1e-3;
std_wv      = 1e-3;
std_ww      = 1e-3;
std_wc      = 0;
std_w       = [std_wu^2, std_wv^2, std_ww^2, std_wc^2]; % standard deviation of system noise
Q           = diag(std_w); 
w_k         = Q * randn(n, N);              % system noise

% Measurement noise statistics:
E_v         = 0;                            % bias of measurement noise
std_va      = 1.5e-3;   
std_vb      = 1.5e-3;  
std_vv      = 1;                            % standard deviation of measurement noise
std_v       = [std_va^2, std_vb^2, std_vv^2];
R           = diag(std_v);                  % variance of measurement noise
v_k         = R * randn(3, N);              % measurement noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial covarience matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given standard deviations of inputs
sigma_V = 5;           % Uncertainty in V(1,1)
sigma_alpha = 0.1;     % Uncertainty in alpha(1,1) [rad]
sigma_beta = 0.1;      % Uncertainty in beta(1,1) [rad]
sigma_Calphaup = 0.5;  % Initial uncertainty for Calpha_up0 

% Compute partial derivatives (symbolic)
syms V11 alpha11 beta11 Calpha_up0 

% Symbolically defining the states
u0 = sqrt(V11^2 / (1 + (1 + tan(alpha11*(1+Calpha_up0))^2) * tan(beta11)^2 + tan(alpha11*(1+Calpha_up0))^2));
w0 = u0 * tan(alpha11*(1+Calpha_up0));
v0 = sqrt(u0^2 + w0^2) * tan(beta11);

% Compuing Jacobian matrix J = [du/dV, du/dalpha, du/dbeta, du/dCalphaup; ...]
J = jacobian([u0; v0; w0; Calpha_up0], [V11, alpha11, beta11, Calpha_up0]); 

% Substituting nominal values (V11 = V(1,1), etc.)
J_numeric = subs(J, {V11, alpha11, beta11, Calpha_up0}, {V(1,1), alpha(1,1), beta(1,1), 0});
J_numeric = double(J_numeric);

% Input covariance matrix
Sigma_inputs = diag([sigma_V^2, sigma_alpha^2, sigma_beta^2, sigma_Calphaup^2]);

% Output covariance matrix P0 getting full covariance matrix! :)
P_0 = J_numeric * Sigma_inputs * J_numeric';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_k             = 0; 
t_k1            = dt;

XX_k1_k1        = zeros(n, N);
PP_k1_k1        = zeros(n, N);
STD_x_cor       = zeros(n, N);
STD_z           = zeros(n, N);
ZZ_pred         = zeros(nm, N);
X_true          = zeros(m, N);

% Initial estimates for states and covariance 
x_k1_k1         = x_0;      
P_k1_k1         = P_0;      

Z_k = Z_k.'; % Inverting so it is [3, 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run the Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic; % start timer

% Run the filter through all N samples
for k = 1:N

    % 1. One step ahead prediction x(k+1|k) (prediction) [4 x 1]
    [t, x_k1_k]     = rk4(@kf_calc_f, x_k1_k1, U_k(k, :).', [t_k, t_k1]); 
    
    % 2. Calc Jacobians Phi(k+1,k) and Gamma(k+1, k) 
    Fx              = zeros(4, 4);             % perturbation of f(x,u,t) [4 x 4]
    Hx              = kf_calc_Hx(0, x_k1_k);   % perturbation of h(x,u,t) [4 x 4]

    % 3. Discretised state transition & input matrices Phi & Gamma [4 x 4]
    [Phi, Psi]    = c2d(Fx, B, dt);   
    [Phi, Gamma]    = c2d(Fx, G, dt);  

    % 4. Calculate the convariance of the state prediction error P(k+1|k) [4 x 4]
    P_k1_k          = Phi*P_k1_k1*Phi' + Gamma*Q*Gamma';  


    % Observation and observation error predictions
    z_k1_k          = kf_calc_h(0, x_k1_k);       % prediction of observation [3 x 1]
    P_zz            = (Hx*P_k1_k * Hx' + R);      % covariance matrix of observation error
    std_z           = sqrt(diag(P_zz));           % standard deviation of observation error (for validation)        


    % 5. Kalman gain calculation K(k+1) (gain) [4 x 3]
    K               = P_k1_k * Hx'*inv(P_zz);

    disp(size(K))


    % 6. Measurement update to calculate optimal state x(k+1|k+1) [4 x 1]
    x_k1_k1         = x_k1_k + K * (Z_k(:,k) - z_k1_k); 
    
   
    % 7. Calculate the covariance of the state estimation of the error
    % P(k+1|k+1) [4 x 4]
    P_k1_k1         = (eye(n) - K*Hx) * P_k1_k; 
    std_x_cor       = sqrt(diag(P_k1_k1));         % standard deviation of state estimation error (for validation)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % True states 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    alpha_m = z_k1_k(1);
    Calpha_up = x_k1_k1(4);
    alpha_true =  alpha_m/(1 + Calpha_up);

    x_true = true_values(Z_k(2:3, k), alpha_true);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updates and storing results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Next step
    t_k             = t_k1; 
    t_k1            = t_k1 + dt;

    p = 0;
    q = n;
    
    % store results
    XX_k1_k1(1:n, k)           = x_k1_k1;           % Estimated states [4 x 1]
    PP_k1_k1(1:n, k+p:k+q-1)   = P_k1_k1;           % Estimated measurements [4 x 4]
    STD_x_cor(1:n, k)          = std_x_cor;         % State std of covariance
    STD_z(1:n-1, k)            = std_z;             % Measurement std of covariance
    ZZ_pred(1, k)              = alpha_true;        % Corrected alpha [1 x 1]
    ZZ_pred(2:n-1, k)          = z_k1_k(2:3);       % Predicted measurements [2 x 1]
    X_true(1:n, k)             = x_true;            % True states u, v, w using corrected alpha [3 x 1]

    % Index update 
    p = p + n;
    q = q + n;  
end

time = toc; % end timer

disp("C_alpha estimate:");
disp(XX_k1_k1(4, 10001));

t = linspace(1, N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - Kalman Filter Results
% UNCOMMENT TO VIEW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %--------------------------------------------
% % STATES
% %--------------------------------------------
% figure; 
% subplot(2, 2, 1)
% plot(t, X_true(1, :) ,'Color', [0.902, 0.624, 0]); 
% hold on;
% plot(t, XX_k1_k1(1, :),'Color', [0, 0, 0]);
% legend('True Measurement', 'Prediction', 'Location', 'southwest');
% xlabel('Steps');
% ylabel('u [m/s]');
% title('True Measurement vs Prediction of u');
% grid on;
% 
% subplot(2, 2, 2)
% plot(t, X_true(2, :),'Color', [0.902, 0.624, 0]);
% hold on;
% plot(t, XX_k1_k1(2, :),'Color', [0, 0, 0]);
% legend('True Measurement', 'Prediction', 'Location', 'southwest');
% xlabel('Steps');
% ylabel('v [m/s]');
% title('True Measurement vs Prediction of v');
% grid on;
% 
% subplot(2, 2, 3)
% plot(t, X_true(3, :),'Color', [0.902, 0.624, 0]);
% hold on;
% plot(t, XX_k1_k1(3, :),'Color', [0, 0, 0]);
% legend('True Measurement', 'Prediction', 'Location', 'northwest');
% xlabel('Steps');
% ylabel('w [m/s]');
% title('True Measurement vs Prediction of w');
% grid on;
% 
% subplot(2, 2, 4)
% plot(t, XX_k1_k1(4, :),'Color', [0, 0, 0]);
% legend('Prediction', 'Location', 'northeast');
% xlabel('Steps');
% ylabel('C_{\alpha_{upwash}}');
% title('C_{\alpha_{upwash}} Prediction');
% grid on; 
% 
% %--------------------------------------------
% % MEASUREMENTS
% %--------------------------------------------
% 
% figure;
% subplot(2, 2, 1)
% plot(t, Z_k(1, :),'Color', [0.902, 0.624, 0]);
% hold on;
% plot(t, ZZ_pred(1, :),'Color', [0, 0, 0]);
% legend('Measurement', 'Corrected prediction', 'Location', 'northwest');
% xlabel('Steps');
% xlabel('Steps');
% ylabel('\alpha [rads]');
% title('True Measurement vs Prediction of \alpha');
% grid on;
% 
% subplot(2, 2, 2)
% plot(t, Z_k(2, :),'Color', [0.902, 0.624, 0]);
% hold on;
% plot(t, ZZ_pred(2, :),'Color', [0, 0, 0]);
% legend('Measurement', 'Prediction', 'Location', 'northwest');
% xlabel('Steps');
% ylabel('\beta [rads]');
% title('True Measurement vs Prediction of \beta');
% grid on;
% 
% subplot(2, 2, 3)
% plot(t, Z_k(3, :),'Color', [0.902, 0.624, 0]);
% hold on;
% plot(t, ZZ_pred(3, :),'Color', [0, 0, 0]);
% legend('Measurement', 'Prediction', 'Location', 'northwest');
% xlabel('Steps');
% ylabel('Velocity V [m/s]');
% title('True Measurement vs Prediction of the Velocity V');
% grid on;
% 
% figure;
% plot(Z_k(1, :), Z_k(2, :), 'Color', [0.902, 0.624, 0]);
% hold on;
% plot(ZZ_pred(1, :), ZZ_pred(2, :), 'Color', [0, 0, 0]);
% xlabel('Alpha [rads]');
% ylabel('Beta [rads]');
% legend('Measurement', 'Prediction');
% title('F16 CM(\alpha_m, \beta_m) Raw vs EKF measurements');
% 
% %--------------------------------------------
% % 3D PLOT
% %--------------------------------------------
% 
% Cm_new = Cm_processing(Cm);
% 
% % creating triangulation (only used for plotting here)
% TRIeval = delaunayn(ZZ_pred([1 2], :).');
% 
% %   viewing angles
% az = 140;
% el = 36;
% 
% figure();
% trisurf(TRIeval, ZZ_pred(1,:).', ZZ_pred(2, :).', Cm_new, 'EdgeColor', 'none'); 
% colormap('cool');
% grid on;
% hold on;
% % plot data points
% plot3(ZZ_pred(1,:).', ZZ_pred(2, :).', Cm_new, '.k'); % note that alpha_m = alpha, beta_m = beta, y = Cm
% view(az, el);
% ylabel('beta [rad]');
% xlabel('alpha [rad]');
% zlabel('C_m [-]');
% title('F16 CM(\alpha_m, \beta_m) Interpolation from EKF Results');
