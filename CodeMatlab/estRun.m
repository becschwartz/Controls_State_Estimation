function [x,y,theta,internalStateOut] = estRun(time, dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
% In this function you implement your estimator. The function arguments
% are:
%  time: current time in [s]
%  dt: current time step [s]
%  internalStateIn: the estimator internal state, definition up to you.
%  steeringAngle: the steering angle of the bike, gamma, [rad]
%  pedalSpeed: the rotational speed of the pedal, omega, [rad/s]
%  measurement: the position measurement valid at the current time step
%
% Note: the measurement is a 2D vector, of x-y position measurement.
%  The measurement sensor may fail to return data, in which case the
%  measurement is given as NaN (not a number).
%
% The function has four outputs:
%  est_x: your current best estimate for the bicycle's x-position
%  est_y: your current best estimate for the bicycle's y-position
%  est_theta: your current best estimate for the bicycle's rotation theta
%  internalState: the estimator's internal state, in a format that can be understood by the next call to this function
% this needs to correspond to your init function:

% INIT
x = internalStateIn.x;
y = internalStateIn.y;
theta = internalStateIn.theta;

if time ~= 0
    Xmv = internalStateIn.Xmv;
end

% GIVEN
gamma = steeringAngle;

% VALUES WITH KNOWN UNCERTAINTY
% Tire Radius
r = 0.425; % tire radius (m); uncertainty is +- 5%
r_low = r - 0.05*r;
r_high = r + 0.05*r;
var_r = ((r_high - r)/3)^2; % GAUSSIAN NOISE: where 3 is the desired z value on z scale
vel = 5*pedalSpeed*r; % expected value, non directional linear speed of bike

% Wheel Base
B = 0.8; % wheel base length (m); uncertainty is +- 10%
B_low = B - 0.1*B;
B_high = B + 0.1*B;
var_B = ((B_high - B)/3)^2; % GAUSSIAN NOISE: where 3 is the desired z value on z scale

% VALUES WITH UNCERTAINTY, EXACT UNCERTAINTY UNKNOWN
% ...

% MEASUREMENT GIVEN (z(k))
isMeas = true; % measurement assumed to be given
if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    x_meas = measurement(1);
    y_meas = measurement(2);
else
    isMeas = false; % no measurement given
    disp('nomeas')
end


%%
% Design UKF for this 3 state nonlin sys and compute its posteriori
% estimate, Xmm = [x;y;theta] = [x1,x2,x3] , given measurements 

% BECCA NOTE
% UKF, Unscented Kalman Filter is the Kalman filter for nonlinear systems. 

% Use the uncented transform to predict E[x(1)] and Var[x(1)]

% UNSCENTED TRANSFORM
% STEP 1 PRIOR UPDATE
% Define 2n sigma points (sx) for i in {0,1,2,...,2n-1}
%	sXm,i = E[x] + (sqrt(n*Var[x])),i
%   sXm,n+1 = E[x] - (sqrt(n*Var[x])),i

% Compute prior stats (with additive noise)
%   sXp,i = q_k-1(sXm,i) for all sigma points

% STEP 2 MEASUREMENT UPDATE
% 

% STEP 3
% Approximate the mean and variance 
%   E[y] ~= sum from 0 to 2n-1 of (1/2n)*sy,i
%   Var[y] ~= sum from 0 to 2n-1 of (1/2n)*(sy,i-E[y])*(sy,i-E[y])'

%   x(k) = q(x(k-1),u(k-1),v(k-1))
%   z(k) = h(x(k),w(k))


% where
%   x(k) = q_k-1(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt + v1(k-1);
%                                  x2(k-1) + vel(k-1)*sin(x3(k-1))*dt + v2(k-1);
%                                  x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt + v3(k-1)]
%   z(k) = hk(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k)) + w1(k);
%                           x2(k)+ 0.5*B*sin(x3(k)) + w1(k);
%                           x3(k) + w3(k)];


%% INIT

% x0 = [x;y;theta];
x1 = x;
x2 = y;
x3 = theta;
x0 = [x1;x2;x3];       % initial mean

% CHANGE 
% NOTE: variances are def higher i think
% for P0, consider variance for B
% keeping var_B as variance for theta just as a placeholder

if time == 0
    P0 = diag([var_r + var_B, var_r + var_B, var_B]);   % initial variance
else
    P0 = Xmv;
    % disp(Xmv)
end


% for V, consider variance of r for measurement values
% keeping var_r as variance for theta just as a placeholder
V = diag([0.01,0.01,0.01]);    % variance of sensor noise

W = diag([0.01,0.01,0.01]);    % variance of process noise
v = [0;0;0];      % initial mean of process noise
w = 0;           % initial mean of measurement noise
n = 2;           % from step x to x+1 is 2 steps

% Initialization of variables
Xpm = zeros(size(x0,1),n);              % Predicted Mean of x
Xmm = zeros(size(x0,1),n);              % measured Mean of x
Xpv = zeros(size(P0,1),size(P0,2),n);   % Predicted Variance of x
Xmv = zeros(size(P0,1),size(P0,2),n);   % Measured Variance of x
K = zeros(size(P0,1),size(P0,2),n);      % Kalman gain
z = zeros(size(x0,1),n);
sXm = zeros(size(x0,1),n);

Xmm(:,1) = x0;              % first measurement is x0
Xmv(:,:,1) = P0(:,:,1);            % first predicted variance is P0
if isMeas == true % if there are measurements 
    z = [x_meas,y_meas,theta];               % given z(1); z(2) in MATLAB
end

%% 1. PRIOR UPDATE
% STEP 1
% Define 2n sigma points (sx) for i in {0,1,2,...,2n-1}
%	sx,i = E[x] + (sqrt(n*Var[x])),i
%   sx,i+1 = E[x] - (sqrt(n*Var[x])),i

% S1, Define 2n sigma points
% only have 6 sigma points since x has x,y,theta (n = 3)

% generate sigma points
sXm(:,1) = Xmm(:,1) + (sqrt(n*Xmv(:,:,1)));
sXm(:,2) = Xmm(:,1) - (sqrt(n*Xmv(:,:,1)));

% compute the prior sigma points
sy(1) = -sx(1) +2*abs(sx(1));
sy(2) = -sx(2) +2*abs(sx(2));

%% 2. MEASUREMENT UPDATE

% S2, Define transform of signma points




%% 1d 
% SOLVE IT YO: xmhat(1) given z(1) = 0.5 


%   q(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt + v1(k-1);
%                       x2(k-1) + vel(k-1)*sin(x3(k-1))*dt + v2(k-1);
%                       x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt + v3(k-1)] 

% Xpm(:,2) = q(Xmm(k-1),v(k-1) = 0)
Xpm(:,2) = [Xmm(1,1) + vel*cos(Xmm(3,1))*dt;
            Xmm(2,1) + vel*sin(Xmm(3,1))*dt;
            Xmm(3,1) + (vel/B)*tan(gamma)*dt]; % Solve for the predicted Mean
        
if isMeas == true % if there are measurement values, update... 
%   A(k-1) = pd_q_k-1_(Xmm(k-1),u(k-1),v(k-1)) / pd(x)      % Jacobian A
     A = [1, 0, 0; 
          0, 1, 0;
         -vel*dt*sin(Xmm(3)), vel*dt*cos(Xmm(3)), 1];
     
%   L(k-1) = pd_q_k-1_(Xmm(k-1),u(k-1),v(k-1)) / pd(v)      % Jacobian L
    L = eye(3);    % JACOBIAN L
    
    Xpv(:,:,2) = A*Xmv(:,:,1)*A + L*V*L'; % Solve for the predicted Variance

    %   H(k) = pd_h_k(Xpm(k),w(k)) / pd(x)
    H = [1, 0, 0;
         0, 1, 0;
        -0.5*B*sin(Xmm(3)), 0.5*B*cos(Xmm(3)), 1]; % JACOBIAN H
    
%   M(k) = pd_h_k(Xpm(k),w(k)) / pd(w)   
    M = eye(3); % JACOBIAN M
    
    K(:,:,2) = (Xpv(:,:,2)*H')*inv(H*Xpv(:,:,2)*H' + M*W*M');    % Kalman Gain
%   h(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k)) + w1(k);
%                   x2(k)+ 0.5*B*sin(x3(k)) + w1(k);
%                   x3(k) + w3(k)];  

    % Xmm(:,k) = Xpm(:,k) + K(k)*(z(k) - hk(Xpm(k),w(k) = 0))    
    Xmm(:,2) = Xpm(:,2) + K(:,:,2)*(z' - [Xpm(1,2)+ 0.5*B*cos(Xpm(3,2)); Xpm(2,2)+ 0.5*B*sin(Xpm(3,2)); Xpm(3,2)]); % Solve for the measured mean
    % how do we use this?
    Xmv(:,:,2) = (eye(size(K,1)) - K(:,:,2)*H)*Xpv(:,:,2);       % solve for measured variance of y
    Xmv = Xmv(:,:,2); 
    x = Xmm(1,2);
    y = Xmm(2,2);
    theta = Xmm(3,2);
else
    Xmv = Xmv(:,:,1);
    x = Xpm(1,2);
    y = Xpm(2,2);
    theta = Xpm(3,2);
end



%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = x;
internalStateOut.y = y;
internalStateOut.theta = theta;
disp(theta)
internalStateOut.Xmv = Xmv;
disp(Xmv)



end

