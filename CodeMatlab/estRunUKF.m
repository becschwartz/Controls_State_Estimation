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

% System
%   x(k) = qk-1(x(k-1) + v(k-1))
%   z(k) = hk(x(k)) + w(k)

% V is variance of process noise 
% W is variance of sensor noise

% UNSCENTED TRANSFORM
% STEP 1 PRIOR UPDATE
% Define 2n sigma points (sx) for i in {0,1,2,...,2n-1}
%	sXm(k-1),i = E[x] + (sqrt(n*Var[x])),i
%   sXm(k-1),n+1 = E[x] - (sqrt(n*Var[x])),i

% Compute prior stats (with additive noise)
%   sXp(k),i = q_k-1(sXm,i) for all sigma points

% Compute prior statistics, Xpm predicted mean, Xpv predicted var
%   Xpm = sum 0 to 2n-1 of (1/2n)*sXp(k),i
%   Xpv = sum 0 to 2n-1 of (1/2n)*(sXp(k),i-Xpm(k))*(sXp(k),i-Xpm(k))' + V


% STEP 2 MEASUREMENT UPDATE
% generate sigma points for the measurements
%   sZp(k),i = h_k-1(sXp(k),i) for all sigma points

% Compute expected measurement, covariance ZZ, cross covariance XZ
%   Zpm(k) = sum 0 to 2n-1 of (1/2n)*sZp(k),i
%   ZZpv(k) = sum 0 to 2n-1 of (1/2n)*(sZp(k),i-Zpm(k))*(sZp(k),i-Zpm(k))' + W
%   XZpv(k) = sum 0 to 2n-1 of (1/2n)*(sXp(k),i-Xpm(k))*(sZp(k),i-Zpm(k))'

% Compute the Kalman gain, then Xmm measured mean and Xmv measured variance
% K(k) = XZpv(k)*inv(ZZpv(k))
% Xmm(k) = Xpm(k) + K(k)*(z(k) - Zpm(k))
% Xmv(k) = Xpv(k) - K(k)*ZZpv(k)*K(k)'

% where
%   x(k) = q_k-1(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt;
%                                  x2(k-1) + vel(k-1)*sin(x3(k-1))*dt;
%                                  x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt]
%   z(k) = hk(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k));
%                           x2(k)+ 0.5*B*sin(x3(k));
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
end


% for V, consider variance of r for measurement values
% keeping var_r as variance for theta just as a placeholder
V = diag([0.1,0.1,0.1]);    % variance of sensor noise

W = diag([0.1,0.1,0.1]);    % variance of process noise
v = [0;0;0];      % initial mean of process noise
w = 0;           % initial mean of measurement noise
nn = 2;           % from step x to x+1 is 2 steps

% Initialization of variables
Xpm = zeros(size(x0,1),1);              % Predicted Mean of x
Xmm = zeros(size(x0,1),nn);              % measured Mean of x
Xpv = zeros(size(P0,1),size(P0,2),1);   % Predicted Variance of x
Xmv = zeros(size(P0,1),size(P0,2),nn);   % Measured Variance of x
K = zeros(size(P0,1),size(P0,2),nn);      % Kalman gain
z = zeros(size(x0,1),nn);
sXm = zeros(size(x0,1),nn);
sXp = zeros(size(x0,1),nn);
Zpm = zeros(size(x0,1),1);
ZZpv = zeros(size(P0,1),size(P0,2),1);
XZpv = zeros(size(P0,1),size(P0,2),1);

Xmm(:,1) = x0;              % first measurement is x0
Xmv(:,:,1) = P0(:,:,1);            % first predicted variance is P0
if isMeas == true % if there are measurements 
    z = [x_meas,y_meas,theta];               % given z(1); z(2) in MATLAB
end

n = 3; % 3 dimensions, x,y,theta


%% 1. PRIOR UPDATE
%   x(k) = q_k-1(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt + v1(k-1);
%                                  x2(k-1) + vel(k-1)*sin(x3(k-1))*dt + v2(k-1);
%                                  x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt + v3(k-1)]

% Define 2n sigma points (sx) for i in {0,1,2,...,2n-1}
%	sXm(k-1),i = E[x] + (sqrt(n*Var[x])),i
%   sXm(k-1),n+1 = E[x] - (sqrt(n*Var[x])),i
sXm(:,1) = Xmm(:,1) + diag(sqrt(n*Xmv(:,:,1))); % [x1;y1;theta1]
sXm(:,2) = Xmm(:,1) - diag(sqrt(n*Xmv(:,:,1))); % [x2;y2;theta2]

% Compute prior stats (with additive noise)
%   sXp(k),i = q_k-1(sXm,i) for all sigma points
sXp(:,1) = [sXm(1,1) + vel*cos(sXm(3,1))*dt; 
            sXm(2,1) + vel*sin(sXm(3,1))*dt;
            sXm(3,1) + (vel/B)*tan(gamma)*dt];      
sXp(:,2) = [sXm(1,2) + vel*cos(sXm(3,2))*dt; 
            sXm(2,2) + vel*sin(sXm(3,2))*dt;
            sXm(3,2) + (vel/B)*tan(gamma)*dt];

% Compute prior statistics, Xpm predicted mean, Xpv predicted var
%   Xpm = sum 0 to 2n-1 of (1/2n)*sXp(k),i
%   Xpv = sum 0 to 2n-1 of (1/2n)*(sXp(k),i-Xpm(k))*(sXp(k),i-Xpm(k))' + V
for i = 1:nn
    Xpm = Xpm + (1/(2*n))*sXp(:,i); 
end
for i = 1:nn
    Xpv =  Xpv + (1/(2*n))*(sXp(:,i)-Xpm)*(sXp(:,i)-Xpm)';
end
Xpv = Xpv + V;

% compute the prior sigma points


%% 2. MEASUREMENT UPDATE

%   z(k) = hk(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k));
%                           x2(k)+ 0.5*B*sin(x3(k));
%                           x3(k) + w3(k)];
if isMeas == true
    % STEP 2 MEASUREMENT UPDATE
    % generate sigma points for the measurements
    %   sZp(k),i = h_k-1(sXp(k),i) for all sigma points
    sZp(:,1) = [sXp(1,1) + 0.5*B*cos(sXp(3,1))*dt; 
                sXp(2,1) + 0.5*B*sin(sXp(3,1))*dt;
                sXp(3,1)];      
    sZp(:,2) = [sXp(1,2) + 0.5*B*cos(sXp(3,2))*dt; 
                sXp(2,2) + 0.5*B*sin(sXp(3,2))*dt;
                sXp(3,2)];

    % Compute expected measurement, covariance ZZ, cross covariance XZ
    for i = 1:nn
        Zpm = Zpm + (1/(2*n))*sXp(:,i); 
    end

    for i = 1:nn
        ZZpv =  ZZpv + (1/(2*n))*(sZp(:,i)-Zpm)*(sZp(:,i)-Zpm)';
    end
    ZZpv = ZZpv + W;

    for i = 1:nn
        XZpv =  XZpv + (1/(2*n))*(sXp(:,i)-Xpm)*(sZp(:,i)-Zpm)';
    end
    

    % Compute the Kalman gain, then Xmm measured mean and Xmv measured variance
    K(:,:,2) = XZpv*inv(ZZpv);
    Xmm(:,2) = Xpm + K(:,:,2)*(z' - Zpm); 
    Xmv(:,:,2) = Xpv - K(:,:,2)*ZZpv*K(:,:,2)';
    Xmv = Xmv(:,:,2); 
    x = Xmm(1,2);
    y = Xmm(2,2);
    theta = Xmm(3,2);
else
    Xmv = Xmv(:,:,1);
    x = Xpm(1);
    y = Xpm(2);
    theta = Xpm(3);
end




%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = x;
disp('x')
disp(x)
internalStateOut.y = y;
internalStateOut.theta = theta;
%disp(theta)
internalStateOut.Xmv = Xmv;
%disp(Xmv)



end

