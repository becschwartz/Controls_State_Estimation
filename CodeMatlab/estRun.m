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

% GIVEN
gamma = steeringAngle;

% VALUES WITH KNOWN UNCERTAINTY
% Tire Radius
r = 0.425; % tire radius (m); uncertainty is +- 5%
r_low = r - 0.05*r;
r_high = r + 0.05*r;
var_r = 0.5*((r_low-r)^2 + (r_high-r)^2); % standard variance 
% NORMAL/ UNIFORM VARIANCE? 
vel = 5*pedalSpeed*r; % expected value, non directional linear speed of bike

% Wheel Base
B = 0.8; % wheel base length (m); uncertainty is +- 10%
B_low = B - 0.1*B;
B_high = B + 0.1*B;
var_B = 0.5*((B_low-r)^2 + (B_high-r)^2); % standard variance
% NORMAL/ UNIFORM VARIANCE? 


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
end


%%
% Design an EKF for this 2 state nonlin sys and compute its posteriori
% estimate, Xmm = [x;y;theta] = [x1,x2,x3] , given measurements 

% BECCA NOTE
% EKF, Extended Kalman Filter is the Kalman filter for nonlinear systems. 

%   x(k) = q(x(k-1),u(k-1),v(k-1))
%   z(k) = h(x(k),w(k))

% note that there is uncertainty in all these 
% where
%   q(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt + v1(k-1);
%                       x2(k-1) + vel(k-1)*sin(x3(k-1))*dt + v2(k-1);
%                       x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt + v3(k-1)]
%   h(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k)) + w1(k);
%                   x2(k)+ 0.5*B*sin(x3(k)) + w1(k);
%                   x3(k) + w3(k)];


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
P0 = diag([var_r + var_B, var_r + var_B, var_B]);        % initial variance
% for V, consider variance of r for measurement values
% keeping var_r as variance for theta just as a placeholder
V = diag([0.1,0.1,0.1]);         % variance of sensor noise

W = diag([0,0,0]);    % variance of process noise
v = [0;0;0];      % initial mean of process noise
w = 0;           % initial mean of measurement noise
n = 2; 

% Initialization of variables
Xpm = zeros(size(x0,1),n);              % Predicted Mean of x
Xmm = zeros(size(x0,1),n);              % measured Mean of x
Xpv = zeros(size(P0,1),size(P0,2),n);   % Predicted Variance of x
Xmv = zeros(size(P0,1),size(P0,2),n);   % Measured Variance of x
K = zeros(size(P0,1),size(P0,2),n);      % Kalman gain
z = zeros(size(x0,1),n);

Xmm(:,1) = x0;              % first measurement is x0
Xmv(:,:,1) = P0;            % first predicted variance is P0
if isMeas == true % if there are measurements 
    z = [x_meas,y_meas,theta];               % given z(1); z(2) in MATLAB
end

%% 1b
% Prior with Jacobian matrices A,L

% THE PRIOR UPDATE
%   Xpm(k) = q(k-1)(Xmm(k-1),u(k-1),v(k-1))                     % predicted mean
%   Xpv(k) = A(k-1)*Xmm(k-1)*A'(k-1) + L(k-1)*V(k-1)*L'(k-1)  % predicted variance
% where pd is partial deriv...
%   A(k-1) = pd_q_k-1_(Xmm(k-1),u(k-1),v(k-1)) / pd(x)      % Jacobian A
%   L(k-1) = pd_q_k-1_(Xmm(k-1),u(k-1),v(k-1)) / pd(v)      % Jacobian L

% Remember, we have
%   x(k) = q(x(k-1),v(k-1))
%   z(k) = h(x(k),w(k))
% where
%   q(x(k-1),v(k-1)) = [x1(k-1) + vel(k-1)*cos(x3(k-1))*dt + v1(k-1);
%                       x2(k-1) + vel(k-1)*sin(x3(k-1))*dt + v2(k-1);
%                       x3(k-1) + (vel(k-1)/B)*tan(gamma(k-1))*dt + v3(k-1)]

% A calculated using partial derivatives... 
% A1 = pd(q(x,v))/x1 = [1; 0; -vel(k-1)*dt*sin((x3(k-1))]
% A2 = pd(q(x,v))/x2 = [0; 1; vel(k-1)*dt*cos((x3(k-1))]
% A3 = pd(q(x,v))/x3 = [0; 0; 1]

% JACOBIAN A MATRIX
% A(k-1) = [1, 0, 0; 
%           0, 1, 0;
%           -vel(k-1)*dt*sin((x3(k-1)), vel(k-1)*dt*cos((x3(k-1)), 1];


% L calculated using partial derivatives...
% L1 = pd(q(x,v))/v1 = [1;0;0]
% L2 = pd(q(x,v))/v2 = [0;1;0]
% L3 = pd(q(x,v))/v3 = [0;0;1]

% JACOBIAN L MATRIX
% L = eye(3)

%% 1c
% Measurement with Jacobian matrices H,M

% THE MEASUREMENT UPDATE
%   K(:,k) = Xpv(:,:,k)*H(:,:,k)'*inv(H(:,:,k)*Xpv(:,:,k)*H(:,:,k)' +
%       M(:,:,k)*W*M(:,:,k)');    % Kalman Gain
%   Xmm(:,k) = Xpm(:,k) + K(:,k)*(z(:,k) - h_k(Xpm(:,k),w(k));   % mean meas
%   Xmv(:,:,k) = (eye(size(K,1)) - K(:,k)*H(:,:,k)*Xpv(:,:,k); % meav var
% where pd is partial deriv...
%   H(k) = pd_h_k(Xpm(k),w(k)) / pd(x)
%   M(k) = pd_h_k(Xpm(k),w(k)) / pd(w)

% Remember, we have
%   x(k) = q(x(k-1),v(k-1))
%   z(k) = h(x(k),w(k))
% where
%   z = h(x(k),w(k)) = [x1(k)+ 0.5*B*cos(x3(k)) + w1(k);
%                       x2(k)+ 0.5*B*sin(x3(k)) + w2(k);
%                       x3(k) + w3(k)];

% H calculated using partial derivatives... 
% H1 = pd(h(xp,w))/x1 = [1; 0; -0.5*B*sin(x3(k))]
% H2 = pd(h(xp,w))/x2 = [0; 1; 0.5*B*cos(x3(k))]
% H3 = pd(h(xp,w))/x3 = [0; 0; 1]

% JACOBIAN H MATRIX
% H(k) = [1, 0, 0;
%         0, 1, 0;
%        -0.5*B*sin(x3(k)), 0.5*B*cos(x3(k)), 1]


% PROBABLY CHANGE
% M calculated using partial derivatives...
% M1 = pd(h(xp,w))/w1 = [1;0;0]
% M2 = pd(h(xp,w))/w2 = [0;1;0]
% M3 = pd(h(xp,w))/w3 = [0;0;1]

% JACOBIAN M MATRIX
% M = eye(3)


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
    x = Xmm(1,2);
    y = Xmm(2,2);
    theta = Xmm(3,2);
else
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

end

