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

% VALUES WITH KNOWN UNCERTAINTY
r = 0.425; % tire radius (m); uncertainty is +- 5%
r_low = 4 - 0.05*r;
r_high = 4 + 0.05*4;
B = 0.8; % wheel base length (m); uncertainty is +- 10%
B_low = B - 0.1*B;
B_high = B + 0.1*B;

x = internalStateIn.x;
y = internalStateIn.y;
theta = internalStateIn.theta;
% tire_radius = r + 0.05*r*randn(1); % FOR NORMAL NOISE
tire_radius = r_low + (r_high - r_low)*rand(1); % FOR UNIFORM NOISE
v_linear = 5*pedalSpeed*tire_radius; % non directional linear speed of bike

% Below, was originally x = x + pedalSpeed; y = y + pedalSpeed
x = x + pedalSpeed; % v_linear*cos(theta)*0.5; 
y = y + pedalSpeed;% v_linear*sin(theta)*0.5;

if ~isnan(measurement(1)) & ~isnan(measurement(2))
    % have a valid measurement
    x = measurement(1);
    y = measurement(2);
    theta = theta; %+ 1;
else
    rand_val_n = randn(1); % use same randn for NORMAL NOISE
    rand_val_u = rand(1); % rand val for UNIFORM NOISE
    
    x=x+0.5*(B + 0.1*B*rand_val_n)*cos(theta); % NORMAL NOISE
    y=y+0.5*(B + 0.1*B*rand_val_n)*sin(theta); % NORMAL NOISE
    x = x + 0.5*(B_low + (B_high-B_low)*rand_val_u) * cos(theta); % UNIFORM NOISE
    y = y + 0.5*(B_low + (B_high-B_low)*rand_val_u) * sin(theta); % UNIFORM NOISE
    theta = theta; %+ 1;
end

% we're unreliable about something or another
% more code here

%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = x;
internalStateOut.y = y;
internalStateOut.theta = theta;

end

