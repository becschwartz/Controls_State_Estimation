function [internalState, studentNames, estimatorType] = estInitialize
% Fill in whatever initialization you'd like here. This function
% generates the internal state of the estimator at time 0. You may do
% whatever you like here, but you must return something that is in the
% format as may be used by your run() function as the first returned variable.
%
% The second returned variable must be a list of student names.
%
% The third return variable must be a string with the estimator type

% we make the internal state a structure, with the first three elements the
% positions x, y; the angle theta; and our favorite colour.

% note that there is *absolutely no prescribed format* for this internal state.
% You can put in it whatever you like. Probably, you'll want to keep the position
% and angle, and probably you'll remove the color.
internalState.x = 0;
internalState.y = 0;
internalState.theta = pi/4; % rad

studentNames = ['Kirsten Biermayer, Rebecca Schwartz']

% replace this with the estimator type. Use one of the following options:
%  'EKF' for Extended Kalman Filter
%  'UKF' for Unscented Kalman Filter
%  'PF' for Particle Filter
%  'OTHER: XXX' if you're using something else, in which case please
%                 replace "XXX" with a (very short) description
estimatorType = 'EKF'

end


