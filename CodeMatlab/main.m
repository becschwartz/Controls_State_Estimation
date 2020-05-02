%%
% This is the main function, which will initilize your estimator, and run
% it using data loaded from a text file.

%% Changes:
% April 05: 
% - print run number to help keep track. 
% - corrected final error output to be in [-pi, pi]


%%
% Provide the index of the experimental run you would like to use. Note
% that using "0" means that you will load the measurement calibration data.

experimentalRun = 0;
fprintf(['Loading the data file #' num2str(experimentalRun) ' \n']);
filename = ['data/run_' num2str(experimentalRun,'%03d') '.csv'];
experimentalData = csvread(filename);

%%
%==========================================================================
% Here, we run your estimator's initialization
%==========================================================================

fprintf('Running the initialization \n');
[internalState, studentNames, estimatorType] = estInitialize();

%%
% Here we will store the estimated position and orientation, for later
% plotting:

numDataPoints = size(experimentalData,1);
estimatedPosition_x = zeros(numDataPoints,1);
estimatedPosition_y = zeros(numDataPoints,1);
estimatedAngle = zeros(numDataPoints,1);

fprintf('Running the system \n');
dt = experimentalData(2,1) - experimentalData(1,1);
loops = 0;
for k = 1:numDataPoints
    loops = loops + 1; 
    %fprintf('This is loop # %d\n',loops) % track number of loops
    t = experimentalData(k,1);
    gamma = experimentalData(k,2);
    omega = experimentalData(k,3);
    measx = experimentalData(k,4);
    measy = experimentalData(k,5);
    
    %run the estimator:
    [x, y, theta, internalState] = estRun(t, dt, internalState, gamma, omega, [measx, measy]);

    %keep track:
    estimatedPosition_x(k) = x;
    estimatedPosition_y(k) = y;
    estimatedAngle(k) = theta;
end    

fprintf('Done running \n');
% make sure the angle is in [-pi,pi]
estimatedAngle = mod(estimatedAngle+pi,2*pi)- pi;

posErr_x = estimatedPosition_x - experimentalData(:,6);
posErr_y = estimatedPosition_y - experimentalData(:,7);
angErr   = mod(estimatedAngle - experimentalData(:,8) + pi, 2*pi) - pi;

fprintf('Final error: \n');
fprintf(['   pos x = ' num2str(posErr_x(end)) ' m \n']);
fprintf(['   pos y = ' num2str(posErr_y(end)) ' m \n']);
fprintf(['   angle = ' num2str(angErr(end)) ' rad \n']);


ax = sum(abs(posErr_x))/numDataPoints;
ay = sum(abs(posErr_y))/numDataPoints;
ath = sum(abs(angErr))/numDataPoints;
score = ax + ay + ath;
if ~isnan(score)
    fprintf('Average error: \n');

    fprintf(['   pos x = ' num2str(ax) ' m \n']);
    fprintf(['   pos y = ' num2str(ay) ' m \n']);
    fprintf(['   angle = ' num2str(ath) ' rad \n']);
    
    % our scalar score
    fprintf(['Average score: ' num2str(score) ' \n'])
end


%%
%==========================================================================
% Make some plots:
%==========================================================================
% Feel free to add additional plots, if you like.
fprintf('Generating plots \n')

%% Figure 1
figure;
hold on;
plot(experimentalData(:,4), experimentalData(:,5), 'rx');
plot(estimatedPosition_x, estimatedPosition_y, 'b-');
plot(experimentalData(:,6), experimentalData(:,7), 'k:.');
hold off

xlabel('x-position [m]');
ylabel('y-position [m]');
legend('meas','est','true');

%% Figure 2
figure;
subplot(5,1,1);
hold on;
plot(experimentalData(:,1), experimentalData(:,6), 'k:.');
plot(experimentalData(:,1), experimentalData(:,4), 'rx');
plot(experimentalData(:,1), estimatedPosition_x, 'b-');
ylabel('Position x [m]');
legend('truth','meas','est');
hold off;
subplot(5,1,2);
hold on;
plot(experimentalData(:,1), experimentalData(:,7), 'k:.');
plot(experimentalData(:,1), experimentalData(:,5), 'rx');
plot(experimentalData(:,1), estimatedPosition_y, 'b-');
ylabel('Position y [m]');
hold off;
subplot(5,1,3);
hold on;
plot(experimentalData(:,1), experimentalData(:,8), 'k:.');
plot(experimentalData(:,1), estimatedAngle, 'b-');
ylabel('Angle theta [rad]');
hold off;
subplot(5,1,4);
plot(experimentalData(:,1), experimentalData(:,2), 'g-');
ylabel('Steering angle gamma [rad]');
subplot(5,1,5);
plot(experimentalData(:,1), experimentalData(:,3), 'g-');
ylabel('Pedal speed omega [rad/s]');
xlabel('Time [s]');

%%
fprintf('Done \n');
