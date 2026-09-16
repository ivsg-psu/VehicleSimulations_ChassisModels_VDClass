% script_test_fcn_VD_kinematicBicycleModelRK4.m
% tests fcn_VD_convertTrajectoryToRelativeTransform.m

% REVISION HISTORY:
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_convertTrajectoryToRelativeTransform
%   % * Wrote the code originally, 
%   % * Using script_test_fcn_VD_kinematicBicycleModelRK4 as starter


% TO-DO:
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - (fill in items here)


%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: basic call
figNum = 10001;
titleString = sprintf('DEMO case: basic call');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set up inputs
steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
simulationTimes = (startTime:deltaT:endTime)';
inputsVsTime = [simulationTimes steering_amplitude_degrees*pi/180*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]

% Set up parameters
clear parameters
parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
parameters.L = 2.5; % wheelbase in meters

% Calculate the trajectory
[stateTrajectory, ~, ~] = ...
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(transformationMatrices));

% Check variable sizes
assert(size(transformationMatrices,1)==length(stateTrajectory(:,1))); 
assert(size(transformationMatrices,2)==1); 

% Check variable values
% Convert the transformationMatrices
stackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(transformationMatrices, (-1));
stackedMatrixNoNans = stackedMatrix(~isnan(stackedMatrix(:,1)),:);

% Fill in copies of initial position
initialStatesHomogenous = [initialStates 1];

% Calculate predicted positions
predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
predictedPositions = predictedPositions_homogenousForm(:,1:3);

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%%
% Compare speeds 

Niterations = 10;

% Do calculation via RK4 solver
tic;
for ith_test = 1:Niterations
    [stateTrajectory, ~, ~] = ...
        fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
        timeInterval, inputsVsTime, parameters, (-1));
end
slow_method = toc;

% Do calculation using transforms
tic;
for ith_test = 1:Niterations
    predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
    predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
    predictedPositions = predictedPositions_homogenousForm(:,1:3);
end
fast_method = toc;

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')
title(sprintf('Ratio is %.1f times faster',slow_method/fast_method))

%% now force to use the GPU
stackedMatrixNoNans_GPU= gpuArray(stackedMatrixNoNans);

% Do calculation using GPU
tic;
for ith_test = 1:Niterations
    predictedPositionsAllPoints = stackedMatrixNoNans_GPU*initialStatesHomogenous';
    predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
    predictedPositions_GPU = predictedPositions_homogenousForm(:,1:3);
    predictedPositions = gather(predictedPositions_GPU);
end
fast_method = toc;

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')
title(sprintf('Ratio is %.1f times faster',slow_method/fast_method))

%% Test cases start here. These are very simple, usually trivial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
figNum = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set up inputs
steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
simulationTimes = (startTime:deltaT:endTime)';
inputsVsTime = [simulationTimes steering_amplitude_degrees*pi/180*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]

% Set up parameters
clear parameters
parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
parameters.L = 2.5; % wheelbase in meters

% Calculate the trajectory
[stateTrajectory, ~, ~] = ...
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(transformationMatrices));

% Check variable sizes
assert(size(transformationMatrices,1)==length(stateTrajectory(:,1))); 
assert(size(transformationMatrices,2)==1); 

% Check variable values
% Convert the transformationMatrices
stackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(transformationMatrices, (-1));
stackedMatrixNoNans = stackedMatrix(~isnan(stackedMatrix(:,1)),:);

% Fill in copies of initial position
initialStatesHomogenous = [initialStates 1];

% Calculate predicted positions
predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
predictedPositions = predictedPositions_homogenousForm(:,1:3);

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set up inputs
steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
simulationTimes = (startTime:deltaT:endTime)';
inputsVsTime = [simulationTimes steering_amplitude_degrees*pi/180*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]

% Set up parameters
clear parameters
parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
parameters.L = 2.5; % wheelbase in meters

% Calculate the trajectory
[stateTrajectory, ~, ~] = ...
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(transformationMatrices));

% Check variable sizes
assert(size(transformationMatrices,1)==length(stateTrajectory(:,1))); 
assert(size(transformationMatrices,2)==1); 

% Check variable values
% Convert the transformationMatrices
stackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(transformationMatrices, (-1));
stackedMatrixNoNans = stackedMatrix(~isnan(stackedMatrix(:,1)),:);

% Fill in copies of initial position
initialStatesHomogenous = [initialStates 1];

% Calculate predicted positions
predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
predictedPositions = predictedPositions_homogenousForm(:,1:3);

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set up inputs
steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
simulationTimes = (startTime:deltaT:endTime)';
inputsVsTime = [simulationTimes steering_amplitude_degrees*pi/180*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]

% Set up parameters
clear parameters
parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
parameters.L = 2.5; % wheelbase in meters

% Calculate the trajectory
[stateTrajectory, ~, ~] = ...
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    transformationMatrices = ...
        fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    transformationMatrices = ...
        fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, ([]));
end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG 

%% Fail conditions
if 1==0
    
end


%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
