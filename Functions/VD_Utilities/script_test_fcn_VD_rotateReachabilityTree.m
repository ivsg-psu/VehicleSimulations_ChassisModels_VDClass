% script_test_fcn_VD_rotateReachabilityTree.m
% tests fcn_VD_rotateReachabilityTree.m

% REVISION HISTORY:
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_rotateReachabilityTree
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
initialStates = [0 0 45*pi/180]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 1;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
vehicleParameters.L = 2.5; % wheelbase in meters

modelIDToUse = 1; % Kinematic bicycle model

% Call the function to get the reachability tree
[stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, (-1));


stateTrajectoriesCellArrayIndices = fcn_DebugTools_breakArrayByNans(stateTrajectories,-1);
NstateTrajectories = length(stateTrajectoriesCellArrayIndices);
NtimePoints = length(stateTrajectoriesCellArrayIndices{1}(:,1));

% Initialize the transform matrix
transformationMatricesReachabilityTree = cell(1,NstateTrajectories);

for ith_trajectory = 1:NstateTrajectories
	thisIndices = stateTrajectoriesCellArrayIndices{ith_trajectory};
	thisTrajectory = stateTrajectories(thisIndices,:);
	
	% Call the function to get transforms for each trajectory
	thisTransformationMatrices = ...
		fcn_VD_convertTrajectoryToRelativeTransform(thisTrajectory, (figNum));
	
	% Convert from cell array into one matrix
	thisStackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(thisTransformationMatrices, (-1));

	% Remove the NaN rows
	thisStackedMatrixNoNans = thisStackedMatrix(~isnan(thisStackedMatrix(:,1)),:);

	% Reshape this matrix into 4x4xM where M is the number of time steps
	thisStackedMatrixReshapedTransposed = reshape(thisStackedMatrixNoNans',4,4,[]);

	% Transpose the matrix back
	% EXAMPLE:
	% % M is m x n x p
	% Mt = permute(M, [2 1 3]);   % Mt is n x m x p, each slice Mt(:,:,k) = M(:,:,k).'
	thisStackedMatrixReshaped = permute(thisStackedMatrixReshapedTransposed, [2 1 3]);

	% Check the number of time points is as expected
	assert(size(thisStackedMatrixReshaped,3)==NtimePoints);

	% Save the results
	transformationMatricesReachabilityTree{1,ith_trajectory} = thisStackedMatrixReshaped;
end

%%%
% Rotate the matrix

% Convert the transformationMatrices by a rotation
rotations = [0 0 45]*pi/180; % radians
translations = [ 0 0 0]; % 1 2 3]; % radians
rotationMatrix = fcn_VD_createTransformMatrix( rotations, translations, (-1));

% Initialize the output
transformationMatricesReachabilityTreeRotated = cell(1,NstateTrajectories);

for ith_trajectory = 1:NstateTrajectories
	% Rotate the matrix
	thisStackedMatrixReshapedRotated = pagemtimes(rotationMatrix, transformationMatricesReachabilityTree{1,ith_trajectory});  
	transformationMatricesReachabilityTreeRotated{1,ith_trajectory} = thisStackedMatrixReshapedRotated;
end

%% Plot the reachability tree
figure(figNum);
clf;

% Grab the angles
anglesInRadians = nan(NstateTrajectories,1);
for ith_trajectory = 1:NstateTrajectories

	thisIndices = stateTrajectoriesCellArrayIndices{1,ith_trajectory};
	thisTrajectory = stateTrajectories(thisIndices,:);
	anglesInRadians(ith_trajectory,1) = steeringAnglesUsed(thisIndices(1,1),1);

end



% Plot the pre-transform version
h_plot = fcn_VD_plotTrajectory(stateTrajectories(:,1:2),(figNum));
set(h_plot,'DisplayName','Calculated Trajectories','LineWidth',5);


% Grab the initial states
initialStates = thisTrajectory(1,:); % global [X Y yaw] in [m m rad]
initialStatesHomogenous = [initialStates 1]';

% Plot the non-rotated version
for ith_trajectory = 1:NstateTrajectories

	% Fill in copies of initial position
	thisStackedMatrixReshaped = transformationMatricesReachabilityTree{1,ith_trajectory};
	thisTrajectoryHomogenous = squeeze(pagemtimes( thisStackedMatrixReshaped, initialStatesHomogenous))';

	h_plot = fcn_VD_plotTrajectory(thisTrajectoryHomogenous(:,1:2),(figNum));
	set(h_plot,'DisplayName',sprintf('%.1f deg',anglesInRadians(ith_trajectory,1)*180/pi), 'LineWidth',3);
end

% Plot the rotated version
for ith_trajectory = 1:NstateTrajectories

	% Fill in copies of initial position
	thisStackedMatrixReshaped = transformationMatricesReachabilityTreeRotated{1,ith_trajectory};
	thisTrajectoryHomogenous = squeeze(pagemtimes( thisStackedMatrixReshaped, initialStatesHomogenous))';

	h_plot = fcn_VD_plotTrajectory(thisTrajectoryHomogenous(:,1:2),(figNum));
	set(h_plot,'DisplayName',sprintf('%.1f deg',anglesInRadians(ith_trajectory,1)*180/pi), 'LineWidth',3);
end

%%

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(transformationMatrices));

% Check variable sizes
assert(size(transformationMatrices,1)==length(stateTrajectory(:,1))); 
assert(size(transformationMatrices,2)==1); 

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%%%%
%  Compare speeds 

Niterations = 10;

% Do calculation via RK4 solver
tic;
for ith_test = 1:Niterations
    [stateTrajectory, ~, ~] = ...
        fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
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

%%%%
%  now force to use the GPU
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

%% DEMO case: basic call at 45 degrees post call
figNum = 10002;
titleString = sprintf('DEMO case: basic call at 45 degrees post call');
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
fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_rotateReachabilityTree(stateTrajectory, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(transformationMatrices));

% Check variable sizes
assert(size(transformationMatrices,1)==length(stateTrajectory(:,1))); 
assert(size(transformationMatrices,2)==1); 

% Convert the transformationMatrices by a rotation
rotations = [0 0 45]*pi/180; % radians
translations = [ 0 0 0]; % 1 2 3]; % radians
rotationMatrix = fcn_VD_createTransformMatrix( rotations, translations, (-1));

transformationMatricesRotated = cellfun(@(X) rotationMatrix*X, transformationMatrices, 'UniformOutput', false);

stackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(transformationMatricesRotated, (-1));
stackedMatrixNoNans = stackedMatrix(~isnan(stackedMatrix(:,1)),:);




% Set the simulation time/state arguments again with 45 degree inputs
initialStates = [0 0 45*pi/180]; % [X Y phi] in [m],[m],[rad]
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
fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));


% Fill in copies of initial position
initialStatesHomogenous = [initialStates 1];

% Calculate predicted positions
predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
predictedPositions = predictedPositions_homogenousForm(:,1:3);

assert(isequal(round(stateTrajectory(:,1:2),4),round(predictedPositions(:,1:2),4)))

h_plot = fcn_VD_plotTrajectory(predictedPositions(:,1:2),(figNum));
set(h_plot,'DisplayName','XY Trajectory (rotated via transform)')

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

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
fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_rotateReachabilityTree(stateTrajectory, ([]));

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
fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

% Call the function
transformationMatrices = ...
    fcn_VD_rotateReachabilityTree(stateTrajectory, (-1));

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
fcn_VD_rotateReachabilityTree(initialStates, deltaT, ...
timeInterval, inputsVsTime, parameters, (-1));

Niterations = 10;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    transformationMatrices = ...
        fcn_VD_rotateReachabilityTree(stateTrajectory, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    transformationMatrices = ...
        fcn_VD_rotateReachabilityTree(stateTrajectory, ([]));
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
