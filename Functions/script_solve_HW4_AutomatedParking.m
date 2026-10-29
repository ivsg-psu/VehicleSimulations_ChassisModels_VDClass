% script_solve_HW4_AutomatedParking.m
% shows how to solve the automated parking problem

% REVISION HISTORY:
%
% 2026_10_26 by Sean Brennan, sbrennan@psu.edu
% - In script_solve_HW4_AutomatedParking.m
%   % * Wrote the code originally,
%   % * Using script_HW2_FirstOrderEngineModelRK4 as starter

% TO-DO:
%
% 2026_10_26 by Sean Brennan, sbrennan@psu.edu
% - (fill in items here)


%% Set up the workspace
close all

%% Code starts here
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

%% DEMO case: simulating a 2013 Prius
figNum = 10999;
figure(figNum); clf;


% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 10;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up vehicle parameters
clear vehicleParameters
vehicleParameters.U = 2;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 5 mph ~= 2* m/s)
vehicleParameters.L = 2.700; % wheelbase in meters of a 2013 Prius, standard class
% How to get the a and b values? (used for plotting)
% Use the weight split!
% FOR PRIUS
% * Total Curb Weight: ~3,042 lbs (1,380 kg)
% * Estimated Front Weight: ~1,855 lbs (61%)
% * Estimated Rear Weight: ~1,187 lbs (39%)
vehicleParameters.a =  vehicleParameters.L * 1187/3042;
vehicleParameters.b =  vehicleParameters.L * 1855/3042;

modelIDToUse = 1; % Kinematic bicycle model

% Call the function
[stateTrajectories, times, steeringAnglesUsed] = ...
	fcn_VD_forwardReachabilityTreeRK4(...
	initialStates, ...
	deltaT, timeInterval, steeringInterval, vehicleParameters, ...
	modelIDToUse, (figNum));

sgtitle('Calculating the maneuver tree', 'Interpreter','none');

%%%% Plot the locations of 90-degrees
all90DegreeXYPoints = nan(length(steeringInterval),2);
for ith_steering = 1:length(steeringInterval)
	thisSteering = steeringInterval(ith_steering); % This is in radians
	thisRadius = vehicleParameters.L/tan(thisSteering); % In meters
	all90DegreeXYPoints(ith_steering,:) = [thisRadius thisRadius];
end
plot(all90DegreeXYPoints(:,1),all90DegreeXYPoints(:,2),'-o','DisplayName','Locations of 90-degree orientation');

sgtitle('Showing the parking points', 'Interpreter','none');

%%%% Highlight the 20-degree solution
% Uses fcn_VD_kinematicBicycleModelRK4

% Set up inputs
steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
steering_amplitude_radians = steering_amplitude_degrees*pi/180;
predictedTurnRadius = vehicleParameters.L/tan(steering_amplitude_radians);


% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
arcLength = 90*pi/180*predictedTurnRadius;
endTime = arcLength/vehicleParameters.U;
timeInterval = [startTime endTime];  % Units are [sec]

simulationTimes = (startTime:deltaT:endTime)';
inputsVsTime = [simulationTimes steering_amplitude_radians*ones(length(simulationTimes(:,1)),1)]; % [times steering angles]

% Call the function
[standardParkingTurnPath, ~, steeringUsed] = ...
	fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
	timeInterval, inputsVsTime, vehicleParameters, (figNum));

fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(1,:)', [steeringUsed(1,1); 0], true, (figNum));
fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(end,:)', [steeringUsed(end,1); 0], true, (figNum));

sgtitle('Highlighting the 20-degree turn', 'Interpreter','none');


%% Step 1 - Generate a standard parking turn path


dataFileName = fullfile(pwd,'Data','parkingData.mat');
if 1==1 && exist(dataFileName,'file')
	load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle'); 
else


	standard_parking_steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
	standard_parking_steering_amplitude_radians = standard_parking_steering_amplitude_degrees*pi/180;

	% Set the simulation time/state arguments
	initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
	deltaT = 0.01; % Units are [sec]
	startTime = 0;
	arcLength = 170*pi/180*predictedTurnRadius;
	endTime = arcLength/vehicleParameters.U;
	timeInterval = [startTime endTime];  % Units are [sec]

	simulationTimes = (startTime:deltaT:endTime)';
	inputsVsTime = [simulationTimes standard_parking_steering_amplitude_radians*ones(length(simulationTimes(:,1)),1)]; % [times steering angles]

    % Call a function to generate a standard parking path
    figNum = 100;
    figure(figNum);clf;
	[standardParkingTurnPath, ~, ~] = ...
		fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
		timeInterval, inputsVsTime, vehicleParameters, (figNum));

	% Call a function to convert the standard path into t- and u-distances
	[tDistancesEachAngle, uDistancesEachAngle] = fcn_VD_parkingSolveUTFromStandardTurn(...
		standardParkingTurnPath, (-1));

	save(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle');

end

%% Generate an example problem, keeping true solution
figNum = 10000;
figure(figNum);clf;

initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueEndOfArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueEndOfArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (figNum)); %#ok<ASGLU>


startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

%% Step 1: Find the x-axis intercept of the given pose
figNum = 10001;
figure(figNum);clf;


[tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(startingPointInXYphi, goalPointInXYphi, (figNum));

disp('Third set:')
disp([tDistance uDistance]);



%% Step 2: Make sure this is in the feasible range

changeInPhiInRadians = goalPointInXYphi(1,3) - startingPointInXYphi(1,3);

% Look up t-value for this goalPhiInRadians value
minimumTfeasible = interp1(standardParkingTurnPath(:,3),tDistancesEachAngle,changeInPhiInRadians,'linear','extrap');

% Look up U-value for this goalPhiInRadians value
maximumUfeasible = interp1(standardParkingTurnPath(:,3),uDistancesEachAngle,changeInPhiInRadians,'linear','extrap');

if tDistance<minimumTfeasible
	error('t-value is not feasible');
end

if uDistance>maximumUfeasible
	error('u-value is not feasible');
end

%% Step 3: Calculate the "approach" distance, the "exit" distance, and the arc angle
figNum = 10003;
figure(figNum);clf;

approachPathDistance = tDistance-minimumTfeasible;
exitPathDistance = maximumUfeasible-uDistance;
endOfArcAngleInRadians = changeInPhiInRadians;

% Plot the results


anglesToPlotTangents = [startingPointInXYphi(1,3); goalPointInXYphi(1,3)] - startingPointInXYphi(1,3);

[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), (figNum));




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

function fcn_INTERNAL_plotUandT(vehicleParameters, startPointXYphi, endPointXYphi, uDistance, tDistance, location, figNum)

endPhiInRadians = endPointXYphi(1,3);
unitEndVector = [cos(endPhiInRadians) sin(endPhiInRadians)];

% Plot scenario
h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startPointXYphi', [0; 0], true, (figNum));
set(h_start{1},'DisplayName','Start');

h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, endPointXYphi', [0; 0], true, (figNum));
set(h_end{1},'DisplayName','Goal');

% Plot orientation projection onto X
quiver(endPointXYphi(1,1),endPointXYphi(1,2),unitEndVector(1)*uDistance, unitEndVector(2)*uDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
textLocation = [endPointXYphi(1,1), endPointXYphi(1,2)] + unitEndVector*uDistance/2;
text(textLocation(1,1), textLocation(1,2),sprintf('u = %.2f meters',uDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');

% Plot x-axis distance
quiver(0,0,1*tDistance, 0*tDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
textLocation = [0, 0.5] + [1 0]*tDistance/2;
text(textLocation(1,1), textLocation(1,2),sprintf('t = %.2f meters',tDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');

% Plot hit
plot(location(1,1), location(1,2),'.','MarkerSize',50);
end


%% fcn_INTERNAL_interpolateAtAngle
function interpolatedPoint = fcn_INTERNAL_interpolateAtAngle(referenceTurnPathXYPhi,angleToInterpolate)
	% Interpolate the last X value in arc
	arcLastX = interp1(referenceTurnPathXYPhi(:,3),referenceTurnPathXYPhi(:,1),angleToInterpolate,'linear');
	arcLastY = interp1(referenceTurnPathXYPhi(:,3),referenceTurnPathXYPhi(:,2),angleToInterpolate,'linear');
	interpolatedPoint = [arcLastX arcLastY angleToInterpolate];
end % Ends fcn_INTERNAL_interpolateAtAngle
