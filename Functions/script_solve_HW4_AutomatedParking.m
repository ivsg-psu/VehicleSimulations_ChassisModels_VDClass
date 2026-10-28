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
figNum = 10001;
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

%% Plot the locations of 90-degrees
all90DegreeXYPoints = nan(length(steeringInterval),2);
for ith_steering = 1:length(steeringInterval)
	thisSteering = steeringInterval(ith_steering); % This is in radians
	thisRadius = vehicleParameters.L/tan(thisSteering); % In meters
	all90DegreeXYPoints(ith_steering,:) = [thisRadius thisRadius];
end
plot(all90DegreeXYPoints(:,1),all90DegreeXYPoints(:,2),'-o','DisplayName','Locations of 90-degree orientation');

sgtitle('Showing the parking points', 'Interpreter','none');

%% Plot the 20-degree solution
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

%% Plot ONLY the 20-degree solution
figNum = 10002;

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

% Call the function to just plot it up to 90 degrees
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
	timeInterval, inputsVsTime, vehicleParameters, (figNum));

sgtitle('Only the 20-degree turn', 'Interpreter','none');


%% Plot approach and exit
approachPath = [-10 0; 0 0];
exitPath = [predictedTurnRadius predictedTurnRadius; predictedTurnRadius predictedTurnRadius+10];
plot(approachPath(:,1),approachPath(:,2),'-','Color',[0 1 0],'LineWidth',5,'DisplayName','Approach Path')
plot(exitPath(:,1),exitPath(:,2),'-','Color',0.7*[0 1 0],'LineWidth',5,'DisplayName','Exit Path')

if 1==0
	fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(1,:)', [steeringUsed(1,1); 0], true, (figNum));
	fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(end,:)', [steeringUsed(end,1); 0], true, (figNum));
end

sgtitle('All extensions of the 20-degree turn', 'Interpreter','none');


%% Plot tangent
tangentIndex = find(standardParkingTurnPath(:,3)>(40*pi/180),1); % Find first location greater than 40 degrees
tangentStates = standardParkingTurnPath(tangentIndex,:);
angleRad = tangentStates(1,3);
rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
tangentPath = [0 0; 20 0]*rotation + [tangentStates(1,1)*ones(2,1) tangentStates(1,2)*ones(2,1)] ;

plot(tangentPath(:,1),tangentPath(:,2),'-','Color',0.7*[1 1 0],'LineWidth',5,'DisplayName','Tangent Path')

sgtitle('All extensions of the 20-degree turn, with tangent shown', 'Interpreter','none');


%% Plot Vehicle "given"
givenPath = [0 0; 3 0]*rotation + [tangentStates(1,1)*ones(2,1) tangentStates(1,2)*ones(2,1)];

givenPosition = [givenPath(end,:)'; angleRad];
fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, givenPosition, [0; 0], true, (figNum));

sgtitle('Given', 'Interpreter','none');


%% Plot Vehicle "given"
givenPath = [0 0; 3 0]*rotation + [tangentStates(1,1)*ones(2,1) tangentStates(1,2)*ones(2,1)];

givenPosition = [givenPath(end,:)'; angleRad];
fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, givenPosition, [0; 0], true, (figNum));

sgtitle('Given', 'Interpreter','none');

%% Example case
figNum = 10003;
figure(figNum);clf;

if 1==0
	startingPointInxy = [0 0 0];
	goalPointInxy     = [12 7 40*pi/180];

	offset = [4 6];
	angleRad = 30*pi/180;
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	newPoints = [startingPointInxy(:,1:2); goalPointInxy(:,1:2)]*rotation + [offset(1,1)*ones(2,1) offset(1,2)*ones(2,1)] ;

	startingPointInXYphi = [newPoints(1,:) angleRad]; %#ok<NASGU>
	goalPointInXYphi = [newPoints(2,:) angleRad+goalPointInxy(1,3)]; %#ok<NASGU>

end

goalPointInXYphi = [11 18 70*pi/180];
startingPointInXYphi = [4 6 30*pi/180];

% Plot scenario
h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi', [0; 0], true, (figNum));
set(h_start{1},'DisplayName','Start');

h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, goalPointInXYphi', [0; 0], true, (figNum));
set(h_end{1},'DisplayName','Goal');

legend('Interpreter','none');
axis([0 20 5 25]);

%% Step 1: Convert the problem into the vehicle's local xy coordinates
figNum = 10004;
figure(figNum);clf;

angleRad = -1*startingPointInXYphi(1,3);
rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];

pointsToMoveXYphi = [startingPointInXYphi; goalPointInXYphi];
pointsToMoveXY = pointsToMoveXYphi(:,1:2);
movedInXY = pointsToMoveXY - ones(2,1)*startingPointInXYphi(1,1:2);

movedRotatedInXY = movedInXY*rotation;
movedStartingPointInXYphi = [movedRotatedInXY(1,:) startingPointInXYphi(1,3)+angleRad];
movedGoalPointInXYphi = [movedRotatedInXY(2,:) goalPointInXYphi(1,3)+angleRad];

% Plot scenario
h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, movedStartingPointInXYphi', [0; 0], true, (figNum));
set(h_start{1},'DisplayName','Start');

h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, movedGoalPointInXYphi', [0; 0], true, (figNum));
set(h_end{1},'DisplayName','Goal');

legend('Interpreter','none');
axis([-6 16 -6 12]);

%% Step 2: Find the x-axis intercept of the given pose
figNum = 10005;
figure(figNum);clf;

% Define a vector based on the goal point phi angle
goalPhiInRadians = movedGoalPointInXYphi(1,3);
unitGoalVector = [cos(goalPhiInRadians) sin(goalPhiInRadians)];
goalSegment = [movedGoalPointInXYphi(1,1:2); movedGoalPointInXYphi(1,1:2)+unitGoalVector];

xAxisSegment = [0 0; 1 0]; % Unit vector in x-axis direction

% Use Path library to find where this segment intersects the x-axis
% FORMAT:
%      [distance, location, path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(figNum))
flagSearchType = 4; % any intersection of the sensor and path lines
[~, location, ~, tDistance, uDistance] = ...
	fcn_Path_findProjectionHitOntoPath(xAxisSegment,goalSegment(1,:),goalSegment(2,:),...
	flagSearchType,-1);

fcn_INTERNAL_plotUandT(vehicleParameters, ...
	movedStartingPointInXYphi, movedGoalPointInXYphi, ...
	uDistance, tDistance, location, figNum);


% % Plot scenario
% h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, movedStartingPointInXYphi', [0; 0], true, (figNum));
% set(h_start{1},'DisplayName','Start');
% 
% h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, movedGoalPointInXYphi', [0; 0], true, (figNum));
% set(h_end{1},'DisplayName','Goal');
% 
% % Plot orientation projection onto X
% quiver(movedGoalPointInXYphi(1,1),movedGoalPointInXYphi(1,2),unitGoalVector(1)*uDistance, unitGoalVector(2)*uDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
% textLocation = [movedGoalPointInXYphi(1,1), movedGoalPointInXYphi(1,2)] + unitGoalVector*uDistance/2;
% text(textLocation(1,1), textLocation(1,2),sprintf('u = %.2f meters',uDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');
% 
% % Plot x-axis distance
% quiver(0,0,1*tDistance, 0*tDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
% textLocation = [0, 0.5] + [1 0]*tDistance/2;
% text(textLocation(1,1), textLocation(1,2),sprintf('t = %.2f meters',tDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');
% 
% % Plot hit
% plot(location(1,1), location(1,2),'.','MarkerSize',50);



%% Step 2.5 - Generate reference x-location curve per each angle
figNum = 10006;
figure(figNum);clf;

dataFileName = fullfile(pwd,'Data','parkingData.mat');
if 1==0 && exist(dataFileName,'file')
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

	% Call a function
	[standardParkingTurnPath, ~, ~] = ...
		fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
		timeInterval, inputsVsTime, vehicleParameters, (figNum));

	% Call the function
	[tDistancesEachAngle, uDistancesEachAngle] = fcn_VD_parkingSolveUTFromStandardTurn(...
		standardParkingTurnPath, (-1));

	save(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle');

end
subplot(1,3,1);
% Call the function just to plot
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
	timeInterval, inputsVsTime, vehicleParameters, (figNum));

h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(1,:)', [standard_parking_steering_amplitude_radians; 0], true, (figNum));
set(h_start{1},'DisplayName','Start');

h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(end,:)', [standard_parking_steering_amplitude_radians; 0], true, (figNum));
set(h_end{1},'DisplayName','Goal');

subplot(1,3,2);
plot(standardParkingTurnPath(:,3)*180/pi, tDistancesEachAngle,'.-');
xlabel('Angle [deg]');
ylabel('t-distances [m]');

subplot(1,3,3);
plot(standardParkingTurnPath(:,3)*180/pi, uDistancesEachAngle,'.-');
hold on;
plot(standardParkingTurnPath(:,3)*180/pi, -1*tDistancesEachAngle,'.-');
xlabel('Angle [deg]');
ylabel('u-distances [m]');


sgtitle('Calculating x-axis projection intercepts (t) and distances (u) for the standard turn', 'Interpreter','none');

%% Step 3: Make sure this is in the feasible range

dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
	load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run a prior section to generate the data file');
end

% Define a vector based on the goal point phi angle
goalPhiInRadians = movedGoalPointInXYphi(1,3);
unitGoalVector = [cos(goalPhiInRadians) sin(goalPhiInRadians)];
goalSegment = [movedGoalPointInXYphi(1,1:2); movedGoalPointInXYphi(1,1:2)+unitGoalVector];

xAxisSegment = [0 0; 1 0]; % Unit vector in x-axis direction

% Use Path library to find where this segment intersects the x-axis
% FORMAT:
%      [distance, location, path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(figNum))
flagSearchType = 4; % any intersection of the sensor and path lines
[distance, location, ~, tDistance, uDistance] = ...
	fcn_Path_findProjectionHitOntoPath(xAxisSegment,goalSegment(1,:),goalSegment(2,:),...
	flagSearchType,-1);

% Look up t-value for this goalPhiInRadians value
minimumTfeasible = interp1(standardParkingTurnPath(:,3),tDistancesEachAngle,goalPhiInRadians,'linear','extrap');

% Look up U-value for this goalPhiInRadians value
maximumUfeasible = interp1(standardParkingTurnPath(:,3),uDistancesEachAngle,goalPhiInRadians,'linear','extrap');

if tDistance<minimumTfeasible
	error('t-value is not feasible');
end

if uDistance>maximumUfeasible
	error('u-value is not feasible');
end

%% Step 4: Calculate the "approach" distance, the "exit" distance, and the arc angle
approachPathDistance = tDistance-minimumTfeasible;
exitPathDistance = maximumUfeasible-uDistance;
endOfArcAngleInRadians = goalPhiInRadians;

% Plot the results

figNum = 10007;
figure(figNum);clf;
anglesToPlotTangents = [standardParkingStartPoint(1,3);standardParkingStartPoint(1,3)+endOfArcAngleInRadians];

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
