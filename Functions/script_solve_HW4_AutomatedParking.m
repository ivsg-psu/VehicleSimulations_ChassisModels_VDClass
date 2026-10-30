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
% [standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle] = ...
% 	fcn_VD_parkingFillStandardParkingTurnPath((-1));



%% Generate an example problem, keeping true solution
figNum = 10000;
figure(figNum);clf;

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run a prior section to generate the data file');
end

% Generate a true solution
initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, ~, exitXYPsi, hPlots] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (figNum));

for ith_plot = 1:3
    set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

% Approximately a Ford Taurus
clear vehicleParameters
vehicleParameters(1).a          = 0.9271; % Distance from front axle to CG, in meters
vehicleParameters(1).b          = 1.5621;  % Distance from rear axle to CG, in meters

fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi', [0; 0], true, (figNum));
fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalPointInXYphi', (true), (figNum));


% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));
assert(isnumeric(approachPathDistance));
assert(isnumeric(changeInArcAngleInRadians));
assert(isnumeric(exitPathDistance));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);
assert(size(approachPathDistance,1)==1);
assert(size(approachPathDistance,2)==1);
assert(size(changeInArcAngleInRadians,1)==1);
assert(size(changeInArcAngleInRadians,2)==1);
assert(size(exitPathDistance,1)==1);
assert(size(exitPathDistance,2)==1);

% Check variable values
assert(isequal(round(trueApproachPathDistance,4),round(approachPathDistance,4)));
assert(isequal(round(trueChangeInArcAngleInRadians,4),round(changeInArcAngleInRadians,4)));
assert(isequal(round(trueExitPathDistance,4),round(exitPathDistance,4)));

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

