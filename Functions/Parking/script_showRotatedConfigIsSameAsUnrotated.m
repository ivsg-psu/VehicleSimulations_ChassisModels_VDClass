figNum = 10001;
figure(figNum);clf;


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

%% Plot 20 degree curve
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
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
	load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run a prior section to generate the data file');
end

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


%% Example case
figNum = 10003;
figure(figNum);clf;

if 1==0
    % Take a given solution, and rotate everything by 30 degrees and offset
    % it by [4 6] in XY

    % Take a good solution
	startingPointInxy = [0 0 0];
	goalPointInxy     = [12 7 40*pi/180];

    % Apply rotation and offset
	offset = [4 6];
	angleRad = 30*pi/180;
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	newPoints = [startingPointInxy(:,1:2); goalPointInxy(:,1:2)]*rotation + [offset(1,1)*ones(2,1) offset(1,2)*ones(2,1)] ;

	startingPointInXYphi = [newPoints(1,:) angleRad]; %#ok<NASGU>
	goalPointInXYphi = [newPoints(2,:) angleRad+goalPointInxy(1,3)]; %#ok<NASGU>

end

% Below are close to above values
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

disp('First set:')
disp([tDistance uDistance]);

fcn_INTERNAL_plotUandT(vehicleParameters, ...
	movedStartingPointInXYphi, movedGoalPointInXYphi, ...
	uDistance, tDistance, location, figNum);

[tDistance, uDistance] = fcn_VD_parkingXYPhiToTU([0 0 0], movedGoalPointInXYphi, (4444));

disp('Second set:')
disp([tDistance uDistance]);


[tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(startingPointInXYphi, goalPointInXYphi, (5555));

disp('Third set:')
disp([tDistance uDistance]);