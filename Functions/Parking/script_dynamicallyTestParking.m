close all;

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


% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

initialXYPhi = [-4 -6 -30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = -40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; -40*pi/180];

[approachXYPsi, curveXYPsi, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (222));

currentSolutionXYPsi = [approachXYPsi; nan nan nan; curveXYPsi; nan nan nan; exitXYPsi];

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);




% Create figure and attach drag callbacks
fig = figure('Name','Draggable Line Segments', ...
    'WindowButtonDownFcn',@startDrag, ...
    'WindowButtonMotionFcn',@drag, ...
    'WindowButtonUpFcn',@stopDrag);


figNum = fig.Number;
[~, vehicleLocationPoints] = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], (true), (figNum));
[~, parkingSpotPoints] = fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalPointInXYphi, (true), (figNum));


grabPointsVehicle = [vehicleLocationPoints(1,1:2); vehicleLocationPoints(2,1:2)];
grabPointsSpot = [parkingSpotPoints(1,1:2); parkingSpotPoints(2,1:2)];



% Set up axes for interactive display
ax = gca;
hold(ax,'on');
axis(ax,'equal');
grid(ax,'on');
axis([-20 20 -20 20]);

% Plot the two movable segments and the parking solution
hVehicle = plot(ax,grabPointsVehicle(:,1),grabPointsVehicle(:,2),'b-','LineWidth',6);
hSpot = plot(ax,grabPointsSpot(:,1),grabPointsSpot(:,2),'r-','LineWidth',6);
hSolution = plot(ax,currentSolutionXYPsi(:,1),currentSolutionXYPsi(:,2), 'b-','LineWidth',1);

% Create hidden distance label at the connector midpoint
ht = text(ax,mean([grabPointsVehicle(2,1) grabPointsSpot(2,1)]),mean([grabPointsVehicle(2,2) grabPointsSpot(2,2)]),'', ...
    'HorizontalAlignment','center','BackgroundColor','w', ...
    'Visible','off');

% Store handles and drag state in figure app data
d.p = {grabPointsVehicle,grabPointsSpot};
d.hPlots = [hVehicle hSpot hSolution];
d.parkingSolution = hSolution;
d.label = ht;
d.active = [];
d.mode = [];
d.offset = [];
d.grabPointsVehicle = grabPointsVehicle;
d.grabPointsSpot = grabPointsSpot;
d.vehicleLocationPoints = vehicleLocationPoints;
d.parkingSpotPoints = parkingSpotPoints;
d.startingPointInXYphi = startingPointInXYphi;
d.goalPointInXYphi = goalPointInXYphi;
d.currentSolutionXYPsi = currentSolutionXYPsi;
guidata(fig,d);

%% startDrag
function startDrag(src,~)
d = guidata(src);
q = src.CurrentAxes.CurrentPoint(1,1:2);

% Find the nearest segment endpoint under the cursor
ends = [d.p{1}; d.p{2}];
[dist,idx] = min(vecnorm(ends-q,2,2));

% Ignore clicks that are too far from any endpoint
threshold = 1;
if dist > threshold
    return
end

% Select the target segment and whether its start or end is moving
d.active = ceil(idx/2);
d.mode = mod(idx-1,2)+1;

% Preserve the grabbed-point offset when dragging whole segments
if d.mode == 1
    % Update the active segment offset from the current query point.
    d.offset = q - d.p{d.active}(1,:);
else
    d.offset = [];
end

guidata(src,d);
end % Ends startDrag


%% drag
function drag(src,~)
d = guidata(src);

% Exit early unless a drag is currently active
if isempty(d.active)
    return
end

q = src.CurrentAxes.CurrentPoint(1,1:2);
p = d.p{d.active};

% Update either the whole segment or just its free endpoint
if d.mode == 1
    delta = q - d.offset - p(1,:);
    p = p + delta;
else
    p(2,:) = q;
end

% Refresh segment, parkingSolution, and distance label

% Update the vehicle
d.p{d.active} = p;
set(d.lines(d.active),'XData',p(:,1),'YData',p(:,2));

% Update the spot
a = d.p{1}(2,:);
b = d.p{2}(2,:);
set(d.parkingSolution,'XData',[a(1) b(1)],'YData',[a(2) b(2)]);

% Update the text box
m = (a+b)/2;
dist = norm(b-a);
set(d.label,'Position',[m 0],'String',sprintf('%.2f',dist),'Visible','on');

guidata(src,d);
drawnow limitrate
end % Ends drag

%% stopDrag
function stopDrag(src,~)
d = guidata(src);
d.active = [];
d.mode = [];
d.offset = [];
guidata(src,d);
end % Ends stopDrag