% script_test_fcn_VD_parkingCalcEntryArcExitPaths
% tests fcn_VD_parkingCalcEntryArcExitPaths.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_parkingCalcEntryArcExitPaths
%
% 2026_10_29 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_parkingCalcEntryArcExitPaths
%   % * First write of script
%   % * Used script_test_fcn_VD_parking+XYPhiToTU as starter
%
% 2026_10_30 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_parkingCalcEntryArcExitPaths
%   % * Created test cases for negative solutions

% TO-DO:
%
% 2026_10_29 by Sean Brennan, sbrennan@psu.edu
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

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
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

fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));
fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalPointInXYphi, (true), (figNum));


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

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% DEMO case: basic call, only plotting core function
figNum = 10002;
titleString = sprintf('DEMO case: basic call, only plotting core function');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, ~, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (-1));

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

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

fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));
fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalPointInXYphi, (true), (figNum));


% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, (figNum));


fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));


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

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: negative angle
figNum = 10003;
titleString = sprintf('DEMO case: negative angle');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

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

[approachXYPsi, ~, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (-1));

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

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

fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));
fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalPointInXYphi, (true), (figNum));


% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, (figNum));


fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));


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

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, ~, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (-1));

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, ([]));

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, ~, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (-1));

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, (-1));

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

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [4 6 30*pi/180];
trueApproachPathDistance = 10;
trueChangeInArcAngleInRadians = 40*pi/180;
trueExitPathDistance = 7;
anglesToPlotTangents = [0; 40*pi/180];

[approachXYPsi, ~, exitXYPsi, hPlots] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	trueApproachPathDistance, trueChangeInArcAngleInRadians, trueExitPathDistance, ...
	(anglesToPlotTangents), (-1));

% for ith_plot = 1:3
%     set(hPlots{ith_plot,1},'LineWidth',10,'Color',[0 0 1]);
% end

startingPointInXYphi = approachXYPsi(1,:);
goalPointInXYphi = exitXYPsi(end,:);


Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations

    % Call the function
    [approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
        standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
        startingPointInXYphi, goalPointInXYphi, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations

    % Call the function
    [approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
        standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
        startingPointInXYphi, goalPointInXYphi, (-1));

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
