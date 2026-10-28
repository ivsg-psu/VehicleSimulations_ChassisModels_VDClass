% script_test_fcn_VD_parkingPlotSolution
% tests fcn_VD_parkingPlotSolution.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_parkingPlotSolution
%
% 2026_10_28 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_parkingPlotSolution
%   % * First write of script
%   % * Used script_test_fcn_VD_plotStickFigure2DOFBicycleModel as starter

% TO-DO:
%
% 2026_10_28 by Sean Brennan, sbrennan@psu.edu
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
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [0 0 0];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 45*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = []; % Default

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), (figNum));
legend('Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));



%% DEMO case: basic call, with nonzero ICs
figNum = 10002;
titleString = sprintf('DEMO case: basic call, with nonzero ICs');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [-7 11 240*pi/180];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 45*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = []; % Default

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), (figNum));
legend('Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, nonzero ICs with tangents
figNum = 10003;
titleString = sprintf('DEMO case: basic call, nonzero ICs with tangents');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [-7 11 240];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 95*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = [0; 40*pi/180; endOfArcAngleInRadians];

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), (figNum));
legend('Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);

% Check variable values
% (too complex to check)

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
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [0 0 0];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 45*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = []; % Default

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), ([]));
legend('Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);

% Check variable values
% (too complex to check)

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
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [0 0 0];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 45*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = []; % Default

% Call the function
[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
	(anglesToPlotTangents), (-1));
legend('Interpreter','none');

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(approachXYPsi));
assert(isnumeric(curveXYPsi));
assert(isnumeric(exitXYPsi));

% Check variable sizes
assert(size(approachXYPsi,1)==2);
assert(size(approachXYPsi,2)==3);
assert(size(curveXYPsi,1)>=2);
assert(size(curveXYPsi,2)==3);
assert(size(exitXYPsi,1)==2);
assert(size(exitXYPsi,2)==3);

% Check variable values
% (too complex to check)

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
	% load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
	load(dataFileName,'standardParkingTurnPath')
else
	error('Need to run a prior section to generate the data file');
end

initialXYPhi = [0 0 0];
approachPathDistance = 2; % meters
endOfArcAngleInRadians = 45*pi/180; % radians
exitPathDistance = 5; % meters
anglesToPlotTangents = []; % Default



Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
	% Call the function
	[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
		standardParkingTurnPath, initialXYPhi, ...
		approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
		(anglesToPlotTangents), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
	% Call the function
	[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
		standardParkingTurnPath, initialXYPhi, ...
		approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
		(anglesToPlotTangents), (-1));
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
