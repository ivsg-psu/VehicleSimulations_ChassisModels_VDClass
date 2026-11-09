% script_test_fcn_VD_parkingGenerateRandomParkingSpot
% tests fcn_VD_parkingGenerateRandomParkingSpot.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_parkingGenerateRandomParkingSpot
%
% 2026_11_08 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_parkingGenerateRandomParkingSpot
%   % * First write of script
%   % * Used script_test_fcn_VD_parking+CalcEntryArcExitPaths as starter

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

rng(0);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [];
roundingLevel = [];

% Call the function
[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
	startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(goalParkingPointInXYphi));
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
% Too difficult to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, different RNG
figNum = 10002;
titleString = sprintf('DEMO case: basic call, different RNG');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

rng(1);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [];
roundingLevel = [];

% Call the function
[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
	startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(goalParkingPointInXYphi));
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
% Too difficult to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: varying input range
figNum = 10003;
titleString = sprintf('DEMO case: varying input range');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

rng(0);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [1 10 45*pi/180 50*pi/180 1 10];
roundingLevel = [];

% Call the function
[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
	startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(goalParkingPointInXYphi));
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
% Too difficult to check

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

rng(0);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [];
roundingLevel = [];

% Call the function
[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
	startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(goalParkingPointInXYphi));
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
% Too difficult to check

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

rng(0);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [];
roundingLevel = [];

% Call the function
[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
	startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(goalParkingPointInXYphi));
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
% Too difficult to check

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

rng(0);

% Set up inputs
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if exist(dataFileName,'file')
    load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle')
else
	error('Need to run fcn_VD_parkingFillStandardParkingTurnPath to generate the data file');
end

startingPointInXYphi = [0 0 0];
flagDifficulty = [];
inputRange = [];
roundingLevel = [];


Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations

	% Call the function
	[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
		standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
		startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations

	% Call the function
	[goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
		standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
		startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (-1));

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
