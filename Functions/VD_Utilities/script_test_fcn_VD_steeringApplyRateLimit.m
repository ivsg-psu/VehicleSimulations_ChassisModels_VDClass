% script_test_fcn_VD_steeringApplyRateLimit.m
% tests fcn_VD_steeringApplyRateLimit.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_steeringApplyRateLimit
%
% 2026_11_12 by Sean Brennan, sbrennan@psu.ed
% - In script_test_fcn_VD_steeringApplyRateLimit
%   % * Wrote the code originally,
%   % * Used script_test_fcn_VD_st2dofDerivatives as starter

% TO-DO:
%
% 2026_11_12 by Sean Brennan, sbrennan@psu.edu
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

%% DEMO case: basic call with default inputs
figNum = 10001;
titleString = sprintf('DEMO case: basic call with default inputs');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = [];

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% DEMO case: basic call, changing rate limit
figNum = 10002;
titleString = sprintf('DEMO case: basic call, changing rate limit');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = 9; % Units are deg/sec

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, uneven time sampling
figNum = 10003;
titleString = sprintf('DEMO case: basic call, uneven time sampling');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 2;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = 20; % Units are deg/sec

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
% Too many to check

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

%% TEST case: making sure ramp matches expected
figNum = 20001;
titleString = sprintf('TEST case: making sure ramp matches expected');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments by creating a steering input so large that it
% has to be rate limited the entire trajectory
timeSteering = (0:0.01:1)';
Nsteering = length(timeSteering);
steeringAngles = 20*ones(Nsteering,1)*pi/180;
steeringAngles(1,1) = 0;
steeringInputsVsTime = [timeSteering steeringAngles];

rateLimitInDegPerSec = 10; % Units are deg/sec

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
assert(isequal(round(rateLimitedSteeringInputsVsTime(end,2)*180/pi,6), round(rateLimitInDegPerSec,6)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% TEST case: making sure ramp matches expected (downward)
figNum = 20002;
titleString = sprintf('TEST case: making sure ramp matches expected (downward)');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments by creating a steering input so large that it
% has to be rate limited the entire trajectory
timeSteering = (0:0.01:1)';
Nsteering = length(timeSteering);
steeringAngles = -20*ones(Nsteering,1)*pi/180;
steeringAngles(1,1) = 0;
steeringInputsVsTime = [timeSteering steeringAngles];

rateLimitInDegPerSec = 10; % Units are deg/sec

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
assert(isequal(-1*round(rateLimitedSteeringInputsVsTime(end,2)*180/pi,6), round(rateLimitInDegPerSec,6)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


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

% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = [];

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
% Too many to check

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = [];

% Call the function
rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rateLimitedSteeringInputsVsTime));

% Check variable sizes
assert(size(rateLimitedSteeringInputsVsTime,1)>=1);
assert(size(rateLimitedSteeringInputsVsTime,2)==2);

% Check variable values
% Too many to check

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);


% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

rateLimitInDegPerSec = [];

Niterations = 500;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
	% Call the function
	rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
	% Call the function
	rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));
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
