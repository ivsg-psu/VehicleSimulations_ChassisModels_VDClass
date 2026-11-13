% script_test_fcn_VD_steeringFillTestCase.m
% tests fcn_VD_steeringFillTestCase.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_steeringFillTestCase
% 
% 2026_11_12 by Sean Brennan, sbrennan@psu.ed
% - In script_test_fcn_VD_steeringFillTestCase
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

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, specifying amplitude
figNum = 10002;
titleString = sprintf('DEMO case: basic call, specifying amplitude');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = [];
testCaseOptions = cell(1,1);
testCaseOptions{1} = 10; % Amplitude in degrees
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, specifying simulationTime
figNum = 10003;
titleString = sprintf('DEMO case: basic call, specifying simulationTime');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = [];
testCaseOptions = [];
simulationTimes = (1:0.4:7)';

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

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

%% TEST case: mode 1 (single sine wave) with default inputs
figNum = 20011;
titleString = sprintf('TEST case: mode 1 (single sine wave) with default inputs');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 1;
testCaseOptions = [];
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 1 (single sine wave) with sinewave input, specifying amplitude, period, and start time
figNum = 20012;
titleString = sprintf('TEST case: mode 1 (single sine wave) with sinewave input, specifying amplitude, period, and start time');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 1;
testCaseOptions = cell(3,1);
testCaseOptions{1} = 10; % Amplitude in degrees
testCaseOptions{2} = 2; % Period in seconds
testCaseOptions{3} = 0.5; % Start time in seconds
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 1 (single sine wave) with sinewave input, specifying simulationTime
figNum = 20013;
titleString = sprintf('TEST case: mode 1 (single sine wave) with sinewave input, specifying simulationTime');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 1;
testCaseOptions = [];
simulationTimes = (0.5:0.4:7)';

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 2 (pseudo-random binary sequence) with default inputs
figNum = 20021;
titleString = sprintf('TEST case: mode 2 (pseudo-random binary sequence) with default inputs');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 2;
testCaseOptions = [];
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 2 (pseudo-random binary sequence) with sinewave input, specifying signal settings
figNum = 20022;
titleString = sprintf('TEST case: mode 2 (pseudo-random binary sequence) with sinewave input, specifying signal settings');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 2;
testCaseOptions = cell(4,1);
%				testCaseOptions{1}: maxSteeringAmplitudeDegrees - the
%				maximum degrees of steering amplitude for input sinewave.
%				Default is 5 degrees.
%
%				testCaseOptions{2}: minSteeringAmplitudeDegrees - the
%				minimum degrees of steering amplitude for input sinewave.
%				Default is -5 degrees.
%
%				testCaseOptions{3}: minDwellTime - the minimum time within
%				a steering state, units are seconds. Default is 0.1
%				seconds.
%
%				testCaseOptions{4}: maxDwellTime - the maximum time within
%				a steering state, units are seconds. Default is 1
%				seconds.
testCaseOptions{1} = 10; % Max amplitude in degrees
testCaseOptions{2} = -20; % Min amplitude in degrees
testCaseOptions{3} = 0.02; % Min dwell time in seconds
testCaseOptions{4} = 0.4; % Max dwell time in seconds

simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 2 (pseudo-random binary sequence) with sinewave input, specifying simulationTime
figNum = 20023;
titleString = sprintf('TEST case: mode 2 (pseudo-random binary sequence) with sinewave input, specifying simulationTime');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 2;
testCaseOptions = [];
simulationTimes = (0.5:0.4:7)';

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% TEST case: mode 3 (continuous sine wave) with default inputs
figNum = 20031;
titleString = sprintf('TEST case: mode 3 (continuous sine wave) with default inputs');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 3;
testCaseOptions = [];
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 3 (continuous sine wave) with sinewave input, specifying amplitude and period
figNum = 20032;
titleString = sprintf('TEST case: mode 3 (continuous sine wave) with sinewave input, specifying amplitude and period');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 3;
testCaseOptions = cell(2,1);
testCaseOptions{1} = 10; % Amplitude in degrees
testCaseOptions{2} = 2; % Period in seconds
simulationTimes = [];

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% TEST case: mode 3 (continuous sine wave) with sinewave input, specifying simulationTime
figNum = 20033;
titleString = sprintf('TEST case: mode 3 (continuous sine wave) with sinewave input, specifying simulationTime');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the input arguments
testCaseID = 3;
testCaseOptions = [];
simulationTimes = (1:0.4:7)';

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

% Check variable values
% Too many to check

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

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

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

% Call the function
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(steeringInputsVsTime));

% Check variable sizes
assert(size(steeringInputsVsTime,1)>=1); 
assert(size(steeringInputsVsTime,2)==2); 

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

Niterations = 500;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations

	% Call the function
	steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations

	% Call the function
	steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
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
