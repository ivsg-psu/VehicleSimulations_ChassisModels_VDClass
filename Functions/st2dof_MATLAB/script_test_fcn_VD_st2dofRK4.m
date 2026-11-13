% script_test_fcn_VD_st2dofRK4.m
% tests fcn_VD_st2dofRK4.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_2DOFBicycleModelRK4
%
% 2026_10_20 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_2DOFBicycleModelRK4
%   % * Used script_test_fcn_VD_kinematicBicycleModelRK4 as starter
%   % * Copied material over from script_test_fcn_VD_bicycle2dofModel
%
% As: script_test_fcn_VD_st2dofRK4
%
% 2026_11_12 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_st2dofRK4
%   % * Renamed script to be consistent with single-track (bicycle) naming
%   % * Used to be called script_test_fcn_VD_2DOFBicycleModelRK4
%


% TO-DO:
%
% 2026_10_20 by Sean Brennan, sbrennan@psu.edu
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

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
rng(0);
testCaseID = 0;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
rateLimitInDegPerSec = [];
inputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

% Call the function
[stateTrajectory, t, steeringUsed] = ...
	fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectory));
assert(isnumeric(t));
assert(isnumeric(steeringUsed));

% Check variable sizes
assert(size(stateTrajectory,1)>=1);
assert(size(stateTrajectory,2)==5);
assert(size(t,1)==size(stateTrajectory,1));
assert(size(t,2)==1);
assert(size(steeringUsed,1)==size(stateTrajectory,1));
assert(size(steeringUsed,2)==1);

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


%% TEST case: compare to kinematic
figNum = 20001;
titleString = sprintf('DEMO case: compare to kinematic');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
rng(0);
testCaseID = 0;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
rateLimitInDegPerSec = [];
inputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

% Call the function
[stateTrajectory, t, steeringUsed] = ...
	fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectory));
assert(isnumeric(t));
assert(isnumeric(steeringUsed));

% Check variable sizes
assert(size(stateTrajectory,1)>=1);
assert(size(stateTrajectory,2)==5);
assert(size(t,1)==size(stateTrajectory,1));
assert(size(t,2)==1);
assert(size(steeringUsed,1)==size(stateTrajectory,1));
assert(size(steeringUsed,2)==1);

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


% Call the function
initialStates = [0 0 0];
vehicleParameters.U = U;
vehicleParameters.L = vehicleParameters.DistancesWheelbase;

[stateTrajectory, t, steeringUsed] = ...
fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
timeInterval, inputsVsTime, vehicleParameters, (figNum)); %#ok<ASGLU>

%% DEMO case: PRBS input
figNum = 20002;
titleString = sprintf('DEMO case: basic call');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.001; % Units are [sec]
startTime = 0;
endTime = 10;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
testCaseID = 2;
testCaseOptions = cell(4,1);
testCaseOptions{1} = 3; % Max amplitude in degrees
testCaseOptions{2} = -3; % Min amplitude in degrees
testCaseOptions{3} = 0.5; % Min dwell time in seconds
testCaseOptions{4} = 2; % Max dwell time in seconds
simulationTimes = (startTime:deltaT:endTime)';

% Call the function
inputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

% Call the function
[stateTrajectory, t, steeringUsed] = ...
	fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, (-1));

% Add noise to outputs
noiseOptions = cell(4,1);
noiseTypes = [1 2 3];
noiseOptions{1} = 0;
noiseOptions{2} = 0.01;
noiseOptions{3} = [0.95 0.005];

% Call the function
noisyStateTrajectory = stateTrajectory;
noisyStateTrajectory(:,1) = fcn_VD_signalAddNoise(stateTrajectory(:,1), (noiseTypes), (noiseOptions), (figNum));
noisyStateTrajectory(:,2) = fcn_VD_signalAddNoise(stateTrajectory(:,2), (noiseTypes), (noiseOptions), (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectory));
assert(isnumeric(t));
assert(isnumeric(steeringUsed));

% Check variable sizes
assert(size(stateTrajectory,1)>=1);
assert(size(stateTrajectory,2)==5);
assert(size(t,1)==size(stateTrajectory,1));
assert(size(t,2)==1);
assert(size(steeringUsed,1)==size(stateTrajectory,1));
assert(size(steeringUsed,2)==1);

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

figure(figNum);
subplot(2,1,1);
plot(t,steeringUsed(:,1),'-','DisplayName','Steering input, road wheel angle (rad)');
legend('Interpreter','none');
ylabel('Steering input [rad]','Interpreter','none');

subplot(2,1,2);
hold on;

randomNumbering = randperm(4, 4);
dataToPlot = [noisyStateTrajectory(:,1:2) stateTrajectory(:,1:2)];
noisyData = [1 1 0 0];

for ith_handle = 1:4
	thisIndex = randomNumbering(ith_handle);
	thisData = dataToPlot(:,thisIndex);
	thisFlag = noisyData(1,thisIndex);

	thisHandle = plot(t,thisData,'-','DisplayName',sprintf('%.0d',ith_handle),'LineWidth',3);
	if 1==thisFlag
		thisColor = get(thisHandle,'Color');
		set(thisHandle,'Color',[thisColor 0.1]);
	end
end

% legend('show');
hLegend = legend('Interpreter','none');
% set(hLegend,'String',{'1' '2' '3' '4'});

xlabel('Time [sec]','Interpreter','none');

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

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
rng(0);
testCaseID = 0;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
rateLimitInDegPerSec = [];
inputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

% Call the function
[stateTrajectory, t, steeringUsed] = ...
	fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectory));
assert(isnumeric(t));
assert(isnumeric(steeringUsed));

% Check variable sizes
assert(size(stateTrajectory,1)>=1);
assert(size(stateTrajectory,2)==5);
assert(size(t,1)==size(stateTrajectory,1));
assert(size(t,2)==1);
assert(size(steeringUsed,1)==size(stateTrajectory,1));
assert(size(steeringUsed,2)==1);

% Check variable values
% (too complex to check)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
rng(0);
testCaseID = 0;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
rateLimitInDegPerSec = [];
inputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

% Call the function
[stateTrajectory, t, steeringUsed] = ...
	fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectory));
assert(isnumeric(t));
assert(isnumeric(steeringUsed));

% Check variable sizes
assert(size(stateTrajectory,1)>=1);
assert(size(stateTrajectory,2)==5);
assert(size(t,1)==size(stateTrajectory,1));
assert(size(t,2)==1);
assert(size(steeringUsed,1)==size(stateTrajectory,1));
assert(size(steeringUsed,2)==1);

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

% Set the simulation time/state arguments
initialStates = [0 0 0 0 0]; % [V r X Y phi] in [m/s], [rad/sec], [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set the input arguments
rng(0);
testCaseID = 0;
testCaseOptions = [];
simulationTimes = [];
steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (-1));
rateLimitInDegPerSec = [];
inputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (-1));

% Set up parameters
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s

Niterations = 5;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
	% Call the function
	[stateTrajectory, t, steeringUsed] = ...
		fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
	% Call the function
	[stateTrajectory, t, steeringUsed] = ...
		fcn_VD_st2dofRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, U, (-1));

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
