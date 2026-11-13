% script_test_fcn_VD_st2dofSlipAngles.m
% tests fcn_VD_st2dofSlipAngles.m

% REVISION HISTORY:
%
% As: script_test_fcn_VD_slipAnglesBicycleModel
%
% 2026_10_22 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_slipAnglesBicycleModel
%   % * First write of script
%   % * Used script_test_fcn_VD_wheelVelocitie+sBicycleModel as starter
%
% As: script_test_fcn_VD_st2dofDerivatives
%
% 2026_11_10 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_st2dofDerivatives
%   % * Renamed function
%   % * Used to be called script_test_fcn_VD_slipAnglesBicycleModel
%
% 2026_11_11 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_st2dofDerivatives
%   % * Updated test cases to match change in formatting for inputs/outputs
%   % * Updated test cases to check new flag input

% TO-DO:
%
% 2026_10_22 by Sean Brennan, sbrennan@psu.edu
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
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));

U = 0;  % U is forward velocity of vehicle CG in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
V = 0; % V is the side velocity of CG, in lateral direction
r = 10*pi/180; % Units are rad/sec
steeringAngles = [10 0]*pi/180; % Steering inputs, units are radians
flagUseLinear = [];

% Call the function
slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(slipAnglesInRadians));

% Check variable sizes
assert(size(slipAnglesInRadians,1)==1);
assert(size(slipAnglesInRadians,2)==2);

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call, linear model
figNum = 10002;
titleString = sprintf('DEMO case: basic call, linear model');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set up inputs
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));

U = 0;  % U is forward velocity of vehicle CG in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
V = 0; % V is the side velocity of CG, in lateral direction
r = 10*pi/180; % Units are rad/sec
steeringAngles = [10 0]*pi/180; % Steering inputs, units are radians
flagUseLinear = [];

% Call the function
slipAnglesInRadiansNonlinear = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (figNum));

% Call the function
flagUseLinear = 1;
slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(slipAnglesInRadians));

% Check variable sizes
assert(size(slipAnglesInRadians,1)==1);
assert(size(slipAnglesInRadians,2)==2);

% Check variable values
assert(~isequal(slipAnglesInRadians,slipAnglesInRadiansNonlinear));

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
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));

U = 10;  % U is forward velocity of vehicle CG in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
V = 0; % V is the side velocity of CG, in lateral direction
r = 300*pi/180; % Units are rad/sec
steeringAngles = [10 0]*pi/180; % Steering inputs, units are radians
flagUseLinear = [];

% Call the function
slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(slipAnglesInRadians));

% Check variable sizes
assert(size(slipAnglesInRadians,1)==1);
assert(size(slipAnglesInRadians,2)==2);

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
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));

U = 10;  % U is forward velocity of vehicle CG in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
V = 0; % V is the side velocity of CG, in lateral direction
r = 300*pi/180; % Units are rad/sec
steeringAngles = [10 0]*pi/180; % Steering inputs, units are radians
flagUseLinear = [];

% Call the function
slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(slipAnglesInRadians));

% Check variable sizes
assert(size(slipAnglesInRadians,1)==1);
assert(size(slipAnglesInRadians,2)==2);

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
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));

U = 10;  % U is forward velocity of vehicle CG in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
V = 0; % V is the side velocity of CG, in lateral direction
r = 300*pi/180; % Units are rad/sec
steeringAngles = [10 0]*pi/180; % Steering inputs, units are radians
flagUseLinear = [];

Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
	% Call the function
	slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
	% Call the function
	slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (-1));
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
