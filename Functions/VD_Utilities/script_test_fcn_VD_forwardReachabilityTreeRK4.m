% script_test_fcn_VD_forwardReachabilityTreeRK4.m
% tests fcn_VD_forwardReachabilityTreeRK4.m

% REVISION HISTORY:
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - In script_test_fcn_VD_forwardReachabilityTreeRK4
%   % * Wrote the code originally, 
%   % * Using script_test_fcn_VD_kinematicPointMassModelRK4 as starter


% TO-DO:
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
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
titleString = sprintf('DEMO case: basic call with kinematic point mass model');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)

modelIDToUse = 0; % Kinematic point mass model

% Call the function
[stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectories));
assert(isnumeric(times));
assert(isnumeric(steeringAnglesUsed));

% Check variable sizes
assert(size(stateTrajectories,1)>=1); 
assert(size(stateTrajectories,2)==3); 
assert(size(times,1)==size(stateTrajectories,1)); 
assert(size(times,2)==1); 
assert(size(steeringAnglesUsed,1)==size(stateTrajectories,1)); 
assert(size(steeringAnglesUsed,2)==1); 

% Check variable values
% (too complex to check)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: basic call
figNum = 10002;
titleString = sprintf('DEMO case: basic call with kinematic bicycle model');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 1;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
vehicleParameters.L = 2.5; % wheelbase in meters

modelIDToUse = 1; % Kinematic bicycle model

% Call the function
[stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectories));
assert(isnumeric(times));
assert(isnumeric(steeringAnglesUsed));

% Check variable sizes
assert(size(stateTrajectories,1)>=1); 
assert(size(stateTrajectories,2)==3); 
assert(size(times,1)==size(stateTrajectories,1)); 
assert(size(times,2)==1); 
assert(size(steeringAnglesUsed,1)==size(stateTrajectories,1)); 
assert(size(steeringAnglesUsed,2)==1); 

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

% Set the simulation time/state arguments
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)

modelIDToUse = 0; % Kinematic point mass model

% Call the function
[stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectories));
assert(isnumeric(times));
assert(isnumeric(steeringAnglesUsed));

% Check variable sizes
assert(size(stateTrajectories,1)>=1); 
assert(size(stateTrajectories,2)==3); 
assert(size(times,1)==size(stateTrajectories,1)); 
assert(size(times,2)==1); 
assert(size(steeringAnglesUsed,1)==size(stateTrajectories,1)); 
assert(size(steeringAnglesUsed,2)==1); 

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
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)

modelIDToUse = 0; % Kinematic point mass model

% Call the function
[stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(stateTrajectories));
assert(isnumeric(times));
assert(isnumeric(steeringAnglesUsed));

% Check variable sizes
assert(size(stateTrajectories,1)>=1); 
assert(size(stateTrajectories,2)==3); 
assert(size(times,1)==size(stateTrajectories,1)); 
assert(size(times,2)==1); 
assert(size(steeringAnglesUsed,1)==size(stateTrajectories,1)); 
assert(size(steeringAnglesUsed,2)==1); 

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
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]
steeringInterval = ((0:2:25)')*pi/180;  % Units are [rad]

% Set up parameters
clear vehicleParameters
vehicleParameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)

modelIDToUse = 0; % Kinematic point mass model


Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [stateTrajectories, times, steeringAnglesUsed] = ...
        fcn_VD_forwardReachabilityTreeRK4(...
        initialStates, ...
        deltaT, timeInterval, steeringInterval, vehicleParameters, ...
        modelIDToUse, ([]));


end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [stateTrajectories, times, steeringAnglesUsed] = ...
        fcn_VD_forwardReachabilityTreeRK4(...
        initialStates, ...
        deltaT, timeInterval, steeringInterval, vehicleParameters, ...
        modelIDToUse, (-1));


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
