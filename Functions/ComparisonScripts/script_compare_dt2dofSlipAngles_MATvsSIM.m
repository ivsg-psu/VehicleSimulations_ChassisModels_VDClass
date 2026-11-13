%% script_compare_st2dofSlipAngles_MATvsSIM.m.m
% This script compares MATLAB function 'fcn_VD_st2dofSlipAngles' with Simulink
% model 'mdl_VD_st2dofSlipAngles_R2025b_SLX.slx'

% REVISION HISTORY:
%
% As: script_compare_dtSlipAngles_MATvsSIM
%
% 2021_07_06 by Satya Prasad, szm888@psu.edu
% - In script_compare_dtSlipAngles_MATvsSIM
%   % * First write of function
%
% 2025_12_29 by Sean Brennan, sbrennan@psu.edu
% - In script_compare_dtSlipAngles_MATvsSIM
%   % * Updated header formatting and comments
%   % * Updated tab stops
%
% As: script_compare_st2dofSlipAngles_MATvsSIM
%
% 2026_11_11 by Sean Brennan, sbrennan@psu.edu
% - In script_compare_dtSlipAngles_MATvsSIM
%   % * Renamed script to be consistent with 2DOF testing
%   % * Updated fcn_VD_dtSlipAngle to new name, fcn_VD_st2dofSlipAngles
%   % * Updated vehicle variable name to vehicleParameters, for consistency

% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Prepare the workspace
close all; % close all the plots

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
initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
deltaT = 0.01; % Units are [sec]
startTime = 0;
endTime = 4.5;
timeInterval = [startTime endTime];  % Units are [sec]

% Set up inputs
inputFrontRoadWheelAngleRadians = 10*pi/180; % radians
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
simulationTimes = (startTime:deltaT:endTime)';
steeringAndTimeInputs = [simulationTimes inputFrontRoadWheelAngleRadians*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]
N_timeSteps = length(simulationTimes(:,1));

% Set up parameters
% USES:
% vehicleParameters.DistancesWheelbase   = 2.6; % wheelbase (m)
% vehicleParameters.DistancesCGFromRearAxle   = 1.3; % length from rear axle to CG (m)
vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters_2013ToyotaPriusFour.mat');
vehicleParameters = fcn_VD_getStudentVehicle((vehicleParametersFile), (-1));
U = 20; % Forward velocity, in m/s
V = 2; % Lateral velocity, in m/s
r = 5*pi/180; % Yaw rate, in rad/sec
flagUseLinear = 1;

% Call the MATLAB function

allSlipsMATLAB = zeros(N_timeSteps,2);
for ith_test = 1:N_timeSteps
	steeringAngles = [steeringAndTimeInputs(ith_test,2) 0];
	slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (-1));
	allSlipsMATLAB(ith_test,:) = slipAnglesInRadians;
end

% Call the Simulink runing RK4 solver
% Which model to use? Depends on the version
currentMATLABversion = version('-release');
modelToUse = cat(2,'mdl_VD_st2dofSlipAngles_R',currentMATLABversion,'_SLX.slx');

% Set the parameters for simulink (must be done since function variables
% are local and not automatically passed into the model)
b = vehicleParameters.DistancesCGFromRearAxle;
a = vehicleParameters.DistancesWheelbase - b;
parameters.a = a;
parameters.b = b;
inputsVsTime = steeringAndTimeInputs;

modelWithParameters = Simulink.SimulationInput(modelToUse); % Copies model from the one given into a "hold" model
modelWithParameters = modelWithParameters.setVariable('U', U);                              % makes U available to model
modelWithParameters = modelWithParameters.setVariable('V', V);                              % makes V available to model
modelWithParameters = modelWithParameters.setVariable('r', r);                              % makes V available to model
modelWithParameters = modelWithParameters.setVariable('initialStates', initialStates);      % makes initialStates available to model
modelWithParameters = modelWithParameters.setVariable('startTime', startTime);              % makes startTime available to model
modelWithParameters = modelWithParameters.setVariable('endTime', endTime);                  % makes endTime available to model
modelWithParameters = modelWithParameters.setVariable('deltaT', deltaT);                    % makes deltaT available to model
modelWithParameters = modelWithParameters.setVariable('inputsVsTime', inputsVsTime);        % makes inputsVsTime available to model
modelWithParameters = modelWithParameters.setVariable('parameters', parameters);            % makes parameters available to model

% Run the simulation in SIMULINK
simout = sim(modelWithParameters);

% Save the results in a big array (for plotting in next part)
% Before saving, we need to check if the full vector is shorter than
% expected length of N_timeSteps
N_simSteps = length(simout.simTime);
if  N_simSteps~= N_timeSteps
	warning('More time was spent (N=%.0f time steps) than expected (expected %.0f time steps) in the simulation. Keeping only the expected time portion.', N_simSteps,N_timeSteps);
end

% Keep the shorter of either the actual length, or expected length:
shorter_index = min(N_timeSteps,N_simSteps);

% Initialize variables
stateTrajectory = nan(shorter_index,3);
t = nan(shorter_index,1);

% Fill in the results to save
t(:,1)                 = simout.simTime(1:shorter_index,1); % Update time
alphaFrontSIMULINK     = simout.alphaFront;
alphaRearSIMULINK      = simout.alphaRear;
allSlipsSIMULINK       = [alphaFrontSIMULINK alphaRearSIMULINK];

% The slip Plots
subplot(2,1,1);
hPlots = fcn_VD_plotTimeSlipAngle(simulationTimes, allSlipsMATLAB, (figNum));
set(hPlots{1},'LineWidth',5,'Color',0.5*[0 0 1],'DisplayName','F, MATLAB');
set(hPlots{2},'LineWidth',5,'Color',0.5*[1 0 0],'DisplayName','R, MATLAB');
hold on;
hPlots = fcn_VD_plotTimeSlipAngle(simulationTimes, allSlipsSIMULINK, (figNum));
set(hPlots{1},'LineWidth',1,'Color',1*[0 1 1],'DisplayName','F, SIMULINK');
set(hPlots{2},'LineWidth',1,'Color',1*[1 1 0],'DisplayName','R, SIMULINK');

subplot(2,1,2);
plot(simulationTimes, abs(allSlipsMATLAB(:,1) - allSlipsSIMULINK(:,1)),'DisplayName','Front slip');
hold on;
plot(simulationTimes, abs(allSlipsMATLAB(:,2) - allSlipsSIMULINK(:,2)),'DisplayName','Rear slip');
xlabel('Time [s]','Interpreter','Latex','Fontsize',18)
ylabel('Simulation errors [rad]','Interpreter','Latex','Fontsize',18)
grid on

% % Check variable types
% assert(isnumeric(stateTrajectory));
% assert(isnumeric(t));
% assert(isnumeric(steeringUsed));
% 
% % Check variable sizes
% assert(size(stateTrajectory,1)>=1); 
% assert(size(stateTrajectory,2)==3); 
% assert(size(t,1)==size(stateTrajectory,1)); 
% assert(size(t,2)==1); 
% assert(size(steeringUsed,1)==size(stateTrajectory,1)); 
% assert(size(steeringUsed,2)==1); 

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

% %% Basic example - NO FIGURE
% figNum = 80001;
% fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
% figure(figNum); close(figNum);
% 
% % Set the simulation time/state arguments
% initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
% deltaT = 0.01; % Units are [sec]
% startTime = 0;
% endTime = 4.5;
% timeInterval = [startTime endTime];  % Units are [sec]
% 
% % Set up inputs
% steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
% Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
% simulationTimes = (startTime:deltaT:endTime)';
% steeringAndTimeInputs = [simulationTimes steering_amplitude_degrees*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]
% 
% % Set up parameters
% clear parameters
% parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
% parameters.L = 2.5; % wheelbase in meters
%
% % Call the function
% [stateTrajectory, t, steeringUsed] = ...
% fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
% timeInterval, steeringAndTimeInputs, U, ([]));
% 
% % sgtitle(titleString, 'Interpreter','none');
% 
% % Check variable types
% assert(isnumeric(stateTrajectory));
% assert(isnumeric(t));
% assert(isnumeric(steeringUsed));
% 
% % Check variable sizes
% assert(size(stateTrajectory,1)>=1); 
% assert(size(stateTrajectory,2)==3); 
% assert(size(t,1)==size(stateTrajectory,1)); 
% assert(size(t,2)==1); 
% assert(size(steeringUsed,1)==size(stateTrajectory,1)); 
% assert(size(steeringUsed,2)==1); 
% 
% % Check variable values
% % (too complex to check)
% 
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==figNum));
% 
% 
% %% Basic fast mode - NO FIGURE, FAST MODE
% figNum = 80002;
% fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
% figure(figNum); close(figNum);
% 
% % Set the simulation time/state arguments
% initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
% deltaT = 0.01; % Units are [sec]
% startTime = 0;
% endTime = 4.5;
% timeInterval = [startTime endTime];  % Units are [sec]
% 
% % Set up inputs
% steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
% Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
% simulationTimes = (startTime:deltaT:endTime)';
% steeringAndTimeInputs = [simulationTimes steering_amplitude_degrees*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]
% 
% % Set up parameters
% clear parameters
% parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
% parameters.L = 2.5; % wheelbase in meters
% 
% % Call the function
% [stateTrajectory, t, steeringUsed] = ...
% fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
% timeInterval, steeringAndTimeInputs, U, (-1));
% 
% % sgtitle(titleString, 'Interpreter','none');
% 
% % Check variable types
% assert(isnumeric(stateTrajectory));
% assert(isnumeric(t));
% assert(isnumeric(steeringUsed));
% 
% % Check variable sizes
% assert(size(stateTrajectory,1)>=1); 
% assert(size(stateTrajectory,2)==3); 
% assert(size(t,1)==size(stateTrajectory,1)); 
% assert(size(t,2)==1); 
% assert(size(steeringUsed,1)==size(stateTrajectory,1)); 
% assert(size(steeringUsed,2)==1); 
% 
% % Check variable values
% % (too complex to check)
% 
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==figNum));
% 
% 
% %% Compare speeds of pre-calculation versus post-calculation versus a fast variant
% figNum = 80003;
% fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
% figure(figNum);
% close(figNum);
% 
% % Set the simulation time/state arguments
% initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
% deltaT = 0.01; % Units are [sec]
% startTime = 0;
% endTime = 4.5;
% timeInterval = [startTime endTime];  % Units are [sec]
% 
% % Set up inputs
% steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
% Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
% simulationTimes = (startTime:deltaT:endTime)';
% steeringAndTimeInputs = [simulationTimes steering_amplitude_degrees*sin((2*pi/Period)*simulationTimes)]; % [times steering angles]
% 
% % Set up parameters
% clear parameters
% parameters.U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 1 mph ~= 2* m/s)
% parameters.L = 2.5; % wheelbase in meters
% 
% Niterations = 50;
% 
% % Do calculation without pre-calculation
% tic;
% for ith_test = 1:Niterations
% 	% Call the function
% 	[stateTrajectory, t, steeringUsed] = ...
% 		fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
% 		timeInterval, steeringAndTimeInputs, U, ([]));
% end
% slow_method = toc;
% 
% % Do calculation with pre-calculation, FAST_MODE on
% tic;
% for ith_test = 1:Niterations
% 	% Call the function
% 	[stateTrajectory, t, steeringUsed] = ...
% 		fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
% 		timeInterval, steeringAndTimeInputs, U, (-1));
% end
% fast_method = toc;
% 
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==figNum));
% 
% % Plot results as bar chart
% figure(373737);
% clf;
% hold on;
% 
% X = categorical({'Normal mode','Fast mode'});
% X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
% Y = [slow_method fast_method ]*1000/Niterations;
% bar(X,Y)
% ylabel('Execution time (Milliseconds)')
% 
% 
% % Make sure plot did NOT open up
% figHandles = get(groot, 'Children');
% assert(~any(figHandles==figNum));


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


