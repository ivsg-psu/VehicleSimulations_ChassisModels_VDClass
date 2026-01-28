%% script_compare_kinematicPointMass_MATvsSIM.m.m
% This script compares MATLAB with Simulink
% for the kinematicPointMass

% REVISION HISTORY:
%
% 2026_01_27 by Sean Brennan, sbrennan@psu.edu
% - First write of function


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Prepare the workspace
close all; % close all the plots

figNum = 777;
figure(figNum); clf;


%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set core simulation inputs: time, the speed, and steering amplitude

% Items used to define steering inputs
steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim, number of time
% points, etc.
TotalTime = 1.5*Period;  % This is how long the simulation will run. Usually 1.5 times the period is enough.
deltaT = 0.01; % This is the time step of the simulation. See the "Model Settings" submenu in Simulink to see where this variable is used.
simulationTimes = (0:deltaT:TotalTime)';
N_timeSteps = length(simulationTimes); % This is the number of time steps we should have

U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: mph ~= 2* m/s)

steeringSIM = steering_amplitude_degrees*sin((2*pi/Period)*(simulationTimes)); % steering angle
steeringMAT = steering_amplitude_degrees*sin((2*pi/Period)*(simulationTimes)); % steering angle

steeringInputs = [simulationTimes steeringSIM];

%% Fill in the vehicle parameters. 
% Use a structure array so we can have several vehicles

% Approximately a Ford Taurus
clear vehicles
vehicles(1).m          = 1031.9; % kg
vehicles(1).Izz        = 1850; % kg-m^2
vehicles(1).a          = 0.9271; % Distance from front axle to CG, in meters
vehicles(1).b          = 1.5621;  % Distance from rear axle to CG, in meters
vehicles(1).Caf        = 77500; % N/rad;
vehicles(1).Car	       = 116250; % N/rad;

% The number of vehicles is the length of the structure array
N_vehicles = length(vehicles);

%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This runs the codes to predict core behavior within SIMULINK and MATLAB,
% looping through each vehicle.
% For this simulation set, the core behavior is
% X: the global X position
% Y: the global Y position
% theta: the global yaw angle

%% Run the simulation in SIMULINK
vehicle_i = 1;

all_X   = nan(N_timeSteps,N_vehicles);
all_Y   = nan(N_timeSteps,N_vehicles);
% all_phi = nan(N_timeSteps,N_vehicles);


mdl = 'mdl_VD_KinematicPointMassModel.slx';           % model name (open or on path)
load_system(mdl);          % optional: load without opening editor

% Run with stop time 10
simOut = sim(mdl, 'StopTime', 'TotalTime');

% Save the results in a big array (for plotting in next part)
% Before saving, we need to check if the full vector is shorter than
% expected length of N_timeSteps
if length(simOut.t) ~= N_timeSteps
	warning('More time was spent than expected in the simulation. Keeping only the expected time portion.')
end

% Keep the shorter of either the actual length, or expected length:
shorter_index = min(N_timeSteps,length(simOut.t));

% Fill in the data arrays
all_X(1:shorter_index,vehicle_i)   = simOut.X(1:shorter_index);
all_Y(1:shorter_index,vehicle_i)   = simOut.Y(1:shorter_index);
% all_phi(1:shorter_index,vehicle_i) = simOut.phi(1:shorter_index);

% Plot results
indexToAttachLabel = (2==t);
% Loop through each of the vehicles
for ith_vehicle=1:N_vehicles
    trajectory = [all_X(:,ith_vehicle),all_Y(:,ith_vehicle)];

    % The XY Plots
	subplot(3,1,1);
    h_plot = fcn_VD_plotTrajectory(trajectory,(figNum));
    text(all_X(indexToAttachLabel,ith_vehicle), all_Y(indexToAttachLabel,ith_vehicle), ...
        ['vehicle' num2str(ith_vehicle)]);
end
set(h_plot,'DisplayName','Simulink XY','LineWidth',5);


%% Run the simulation in MATLAB
% Initialize all the arrays with Not-a-Number (nan). When plotting, any
% values that are nan will be left empty.
all_Xmat   = nan(N_timeSteps,N_vehicles);
all_Ymat   = nan(N_timeSteps,N_vehicles);
all_phi = nan(N_timeSteps,N_vehicles);
t       = nan(N_timeSteps,N_vehicles);

% Loop through all the vehicles, simulating the trajectory of each within
% the for loop
ith_vehicle = 1;
% Print to the console which vehicle we are working on
fprintf(1,'Working on vehicle: %d\n', ith_vehicle);

% RK4 in MATLAB Script
% Set initial conditions
initialStates = [0; 0; 0];
currentStates = initialStates;

for ith_time = 1:N_timeSteps
	thisTime = simulationTimes(ith_time);
	t(ith_time,ith_vehicle)       = simulationTimes(ith_time); % Update time
	inputOmega = steeringMAT(ith_time,1);

	% Fill in the results to save
	all_Xmat(ith_time,ith_vehicle)   = currentStates(1);
	all_Ymat(ith_time,ith_vehicle)   = currentStates(2);
	all_phi(ith_time,ith_vehicle)    = currentStates(3);

	% Use Runga-Kutta to predict next position
	[~, y] = fcn_VD_RungeKutta(...
		@(t,y) fcn_VD_kinematicPointMassModel(y,inputOmega, U), ...
		currentStates, thisTime, deltaT);

	currentStates = y;
end

%% Plot the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____  _       _   _   _             
%  |  __ \| |     | | | | (_)            
%  | |__) | | ___ | |_| |_ _ _ __   __ _ 
%  |  ___/| |/ _ \| __| __| | '_ \ / _` |
%  | |    | | (_) | |_| |_| | | | | (_| |
%  |_|    |_|\___/ \__|\__|_|_| |_|\__, |
%                                   __/ |
%                                  |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the plotting. Both the plots for each situation should agree with 
% each other (otherwise the equations are not consistent)

indexToAttachLabel = (2==t);
% Loop through each of the vehicles
trajectory = [all_Xmat(:,ith_vehicle),all_Ymat(:,ith_vehicle)];

% The XY Plots
subplot(3,1,1);
h_plot = fcn_VD_plotTrajectory(trajectory,(figNum));
text(all_Xmat(indexToAttachLabel,ith_vehicle), all_Ymat(indexToAttachLabel,ith_vehicle), ...
	['vehicle' num2str(ith_vehicle)]);
set(h_plot,'DisplayName','MATLAB XY','LineWidth',5);

subplot(3,1,2);
plot(simulationTimes, abs(all_Xmat - all_X));
xlabel('Time [s]','Interpreter','Latex','Fontsize',18)
ylabel('X-error [m]','Interpreter','Latex','Fontsize',18)
grid on

subplot(3,1,3);
plot(simulationTimes, abs(all_Ymat - (all_Y+100*eps)));
xlabel('Time [s]','Interpreter','Latex','Fontsize',18)
ylabel('Y-error [m]','Interpreter','Latex','Fontsize',18)
grid on

% % The Steering Input
% h4 = figure(66);
% hold on;
% set(h4,'Name','Steering Input')
% plot(t(:,ith_vehicle),df,'b');
% grid on;
% xlabel('Time (sec)'); 
% ylabel('Steering input (rad)');
% title('Steering Input');


