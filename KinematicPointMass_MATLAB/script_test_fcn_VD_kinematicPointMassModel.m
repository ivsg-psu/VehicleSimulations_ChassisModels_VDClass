%% script_test_fcn_VD_kinematicPointMassModel.m
% This scripts test the model 'fcn_VD_kinematicPointMassModel'

% REVISION HISTORY:
%
% 2026_01_26 by Sean Brennan, sbrennan@psu.edu
% - First write of function, using fcn_VD_bicycle2dofModel as starter


% TO-DO:
% - 2026_01_26 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Prep the workspace
close all;  % close all the plots

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
% Set core simulation inputs: the speed, and steering amplitude
U = 20;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: mph ~= 2* m/s)

% Items used to define steering inputs
steering_amplitude_degrees = 2; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim, number of time
% points, etc.
TotalTime = 1.5*Period;  % This is how long the simulation will run. Usually 1.5 times the period is enough.
deltaT = 0.01; % This is the time step of the simulation. See the "Model Settings" submenu in Simulink to see where this variable is used.
N_timeSteps = floor(TotalTime/deltaT) + 1; % This is the number of time steps we should have

%% Fill in the vehicle parameters. 
% Use a structure array so we can have several vehicles

% Approximately a Ford Taurus
clear vehicles
vehicles(1).m          = 1031.9; % kg
vehicles(1).Izz         = 1850; % kg-m^2
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
% V: lateral speed
% r: rotational rate of vehicle around z-axis

% Initialize all the arrays with Not-a-Number (nan). When plotting, any
% values that are nan will be left empty.
all_X   = nan(N_timeSteps,N_vehicles);
all_Y   = nan(N_timeSteps,N_vehicles);
all_phi = nan(N_timeSteps,N_vehicles);
t       = nan(N_timeSteps,N_vehicles);

% Loop through all the vehicles, simulating the trajectory of each within
% the for loop
for ith_vehicle = 1:N_vehicles
    % Print to the console which vehicle we are working on
    fprintf(1,'Working on vehicle: %d\n', ith_vehicle);
    
    %% RK4 in MATLAB Script
    % Set initial conditions
    input_states = [0; 0; 0];
    
    counter = 1; % Initialize the counter
    for time = 0:deltaT:TotalTime-deltaT
        t(counter,ith_vehicle)       = time; % Update time
        inputOmega = steering_amplitude_degrees*sin((2*pi/Period)*time); % steering angle

        all_X(counter,ith_vehicle)   = input_states(1);
        all_Y(counter,ith_vehicle)   = input_states(2);
        all_phi(counter,ith_vehicle) = input_states(3);

		[~, y] = fcn_VD_RungeKutta(...
                @(t,y) fcn_VD_kinematicPointMassModel(y,inputOmega, U), ...
                input_states, time, deltaT);
        input_states = y;
        clear y;
        
        counter = counter+1; % increase the counter
    end
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

index_to_label = (2==t);
% Loop through each of the vehicles
for ith_vehicle=1:N_vehicles
    
    % The XY Plots
    h3 = figure(77);
    hold on;
    set(h3,'Name','XYposition')
    plot(all_X(:,ith_vehicle),all_Y(:,ith_vehicle),'b');
    grid on;
    text(all_X(index_to_label,ith_vehicle), all_Y(index_to_label,ith_vehicle), ...
        ['vehicle' num2str(ith_vehicle)]);
    xlabel('X position [m]');
    ylabel('Y position [m]');
    title('Position');
end
% 
% % The Steering Input
% h4 = figure(66);
% hold on;
% set(h4,'Name','Steering Input')
% plot(t(:,ith_vehicle),df,'b');
% grid on;
% xlabel('Time (sec)'); 
% ylabel('Steering input (rad)');
% title('Steering Input');