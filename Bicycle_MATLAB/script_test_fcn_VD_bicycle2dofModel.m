%% script_test_fcn_VD_bicycle2dofModel.m
% This scripts test the model 'fcn_VD_bicycle2dofModel'
%
% Author: Satya Prasad on 2021/05/16
% Questions or comments? szm888@psu.edu

%% Prep the workspace
close all;  % close all the plots
clear
clc

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
vehicle(1).m          = 1031.9; % kg
vehicle(1).Iz         = 1850; % kg-m^2
vehicle(1).a          = 0.9271; % Distance from front axle to CG, in meters
vehicle(1).b          = 1.5621;  % Distance from rear axle to CG, in meters
vehicle(1).Caf        = 77500; % N/rad;
vehicle(1).Car	      = 116250; % N/rad;

% The number of vehicles is the length of the structure array
N_vehicles = length(vehicle);

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
all_r   = nan(N_timeSteps,N_vehicles);
all_V   = nan(N_timeSteps,N_vehicles);
t       = nan(N_timeSteps,N_vehicles);
df      = nan(N_timeSteps,N_vehicles);

% Loop through all the vehicles, simulating the trajectory of each within
% the for loop
for vehicle_i=1:N_vehicles
    % Print to the console which vehicle we are working on
    fprintf(1,'Working on vehicle: %d\n', vehicle_i);
    
    %% RK4 in MATLAB Script
    % Fill in the parameters needed by the simulation for it to run
    m = vehicle(vehicle_i).m;
    Iz = vehicle(vehicle_i).Iz;
    a = vehicle(vehicle_i).a;
    b = vehicle(vehicle_i).b;
    L = a+b;
    Caf = vehicle(vehicle_i).Caf;
    Car = vehicle(vehicle_i).Car;
    
    % Set initial conditions
    input_states = [0; 0; 0; 0; 0];
    
    counter = 1; % Initialize the counter
    for time = 0:deltaT:TotalTime-deltaT
        t(counter,vehicle_i)       = time; % Update time
        delta_f = (1-1*(0<time-Period))*...
                  (pi/180)*steering_amplitude_degrees*sin((2*pi/Period)*time); % steering angle
        df(counter,vehicle_i)      = delta_f;
        all_V(counter,vehicle_i)   = input_states(1);
        all_r(counter,vehicle_i)   = input_states(2);
        all_X(counter,vehicle_i)   = input_states(3);
        all_Y(counter,vehicle_i)   = input_states(4);
        all_phi(counter,vehicle_i) = input_states(5);
        
        [~, y] = fcn_VD_RungeKutta(...
                @(t,y) fcn_VD_bicycle2dofModel(t,y,steering_amplitude_degrees,Period,U,vehicle), ...
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

ind = (2==t);
% Loop through each of the vehicles
for vehicle_i=1:N_vehicles
    % Plot the yaw rate
    h1 = figure(99);
    hold on;
    set(h1,'Name','Yawrate')
    plot(t(:,vehicle_i),all_r(:,vehicle_i),'b');
    grid on;
    text(2, all_r(ind,vehicle_i), ['vehicle' num2str(vehicle_i)]);
    xlabel('Time (sec)');
    ylabel('Yawrate (rad/sec)');
    title('Yawrate');
    
    % Plot the lateral velocity
    h2 = figure(88);
    hold on;
    set(h2,'Name','LatVel')
    plot(t(:,vehicle_i),all_V(:,vehicle_i),'b');
    grid on;
    text(2, all_V(ind,vehicle_i), ['vehicle' num2str(vehicle_i)]);
    xlabel('Time (sec)');
    ylabel('Lateral Velocity (m/sec)');
    title('Lateral Velocity');
    
    % The XY Plots
    h3 = figure(77);
    hold on;
    set(h2,'Name','XYposition')
    plot(all_X(:,vehicle_i),all_Y(:,vehicle_i),'b');
    grid on;
    text(all_X(ind,vehicle_i), all_Y(ind,vehicle_i), ...
        ['vehicle' num2str(vehicle_i)]);
    xlabel('X position [m]');
    ylabel('Y position [m]');
    title('Position');
end

% The Steering Input
h4 = figure(66);
hold on;
set(h4,'Name','Steering Input')
plot(t(:,vehicle_i),df,'b');
grid on;
xlabel('Time (sec)'); 
ylabel('Steering input (rad)');
title('Steering Input');