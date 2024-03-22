%% script_test_fcn_VD_stNormalForce.m
% This script tests MATLAB function 'fcn_VD_stNormalForce'
%
% Author: Satya Prasad on 2021/08/13
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
close all; % close all the plots
clear all %#ok<CLALL>
clc

%% Add path
addpath('../VD_Utilities')
addpath('../VD_Utilities/Bicycle')

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
%% Define vehicle properties
% Define a MATLAB structure that specifies the physical values for a vehicle.
% For convenience, we ask that you call this stucture 'vehicle'.
vehicle.m   = 1600; % mass (kg)
vehicle.Izz = 2500; % mass moment of inertia (kg m^2)
vehicle.Iw  = 1.2; % mass moment of inertia of a wheel (kg m^2)
vehicle.Re  = 0.32; % effective radius of a wheel (m)
vehicle.a   = 1.3; % length from front axle to CG (m)
vehicle.L   = 2.6; % wheelbase (m)
vehicle.b   = 1.3; % length from rear axle to CG (m)
vehicle.d   = 1.5; % track width (m)
vehicle.h_cg = 0.42; % height of the cg (m)
vehicle.Ca  = 2*[95000; 110000]; % wheel cornering stiffnesses
vehicle.Cx  = 2*[65000; 65000]; % longitudinal stiffnesses

%% Define load transfer conditions
type_of_transfer = 'longitudinal';
% type_of_transfer = 'default';

%% Define inputs to the vehicle model
road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
acceleration = [2, 0.5]; % acceleration
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim
TotalTime = 5; % This is how long the simulation will run.
deltaT = 0.01;
N_timeSteps = floor(TotalTime/deltaT)+1; % This is the number of time steps we should have

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
%% Run the simulation in MATLAB
% variables to store outputs of Matlab simulation
matlab_Fz = nan(N_timeSteps,2);
matlab_time = nan(N_timeSteps,1);

counter = 1;
for t = 0:deltaT:TotalTime
matlab_time(counter) = t;

%% Inputs
ax = acceleration(1)*sin(pi/Period*t);
ay = acceleration(2)*sin(pi/Period*t);

%% Normal Forces
normal_force = fcn_VD_stNormalForce([ax;ay],vehicle,road_properties,...
                type_of_transfer);
matlab_Fz(counter,:) = normal_force';

counter = counter+1;
end

%% Plot results
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
fcn_VD_plotTimeNormalForce(matlab_time,matlab_Fz);
