%% script_testfcn_VD_stWheelSlip.m
% This script tests MATLAB function 'fcn_VD_stWheelSlip'
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

%% Define inputs to the vehicle model
U = 24.59; % longitudinal velocity [m/s]
V = 0.0; % lateral velocity [m/s]
r = 0.5; % yaw rate [rad/s]
omega = 0.98*U/vehicle.Re*ones(2,1); % angular velocity of wheel [rad/s]
steering_amplitude = 2*pi/180; % front steering angle [rad]
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
matlab_kappa = nan(N_timeSteps,2);
matlab_time = nan(N_timeSteps,1);

counter = 1;
for t = 0:deltaT:TotalTime
matlab_time(counter) = t;

%% Inputs
delta_f = (1-1*(0<t-Period))*steering_amplitude*sin((2*pi/Period)*t); % front steering angle
steering_angle = [delta_f; 0];

Vx = abs(U*sin((0.5*pi/Period)*t + pi/2));
Vy = V*sin((pi/Period)*t);
yaw_rate = r*sin((2*pi/Period)*t);

%% Slips
% Wheel Slip/Longitudinal Slip
wheel_slip = fcn_VD_stWheelSlip(Vx,Vy,yaw_rate,omega,steering_angle,vehicle);
matlab_kappa(counter,:) = wheel_slip';

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
fcn_VD_plotTimeWheelSlip(matlab_time,matlab_kappa);
