%% script_compare_dtSlipAngle_MATvsSIM.m
% This script compares MATLAB function 'fcn_VD_dtSlipAngle' with Simulink
% model 'mdl_VD_dtSlipAngle.slx'

% REVISION HISTORY:
%
% 2021_07_06 by Satya Prasad, szm888@psu.edu
% - First write of function
%
% 2025_12_29 by Sean Brennan, sbrennan@psu.edu
% - Updated header formatting and comments
% - Updated tab stops


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Prepare the workspace
close all; % close all the plots

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
vehicle.Ca  = [95000; 95000; 110000; 110000]; % wheel cornering stiffnesses
vehicle.Cx  = [65000; 65000; 65000; 65000]; % longitudinal stiffnesses

%% Define inputs to the vehicle model
U = 24.59; % longitudinal velocity [m/s]
V = 0.0; % lateral velocity [m/s]
r = 0.5; % yaw rate [rad/s]
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
%% Run the simulation in SIMULINK
sim('mdl_VD_dtSlipAngle.slx', TotalTime);

%% Run the simulation in MATLAB
% variables to store outputs of Matlab simulation
matlab_alpha = nan(N_timeSteps,4);
matlab_time  = nan(N_timeSteps,1);

counter = 1;
for t = 0:deltaT:TotalTime
    matlab_time(counter) = t;

    %% Inputs
    delta_f = (1-1*(0<t-Period))*steering_amplitude*sin((2*pi/Period)*t); % front steering angle
    steering_angle = [delta_f; delta_f; 0; 0];

    Vx = abs(U*sin((0.5*pi/Period)*t + pi/2));
    Vy = V*sin((pi/Period)*t);
    yaw_rate = r*sin((2*pi/Period)*t);

    %% Slips
    % Slip Angle/Lateral Slip
    slip_angle = fcn_VD_dtSlipAngle(Vx,Vy,yaw_rate,steering_angle,vehicle);
    matlab_alpha(counter,:) = slip_angle';

    counter = counter+1;
end

%% Plots to compare MATLAB simulation with Simulink simulation
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
fcn_VD_plotCompareSlipAngle(time,alpha,'Simulink',...
    matlab_time,matlab_alpha,'Matlab');
