%% script_compare_st5dofModel_MATvsSIM.m
% This script compares MATLAB function 'fcn_VD_st5dofModel' with Simulink 
% model 'mdl_VD_st5dofModel'
%
% Author: Satya Prasad on 2021/08/13
% Questions or comments? szm888@psu.edu

%% Prepare the workspace
close all; % close all the plots
clear all %#ok<CLALL>
clc

%% Add path
addpath('../Bicycle_Simulink')
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

vehicle.contact_patch_length = 0.15;
vehicle.friction_ratio = 1; % [No units]

%% Define initial conditions
initial.longitudinalSpeed = 25; % longitudinal velocity of vehicle [m/s]
initial.wheelSpeeds = initial.longitudinalSpeed*ones(1,2)/vehicle.Re; % angular velocity of wheel [rad/s]
initial.east    = 0; % initial pose
initial.north   = 0;
initial.heading = 0; % [rad]
% Parameters and initial conditions for matlab model
U = initial.longitudinalSpeed; % longitudinal velocity [m/s]
V = 0; % lateral velocity [m/s]
r = 0; % yaw rate [rad/s]
omega = initial.wheelSpeeds'; % angular velocity of wheel [rad/s]
east = 0; north = 0; heading = 0; % initial pose
road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
friction_coefficient = [0.9, 0.9];

%% Define load transfer conditions
vdParam.longitudinalTransfer = 0;
if vdParam.longitudinalTransfer
    type_of_transfer = 'longitudinal';
else
    type_of_transfer = 'default';
end

%% Define inputs to the vehicle model
% Items used to define steering input
steering_amplitude = 2*pi/180; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
wheel_torque = [0, 5]; % wheel torque [Nm]

% Define items used to determine how long to run sim
TotalTime = 1.5*Period; % This is how long the simulation will run.
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
sim('mdl_VD_st5dofModel.slx', TotalTime);

%% Run the simulation in MATLAB
% variables to store outputs of Matlab simulation
matlab_time = nan(N_timeSteps,1);
matlab_States = nan(N_timeSteps,7); matlab_pose = nan(N_timeSteps,3);
matlab_alpha = nan(N_timeSteps,2); matlab_kappa = nan(N_timeSteps,2);
matlab_Fz = nan(N_timeSteps,2);  matlab_Mz = nan(N_timeSteps,2);
matlab_Fx = nan(N_timeSteps,2); matlab_Fy = nan(N_timeSteps,2);

global flag_update global_acceleration
global_acceleration = zeros(5,1);
input_states = [U;V;r;omega;east;north;heading]; % initial conditions
counter = 1;
for t = 0:deltaT:TotalTime
matlab_time(counter) = t;
matlab_States(counter,1:5) = input_states(1:5)';
matlab_pose(counter,:)     = input_states(6:8)';

%% Estimate Steering Input for time 't'
delta_f = (1-1*(0<t-Period))*steering_amplitude*sin((2*pi/Period)*t); % front steering angle [rad]
steering_angle = [delta_f; 0];

%% Estimate Slips for time 't'
% Slip Angle/Lateral Slip
slip_angle = fcn_VD_stSlipAngle(U,V,r,steering_angle,vehicle);
% Wheel Slip/Longitudinal Slip
wheel_slip = fcn_VD_stWheelSlip(U,V,r,omega,steering_angle,vehicle);
matlab_alpha(counter,:) = slip_angle';
matlab_kappa(counter,:) = wheel_slip';

%% Call to 7-DoF Vehicle Model
flag_update = true; % set it to to true before every call to RK4 method
[~,y] = fcn_VD_RungeKutta(@(t,y) fcn_VD_st5dofModel(t,y,...
    steering_amplitude,Period,wheel_torque',...
    vehicle,road_properties,friction_coefficient',type_of_transfer),...
    input_states,t,deltaT);
U = y(1); V = y(2); r = y(3); omega = y(4:5);
input_states = y; clear y;
matlab_States(counter,6:7) = global_acceleration(1:2)';

%% Estimate Normal Forces for time 't'
if 1==counter
    normal_force = fcn_VD_stNormalForce([0;0],vehicle,road_properties,...
        type_of_transfer);
else
    normal_force = fcn_VD_stNormalForce(matlab_States(counter-1,6:7)',vehicle,...
        road_properties,type_of_transfer);
end
matlab_Fz(counter,:) = normal_force';

%% Estimate Tire forces for time 't'
tire_force = fcn_VD_stTireForceBrush(slip_angle,wheel_slip,normal_force,...
    friction_coefficient',vehicle);
matlab_Fx(counter,:) = tire_force(:,1)';
matlab_Fy(counter,:) = tire_force(:,2)';

%% Aligning Moment
aligning_moment = fcn_VD_stAligningMomentBrush(slip_angle,normal_force,...
                    friction_coefficient',vehicle);
matlab_Mz(counter,:) = aligning_moment';

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
fcn_VD_plotCompareWheelSlip(time,kappa,'Simulink',...
    matlab_time,matlab_kappa,'Matlab');

fcn_VD_plotCompareNormalForce(time,Fz,'Simulink',...
    matlab_time,matlab_Fz,'Matlab');
fcn_VD_plotCompareLongitudinalTireForce(time,Fx,'Simulink',...
    matlab_time,matlab_Fx,'Matlab');
fcn_VD_plotCompareLateralTireForce(time,Fy,'Simulink',...
    matlab_time,matlab_Fy,'Matlab');
fcn_VD_plotCompareAligningMoment(matlab_time,matlab_Mz,'Simulink',...
    matlab_time,matlab_Mz,'Matlab');

fcn_VD_plotCompareLongitudinalAcceleration(time,States(:,6),'Simulink',...
    matlab_time,matlab_States(:,6),'Matlab');
fcn_VD_plotCompareLateralAcceleration(time,States(:,7),'Simulink',...
    matlab_time,matlab_States(:,7),'Matlab');

fcn_VD_plotCompareLongitudinalVelocity(time,States(:,1),'Simulink',...
    matlab_time,matlab_States(:,1),'Matlab');
fcn_VD_plotCompareLateralVelocity(time,States(:,2),'Simulink',...
    matlab_time,matlab_States(:,2),'Matlab');
fcn_VD_plotCompareYawRate(time,States(:,3),'Simulink',...
    matlab_time,matlab_States(:,3),'Matlab');
fcn_VD_plotCompareWheelSpeed(time,States(:,(4:5)),'Simulink',...
    matlab_time,matlab_States(:,(4:5)),'Matlab');

fcn_VD_plotCompareTrajectory(pose(:,[1,2]),'Simulink',...
    matlab_pose(:,[1,2]),'Matlab');
fcn_VD_plotCompareYaw(time,pose(:,3),'Simulink',...
    matlab_time,matlab_pose(:,3),'Matlab');
