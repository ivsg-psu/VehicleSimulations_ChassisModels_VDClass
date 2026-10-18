%% script_test_mdl_VD_st5dofModelFriction_ConstantSpeed.m
% This scripts tests the simulink model 'mdl_VD_st5dofModelFriction_ConstantSpeed'
% It uses bicycle/single-track vehicle model with brush tire model.
% Zero wheel torque, Zero wheel slip, Longitudinal velocity is constant
% Best way would be using a constant velocity controller.

% REVISION HISTORY:
%
% 2021_07_26 by Satya Prasad, szm888@psu.edu
% - First write of function
%
% 2025_12_29 by Sean Brennan, sbrennan@psu.edu
% - Updated header formatting and comments
% - Updated tab stops


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Prep the workspace
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
vehicle.Ca  = [95000; 110000]; % wheel cornering stiffnesses
vehicle.Cx  = [65000; 65000]; % longitudinal stiffnesses

vdParam.contactPatchLength = 0.15;
vdParam.frictionCoefficientRatio = 1;
vdParam.sampling_time_gps = 0.01; % [seconds]
vdParam.sampling_time_imu = 0.01;

%% Define initial conditions
initial.longitudinalSpeed = 25; % longitudinal velocity of vehicle [m/s]
initial.wheelSpeeds = initial.longitudinalSpeed*ones(1,2)/vehicle.Re; % angular velocity of wheel [rad/s]
initial.east    = 0; % initial pose
initial.north   = 0;
initial.heading = 0; % [rad]
road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
friction_coefficient = [0.9, 0.9];

%% Define load transfer conditions
vdParam.longitudinalTransfer = 1;

%% Define inputs to the vehicle model
% Items used to define steering input
steering_amplitude = 2*pi/180; % 2 degrees of steering amplitude for input sinewave
Period = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements

% Define items used to determine how long to run sim
TotalTime = 1.5*Period; % This is how long the simulation will run.
deltaT = 0.01; % simulation step size

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
% Run the simulation in SIMULINK
sim('mdl_VD_st5dofModelFriction_ConstantSpeed.slx', TotalTime);

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
fcn_VD_plotTimeSteeringAngle(time, delta);

fcn_VD_plotTimeSlipAngle(time, alpha);
fcn_VD_plotTimeWheelSlip(time, kappa);
fcn_VD_plotTimeLongitudinalTireForce(time, Fx);
fcn_VD_plotTimeLateralTireForce(time, Fy);
fcn_VD_plotTimeLongitudinalAcceleration(time, States(:,6));
fcn_VD_plotTimeLateralAcceleration(time, States(:,7));

fcn_VD_plotTimeLongitudinalVelocity(time, States(:,1));
fcn_VD_plotTimeLateralVelocity(time, States(:,2));
fcn_VD_plotTimeYawRate(time, States(:,3));
fcn_VD_plotTimeWheelSpeed(time, States(:,(4:5)));
fcn_VD_plotTrajectory(pose(:,[1,2]));
fcn_VD_plotTimeYaw(time, pose(:,3));

fcn_VD_plotTimeAligningMoment(time, Mz);
fcn_VD_plotTimeFriction(time, friction_coefficients_estimate);