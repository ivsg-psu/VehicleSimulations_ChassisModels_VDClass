%% script_compare_dt7dofModelWithSteeringControl_MATvsSIM.m
% This script tests steering control of simulink model with matlab model

% REVISION HISTORY:
%
% 2021_07_20 by Satya Prasad, szm888@psu.edu
% - First write of function
%
% 2024_03_21 - S. Brennan
% - removed workspace prep and moved this into the main demo code area
%
% 2025_12_29 by Sean Brennan, sbrennan@psu.edu
% - Updated header formatting and comments
% - Updated tab stops


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)


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
%% Define vehicle and controller properties
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

vehicle.contact_patch_length = 0.15; % [meters]
vehicle.friction_ratio = 1; % [No units]

controller.look_ahead_distance = 20; % look-ahead distance [meters]
controller.steering_Pgain = 0.1; % P gain for steering control

% Load the reference_traversal from saved data. This was created out of the
% Path library using fcn_Path_convertPathToTraversalStructure
load('reference_traversal.mat','reference_traversal');
if 1==0
    fig_num = 378337;
    data.traversal{1} = reference_traversal;
    fcn_Path_plotTraversalsXY(data,fig_num);
end

% NOTE: this looks odd - usually station is first
inputTrajectory = [reference_traversal.X,reference_traversal.Y,...
    reference_traversal.Yaw,reference_traversal.Station];

vdParam.fieldsTrajectory.east = 1;
vdParam.fieldsTrajectory.north = 2;
vdParam.fieldsTrajectory.yaw = 3;
vdParam.fieldsTrajectory.station = 4;

vdParam.searchDistance = 4; % [meters]
vdParam.trajectorySize = size(inputTrajectory);
vdParam.sampling_time_gps = 0.01; % [seconds]
vdParam.sampling_time_imu = 0.01;

%% Define initial conditions
% Parameters and initial conditions for simulink model
initial.longitudinalSpeed = 25; % longitudinal velocity of vehicle [m/s]
initial.wheelSpeeds = initial.longitudinalSpeed*ones(1,4)/vehicle.Re; % angular velocity of wheel [rad/s]
initial.east    = 0; % initial pose
initial.north   = 0;
initial.heading = 0; % [rad]

% Parameters and initial conditions for matlab model
U = initial.longitudinalSpeed; % longitudinal velocity [m/s]
V = 0; % lateral velocity [m/s]
r = 0; % yaw rate [rad/s]
omega = initial.wheelSpeeds'; % angular velocity of wheel [rad/s]
east = 0; north = 0; heading = 0; % initial pose

% Parameters and initial conditions for simulink and matlab model
road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
friction_coefficient = [0.9, 0.9, 0.9, 0.9];

%% Define load transfer conditions
vdParam.longitudinalTransfer = 1;
if vdParam.longitudinalTransfer
    vdParam.lateralTransfer = 1;
    type_of_transfer = 'both';
else
    vdParam.lateralTransfer = 0;
    type_of_transfer = 'default';
end

%% Define inputs to the vehicle model
wheel_torque = [0, 0, 5, 5]; % wheel torque [Nm]

% Define items used to determine how long to run sim
TotalTime = 6; % This is how long the simulation will run.
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
sim('mdl_VD_dt7dofModelPathFollowing.slx', TotalTime);

%% Plots to check MATLAB simulation with Simulink
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

fcn_VD_plotCompareLongitudinalAcceleration(time,States(:,8),'Simulink',...
    matlab_time,matlab_States(:,8),'Matlab');
fcn_VD_plotCompareLateralAcceleration(time,States(:,9),'Simulink',...
    matlab_time,matlab_States(:,9),'Matlab');

fcn_VD_plotCompareLongitudinalVelocity(time,States(:,1),'Simulink',...
    matlab_time,matlab_States(:,1),'Matlab');
fcn_VD_plotCompareLateralVelocity(time,States(:,2),'Simulink',...
    matlab_time,matlab_States(:,2),'Matlab');
fcn_VD_plotCompareYawRate(time,States(:,3),'Simulink',...
    matlab_time,matlab_States(:,3),'Matlab');
fcn_VD_plotCompareWheelSpeed(time,States(:,(4:7)),'Simulink',...
    matlab_time,matlab_States(:,(4:7)),'Matlab');

fcn_VD_plotCompareTrajectory([reference_traversal.X, reference_traversal.Y],...
    'Input',simulink_pose(:,[1,2]),'Simulink Output');
fcn_VD_plotCompareTrajectory([reference_traversal.X, reference_traversal.Y],...
    'Input',matlab_pose(:,[1,2]),'Matlab Output');
fcn_VD_plotCompareTrajectory(simulink_pose(:,[1,2]),'Simulink',...
    matlab_pose(:,[1,2]),'Matlab');
fcn_VD_plotCompareYaw(time,simulink_pose(:,3),'Simulink',...
    matlab_time,matlab_pose(:,3),'Matlab');

fcn_VD_plotTimeFriction(matlab_time,matlab_friction);
fcn_VD_plotTimeFriction(time,friction_coefficients_estimate);