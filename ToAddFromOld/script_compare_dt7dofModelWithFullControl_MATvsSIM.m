%% script_compare_dt7dofModelWithFullControl_MATvsSIM.m
% This script tests full control of simulink model with matlab model

% REVISION HISTORY:
%
% 2021_07_20 by Satya Prasad, szm888@psu.edu
% - First write of function
%
% 2025_12_29 by Sean Brennan, sbrennan@psu.edu
% - Updated header formatting and comments
% - Updated tab stops


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)


%% Prepare the workspace
close all % close all the plots

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
controller.velocity_Pgain = 200; % P gain for steering control

load('reference_traversal.mat');
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
initial.longitudinalSpeed = 23; % longitudinal velocity of vehicle [m/s]
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

desired_U = 25; % desired longitudinal velocity [m/s]

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
sim('mdl_VD_dt7dofModelWithFullControl.slx', TotalTime);

%% Run the simulation in MATLAB
% variables to store outputs of Matlab simulation
matlab_time = nan(N_timeSteps,1);
matlab_States = nan(N_timeSteps,9); matlab_pose = nan(N_timeSteps,3);
matlab_alpha = nan(N_timeSteps,4); matlab_kappa = nan(N_timeSteps,4);
matlab_Fz = nan(N_timeSteps,4); matlab_Mz = nan(N_timeSteps,4);
matlab_Fx = nan(N_timeSteps,4); matlab_Fy = nan(N_timeSteps,4);
matlab_friction = nan(N_timeSteps,4);

global flag_update global_acceleration
global_acceleration = zeros(7,1);
input_states = [U;V;r;omega;east;north;heading]; % initial conditions
counter = 1;
for t = 0:deltaT:TotalTime
    matlab_time(counter) = t;
    matlab_States(counter,1:7) = input_states(1:7)';
    matlab_pose(counter,:)     = input_states(8:10)';

    %% Controller: Steering
    pose = matlab_pose(counter,:)';
    target_lookAhead_pose = fcn_VD_snapLookAheadPoseOnToTraversal(pose,...
        reference_traversal,controller);
    steering_angle = fcn_VD_lookAheadLatController(pose,target_lookAhead_pose,...
        controller);
    wheel_torque   = fcn_VD_velocityController(U,desired_U,controller);

    %% Estimate Slips for time 't'
    % Slip Angle/Lateral Slip
    slip_angle = fcn_VD_dtSlipAngle(U,V,r,steering_angle,vehicle);
    % Wheel Slip/Longitudinal Slip
    wheel_slip = fcn_VD_dtWheelSlip(U,V,r,omega,steering_angle,vehicle);
    matlab_alpha(counter,:) = slip_angle';
    matlab_kappa(counter,:) = wheel_slip';

    %% 7-DoF Vehicle Model
    flag_update = true; % set it to to true before every call to RK4 method
    [~,y] = fcn_VD_RungeKutta(@(t,y) fcn_VD_dt7dofModelForController(t,y,...
        steering_angle,wheel_torque,...
        vehicle,road_properties,friction_coefficient',type_of_transfer),...
        input_states,t,deltaT);
    U = y(1); V = y(2); r = y(3); omega = y(4:7);
    input_states = y; clear y;
    matlab_States(counter,8:9) = global_acceleration(1:2)';

    %% Estimate Normal Forces for time 't'
    if 1==counter
        normal_force = fcn_VD_dtNormalForce([0;0],vehicle,road_properties,...
            type_of_transfer);
    else
        normal_force = fcn_VD_dtNormalForce(matlab_States(counter-1,8:9)',vehicle,...
            road_properties,type_of_transfer);
    end
    matlab_Fz(counter,:) = normal_force';

    %% Estimate Tire forces for time 't'
    tire_force = fcn_VD_dtTireForceBrush(slip_angle,wheel_slip,normal_force,...
        friction_coefficient',vehicle);
    matlab_Fx(counter,:) = tire_force(:,1)';
    matlab_Fy(counter,:) = tire_force(:,2)';

    %% Aligning Moment
    aligning_moment = fcn_VD_dtAligningMomentBrush(slip_angle,normal_force,...
        friction_coefficient',vehicle);
    matlab_Mz(counter,:) = aligning_moment';

    %% Friction Estimation
    estimated_vel = matlab_States(counter,1:3)'+...
        [normrnd(0,1.4142*0.25)*vdParam.sampling_time_gps; ...
        normrnd(0,1.4142*0.25)*vdParam.sampling_time_gps; ...
        normrnd(0,1.4142*0.009)*vdParam.sampling_time_gps];
    if 1==counter
        estimated_accel = [0; 0]+...
            [normrnd(0,0.5)*(vdParam.sampling_time_imu^2); ...
            normrnd(0,0.5)*(vdParam.sampling_time_imu^2)];
    else
        estimated_accel = matlab_States(counter-1,8:9)'+...
            [normrnd(0,0.5)*(vdParam.sampling_time_imu^2); ...
            normrnd(0,0.5)*(vdParam.sampling_time_imu^2)];
    end
    slip_angle = fcn_VD_dtSlipAngle(estimated_vel(1),estimated_vel(2),...
        estimated_vel(3),steering_angle,vehicle);
    normal_force = fcn_VD_dtNormalForce(estimated_accel,vehicle,road_properties,...
        type_of_transfer);
    friction_estimate = fcn_VD_estimateFrictionCoefficient(slip_angle,...
        normal_force,aligning_moment,vehicle);
    matlab_friction(counter,:) = friction_estimate';

    counter = counter+1;
end

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