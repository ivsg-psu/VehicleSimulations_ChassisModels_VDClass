%% script_compare_dtTireForceBrush_MATvsSIM.m
% This script compares MATLAB function 'fcn_VD_dtTireForceBrush' with
% Simulink model 'mdl_VD_dtTireForceBrush.slx'

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
U = 24.59; % longitudinal velocity [m/s]
V = 0; % lateral velocity [m/s]
r = 0.5; % yaw rate [rad/s]
omega = 0.98*U/vehicle.Re*ones(1,4); % angular velocity of wheel [rad/s]
steering_amplitude = 2*pi/180; % front steering angle [rad]
road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
acceleration = [2, 0.5]; % acceleration
friction_coefficient = [0.9, 0.9, 0.9, 0.9];
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
sim('mdl_VD_dtTireForceBrush.slx', TotalTime);

%% Run the simulation in MATLAB
% variables to store outputs of Matlab simulation
matlab_Fx = nan(N_timeSteps,4); matlab_Fy = nan(N_timeSteps,4);
matlab_time = nan(N_timeSteps,1);

counter = 1;
for t = 0:deltaT:TotalTime
    matlab_time(counter) = t;

    %% Inputs
    delta_f = (1-1*(0<t-Period))*steering_amplitude*sin((2*pi/Period)*t); % front steering angle
    steering_angle = [delta_f; delta_f; 0; 0];

    Vx = abs(U*sin((0.5*pi/Period)*t + pi/2));
    Vy = V*sin((pi/Period)*t);
    yaw_rate = r*sin((2*pi/Period)*t);

    ax = acceleration(1)*sin(pi/Period*t);
    ay = acceleration(2)*sin(pi/Period*t);

    %% Slips
    % Slip Angle/Lateral Slip
    slip_angle = fcn_VD_dtSlipAngle(Vx,Vy,yaw_rate,steering_angle,vehicle);
    % Wheel Slip/Longitudinal Slip
    wheel_slip = fcn_VD_dtWheelSlip(Vx,Vy,yaw_rate,omega',steering_angle,vehicle);
    combined_slip = fcn_VD_dtCombinedSlip(slip_angle,wheel_slip); % combined slip

    %% Normal Forces
    normal_force = fcn_VD_dtNormalForce([ax;ay],vehicle,road_properties,...
        type_of_transfer);

    %% Tire Forces
    tire_force = fcn_VD_dtTireForceBrush(slip_angle,wheel_slip,normal_force,...
        friction_coefficient',vehicle);
    matlab_Fx(counter,:) = tire_force(:,1)';
    matlab_Fy(counter,:) = tire_force(:,2)';

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
fcn_VD_plotCompareLongitudinalTireForce(time,Fx,'Simulink',...
    matlab_time,matlab_Fx,'Matlab');
fcn_VD_plotCompareLateralTireForce(time,Fy,'Simulink',...
    matlab_time,matlab_Fy,'Matlab');
