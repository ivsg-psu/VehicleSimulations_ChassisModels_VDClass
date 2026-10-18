function wheel_slip = fcn_VD_stWheelSlip(U,V,r,wheel_angular_velocity,...
                        steering_angle,vehicle)
%% fcn_VD_stWheelSlip
%   This function computes Wheel Slip/Longitudinal Slip for front and rear wheels.
%   It uses single-track/bicycle vehicle model.
%
%   Coordinate System: ISO
%
% FORMAT:
%
%   wheel_slip = fcn_VD_stWheelSlip(U,V,r,wheel_angular_velocity,...
%                   steering_angle,vehicle)
%
% INPUTS:
%
%   U: Longitudinal velocity [m/s]
%   V: Lateral velocity [m/s]
%   r: Yaw rate [rad/s]
%   wheel_angular_velocity: A 2x1 vector of wheel angular velocities [rad/s]
%   [Front; Rear]
%   steering_angle: A 2x1 vector of steering angles [rad]
%   [Front; Rear]
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   wheel_slip: A 2x1 vector of wheel slip [No Units]
%   [Front; Rear]
%
% This function was written on 2021/08/13 by Satya Prasad
% Questions or comments? szm888@psu.edu

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1, 'STARTING function: %s, in file: %s\n', st(1).name, st(1).file);
end

%% Check input arguments
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
if flag_check_inputs
    % Are there the right number of inputs?
    if 6~=nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the inputs
    fcn_VD_checkInputsToFunctions(U,'non negative');
    fcn_VD_checkInputsToFunctions(V,'number');
    fcn_VD_checkInputsToFunctions(r,'number');
    fcn_VD_checkInputsToFunctions(wheel_angular_velocity,'vector2');
    fcn_VD_checkInputsToFunctions(steering_angle,'vector2');
end

%% Calculate Longitudinal Slip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wheel_slip = nan(2,1); % Initialize a variable

wheel_velocity = fcn_VD_stWheelVelocity(U,V,r,vehicle); % calculate wheel velocities
wheel_longitudinal_velocity = wheel_velocity(:,1).*cos(steering_angle)+...
                              wheel_velocity(:,2).*sin(steering_angle); % velocity of wheel along the wheel direction
wheel_longitudinal_velocity(0==wheel_longitudinal_velocity) = eps;

% Wheel-Slip calculation
wheel_slip(0>wheel_longitudinal_velocity) = 0;
% wheel slip is calculated only when the tire is moving in the direction it
% is pointed
target_indices = 0<wheel_longitudinal_velocity;
wheel_slip(target_indices) = (vehicle.Re.*wheel_angular_velocity(target_indices)-...
                              wheel_longitudinal_velocity(target_indices))./...
                              wheel_longitudinal_velocity(target_indices);

% Saturation
wheel_slip(1<wheel_slip)  = 1; % maximum wheel slip
wheel_slip(-1>wheel_slip) = -1; % minimum wheel slip

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end