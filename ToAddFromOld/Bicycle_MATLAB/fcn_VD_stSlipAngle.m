function slip_angles = fcn_VD_stSlipAngle(U, V, r, steering_angle, vehicle)
%% fcn_VD_stSlipAngle
%   This function computes slip-angles for front and rear wheels. It
%   uses a single-track/bicycle model.
%
%   Coordinate System: ISO
%
% FORMAT:
%
%   slip_angles = fcn_VD_stSlipAngle(U, V, r, steering_angle, vehicle)
%
% INPUTS:
%
%   U: Longitudinal velocity [m/s]
%   V: Lateral velocity [m/s]
%   r: Yaw rate [rad/s]
%   steering_angle: A 2x1 vector of steering angles [rad] [Front; Rear]
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   slip_angles: A 2x1 vector containing slip-angles of front and rear wheels [rad]
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
    if 5 ~= nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the inputs
    fcn_VD_checkInputsToFunctions(U,'non negative');
    fcn_VD_checkInputsToFunctions(V,'number');
    fcn_VD_checkInputsToFunctions(r,'number');
    fcn_VD_checkInputsToFunctions(steering_angle,'vector2');
end

%% Calculate Slip-Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wheel_velocity = fcn_VD_stWheelVelocity(U,V,r,vehicle); % calculate wheel velocities
slip_angles    = atan(wheel_velocity(:,2)./wheel_velocity(:,1))-steering_angle; % Slip angle

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