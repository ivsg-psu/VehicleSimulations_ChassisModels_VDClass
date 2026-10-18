function dydt = fcn_VD_bicycle2dofModel(t,y,steering_amplitude,Period,U,vehicle)
%% fcn_VD_bicycle2dofModel
%   This function uses Bicycle model with Linear Tire model
%
% FORMAT:
%
%   dydt = fcn_VD_bicycle2dofModel(t,y,steering_amplitude,Period,U,vehicle)
%
% INPUTS:
%
%   y: A 5x1 vector of velocities and global pose [V; r; X; Y; Phi]
%   steering_amplitude:
%   Period:
%   U: Longitudinal velocity [m/s]
%   vehicle:
%
% OUTPUTS:
%
%   dydt: A 5x1 vector of accelerations and velocities
%
% This function was written on 2021_07_09 by Satya Prasad
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
    fcn_VD_checkInputsToFunctions(t,'non negative');
    fcn_VD_checkInputsToFunctions(y,'vector5');
end

%% Implements Bicycle Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = y(1); % lateral velocity [m/s]
r = y(2); % yaw rate [rad/s]
pose = y(3:5); % pose of the vehicle

%% Inputs
delta_f = (1-1*(0<t-Period))*(pi/180)*steering_amplitude*sin((2*pi/Period)*t);
steering_angle = [delta_f; 0];

%% Slips
slip_angles = fcn_VD_stSlipAngle(U,V,r,steering_angle,vehicle);

%% Tire Forces
lateral_forces = fcn_VD_stTireForceLinear(slip_angles, vehicle);

%% Newtonian Dynamics
DvelDt = fcn_VD_stLateralDynamics(t, [V; r], U, lateral_forces, vehicle);

%% Body2Global Coordinates
DposeDt = fcn_VD_Body2GlobalCoordinates(t, pose, U, V, r);

%% Output
dydt = [DvelDt; DposeDt];

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