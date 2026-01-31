function dydt = fcn_VD_stLateralDynamics(~, y, U, lateral_forces, vehicle)
%% fcn_VD_stLateralDynamics
%   This function calculates lateral acceleration and yaw acceleration.
%   It uses single-track/bicycle model
%
% FORMAT:
%
%   dydt = fcn_VD_stLateralDynamics(~, y, U, Fy, vehicle)
%   dydt ~ [Vdot; rdot]
%
% INPUTS:
%
%   y: A 2x1 vector of lateral velocity and yaw rate [m/s; rad/s]
%   U: Longitudinal velocity [m/s]
%   lateral_forces: A 2x1 vector containing lateral-forces at front and rear wheels [Newton]
%   vehicle: MATLAB structure containing vehicle properties
%
% OUTPUTS:
%
%   dydt: A 2x1 vector of lateral acceleration and yaw acceleration [m/s^2 rad/s^2]
%
% This function was written on 2021/05/16 by Satya Prasad
% Questions or comments? szm888@psu.edu
%

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
    

	% Check the y input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(y, '2column_of_numbers',[1 2]);

	% Check the U input to be sure it has 1 col, 1 row, positive
	fcn_DebugTools_checkInputsToFunctions(U, 'positive_1column_of_numbers',[1 1]);

	% Check the lateral_forces input to be sure it has 2 col, 1+ rows
	fcn_DebugTools_checkInputsToFunctions(lateral_forces, '2column_of_numbers',[1 2]);

end

%% Calculate Lateral acceleration and Yaw acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = y(2);
dVdt = (1/vehicle.m)*(lateral_forces(1)+lateral_forces(2))-r*U;
drdt = (1/vehicle.Iz)*(vehicle.a*lateral_forces(1)-vehicle.b*lateral_forces(2));
dydt = [dVdt; drdt];

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