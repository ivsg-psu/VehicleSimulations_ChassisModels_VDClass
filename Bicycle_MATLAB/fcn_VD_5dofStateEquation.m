function dydt = fcn_VD_5dofStateEquation(~,y,acceleration)
%% fcn_VD_5dofStateEquation
%   This function is a differential equation for longitudinal velocity,
%   lateral velocity, yaw rate, and angular velocities of wheels.
%
% FORMAT:
%
%   dydt = fcn_VD_5dofStateEquation(~,y,acceleration)
%   dydt ~ [dUdt; dVdt; drdt; domegadt];
%
% INPUTS:
%
%   y: A 5x1 vector of velocities in body coordinates.
%   [longitudinal velocity; lateral velocity; yaw rate; wheel angular
%   velocity[Front; Rear]]
%   acceleration: A 5x1 vector of accelerations.
%   [longitudinal acceleration; lateral acceleration; yaw acceleration;
%   wheel angular acceleration[Front; Rear]]
%
% OUTPUTS:
%
%   dydt: A 5x1 vector of accelerations.
%   [longitudinal acceleration; lateral acceleration; yaw acceleration;
%   wheel angular acceleration[Front; Rear]]
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
    if 3~=nargin
        error('Incorrect number of input arguments.')
    end
    
    % Check the inputs
    fcn_VD_checkInputsToFunctions(y,'vector5');
    fcn_VD_checkInputsToFunctions(acceleration,'vector5');
end

%% Differential equation for 5-DOF model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dUdt = acceleration(1)+y(3)*y(2); % Vxdot
dVdt = acceleration(2)-y(3)*y(1); % Vydot
drdt = acceleration(3); % yaw acceleration

domegadt = acceleration(4:5); % angular acceleration of wheel

dydt = [dUdt; dVdt; drdt; domegadt];

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