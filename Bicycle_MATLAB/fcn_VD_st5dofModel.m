function dydt = fcn_VD_st5dofModel( t, y,...
    steering_amplitude ,Period, wheel_torque,...
    vehicle, road_properties, friction_coefficient, type_of_transfer, varargin)

%% fcn_VD_st5dofModel
%   fcn_VD_st5dofModel sets-up a 5-DoF Single-Track vehicle model using
%   Brush tire model. This works without a controller but with
%   time-dependent inputs.
%
%
% FORMAT:
%
%       dydt = fcn_VD_st5dofModel( t, y,...
%       steering_amplitude, Period, wheel_torque,...
%       vehicle, road_properties, friction_coefficient, type_of_transfer, (figNum))
% with
%       dydt ~ [dUdt; dVdt; drdt; domegadt1 (to 4); dXdt; dYdt; dPhidt];
%
% INPUTS:
%
%     t: A number indicating time corresponding to y.
%
%     y: A 10x1 vector of velocities and pose consisting of:
%           [U; 
%            V; 
%            r; 
%            omega1; 
%            omega2; 
%            omega3; 
%            omega4; 
%            X; 
%            Y; 
%            Phi]
% 
%     steering_amplitude: the amplitude of steering input
%
%     Period: the period of the sinewave steering input
%
%     wheel_torque: the input torque to the wheels
%
%     vehicle: MATLAB structure containing vehicle properties.
%
%     road_properties: MATLAB structure containing road properties.
%
%     friction_coefficient: A 2x1 vector of friction coefficients.
%     [Front; Rear]
%
%     type_of_transfer: To decide type of load transfer, as one of
%     following:
%         'longitudinal': Only longitudinal weight transfer
%         'default': No weight transfer. Any non-matching string argument
%                    will give this result.
%
%     (OPTIONAL INPUTS)
%
%     figNum: a figure number to plot results.
%
% OUTPUTS:
%
%      dydt: A 8x1 vector of accelerations and velocities.
%           [longitudinal acceleration; 
%            lateral acceleration; 
%            yaw acceleration;
%            wheel angular acceleration_Front; 
%            wheel_angular_acceleration_Rear] 
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_st5dofModel
%     for a full test suite.
%
% This function was written on 2021_08_13 by Satya Prasad and is maintained
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% 2021_08_13 by Satya Prasad, szm888@psu.edu
% - Wrote the code originally
% 
% 2025_12_30 by Sean Brennan, sbrennan@psu.edu
% - Fixed header formatting to standard form
% - Fixed input checking to standard form

% TO-DO:
%
% 2025_11_21 by Sean Brennan, sbrennan@psu.edu
% - Need to update the output argument list correctly (with clarity)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 10; % The largest Number of argument inputs to the function
flag_max_speed = 0; % The default. This runs code with all error checking
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_VD_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_VD_FLAG_CHECK_INPUTS");
    MATLABFLAG_VD_FLAG_DO_DEBUG = getenv("MATLABFLAG_VD_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VD_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VD_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VD_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VD_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug % If debugging is on, print on entry/exit to the function
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_figNum = 999978; %#ok<NASGU>
else
    debug_figNum = []; %#ok<NASGU>
end

%% check input arguments?
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
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(9,MAX_NARGIN);

        % Check the t input to be sure it has 1 column and 1 row, positive
        fcn_DebugTools_checkInputsToFunctions(y, 'positive_1column_of_numbers',[1 1]);

        % Check the y input to be sure it has 1 column, 8 rows
        fcn_DebugTools_checkInputsToFunctions(y, '1column_of_numbers',[8 8]);

        % Check the friction_coefficient input to be sure it has 1 column, 2 rows
        fcn_DebugTools_checkInputsToFunctions(friction_coefficient, '1column_of_numbers',[2 2]);

        % Check the type_of_transfer input to be sure it has 1 column, 2 rows
        fcn_DebugTools_checkInputsToFunctions(type_of_transfer, 'URHERE',[2 2]);

        % % Check the inputs
        % fcn_VD_checkInputsToFunctions(t,'non negative');
        % fcn_VD_checkInputsToFunctions(y,'vector8');
        % fcn_VD_checkInputsToFunctions(friction_coefficient,'vector2');
        % fcn_VD_checkInputsToFunctions(type_of_transfer,'string');

    end
end


% 
%   Set the start values
% [flag_start_is_a_point_type, start_zone_definition] = fcn_Laps_checkZoneType(start_zone_definition, 'start_definition', -1);
% 
% 
%   The following area checks for variable argument inputs (varargin)
% 
%   Does the user want to specify the end_definition?
%   Set defaults first:
% end_zone_definition = start_zone_definition; % Default case
% flag_end_is_a_point_type = flag_start_is_a_point_type; % Inheret the start case
%   Check for user input
% if 3 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         % Set the end values
%         [flag_end_is_a_point_type, end_zone_definition] = fcn_Laps_checkZoneType(temp, 'end_definition', -1);
%     end
% end
% 
%   Does the user want to specify excursion_definition?
% flag_use_excursion_definition = 0; % Default case
% flag_excursion_is_a_point_type = 1; % Default case
% if 4 <= nargin
%     temp = varargin{2};
%     if ~isempty(temp)
%         % Set the excursion values
%         [flag_excursion_is_a_point_type, excursion_definition] = fcn_Laps_checkZoneType(temp, 'excursion_definition',-1);
%         flag_use_excursion_definition = 1;
%     end
% end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; %#ok<NASGU>
        flag_do_plots = 1;
    end
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

% flag_update: A global variable to decide the value of acceleration at a
% time-step 't'
% global_acceleration: A global variable that stores accelerations at a
% time-step 't'
global flag_update global_acceleration % change the variable name based on the matlab script
persistent delayed_acceleration % variable to store acceleration in the previous time-step

   
URHERE

%% Implement 5-DoF Vehicle Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = y(1); % longitudinal velocity [m/s]
V = y(2); % lateral velocity [m/s]
r = y(3); % yaw rate [rad/s]
omega = y(4:5); % angular velocity of wheels [rad/s]
pose  = y(6:8); % position and orientation of the vehicle

%% Steering Input
delta_f = (1-1*(0<t-Period))*steering_amplitude*sin((2*pi/Period)*t); % front steering angle [rad]
steering_angle = [delta_f; 0];

%% Slips
% Slip Angle/Lateral Slip
slip_angle = fcn_VD_stSlipAngle(U,V,r,steering_angle,vehicle);
% Wheel Slip/Longitudinal Slip
wheel_slip = fcn_VD_stWheelSlip(U,V,r,omega,steering_angle,vehicle);

%% Normal Forces
if flag_update
    % load transfer is calculated based on the acceleration in the previous
    % time-step 't-deltaT'
    delayed_acceleration = global_acceleration;
end
normal_force = fcn_VD_stNormalForce(delayed_acceleration(1:2),vehicle,...
    road_properties,type_of_transfer);

%% Tire Forces
tire_force = fcn_VD_stTireForceBrush(slip_angle,wheel_slip,normal_force,...
    friction_coefficient,vehicle);

%% Newtonian Dynamics
acceleration = fcn_VD_st5dofForceEquation(tire_force,wheel_torque,...
    steering_angle,vehicle,road_properties);
DvelDt = fcn_VD_5dofStateEquation(t,y(1:5),acceleration);
if flag_update
    % 'global_acceleration' will be updated only in the first call to the
    % function in RK4 method i.e., while computing k1 in 'fcn_VD_RungeKutta'
    global_acceleration = acceleration;
    flag_update = false;
end

%% Body2Global coordinates
DposeDt = fcn_VD_Body2GlobalCoordinates(t,pose,U,V,r);

%% Write to output
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