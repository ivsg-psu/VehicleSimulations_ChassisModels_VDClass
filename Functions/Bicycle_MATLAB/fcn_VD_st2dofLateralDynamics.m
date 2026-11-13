function dXdt = fcn_VD_st2dofLateralDynamics(X, U, lateralForces, vehicleParameters, varargin)
%% fcn_VD_stLateralDynamics
%   This function calculates lateral acceleration and yaw acceleration for
%   the 2DOF single-track "bicycle" model.
%
% FORMAT:
%
%   dXdt = fcn_VD_stLateralDynamics(X, U, lateralForces, vehicleParameters, (figNum))
%
%   dXdt ~ [Vdot; rdot]
%
% INPUTS:
%
%      X: A 1x2 vector of lateral velocity and yaw rate [m/s rad/s]
%
%      U: Longitudinal velocity [m/s]
%
%      lateralForces: A 1x2 vector containing lateral-forces at front and
%      rear wheels [Newtons]
%
%      vehicleParameters: a structure with the following subfields
%          vehicleParameters.DistancesWheelbase
%          vehicleParameters.DistancesCGFromRearAxle
%          vehicleParameters.Mass
%          vehicleParameters.Izz
%          vehicleParameters.TireFrontCorneringCoefficient
%          vehicleParameters.TireRearCorneringCoefficient
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%        dXdt: A 1x2 vector of lateral acceleration and yaw acceleration [m/s^2 rad/s^2]
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_wheelVelocitiesBicycleModel
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_st2dofLateralDynamics
%     for a full test suite.
%
% This function was originally written on 2021/05/16 by Satya Prasad. It is
% maintained by S. Brennan. Questions or comments? sbrennan@psu.edu


% REVISION HISTORY:
%
% As: fcn_VD_stLateralDynamics
%
% 2021_05_16 by Satya Prasad, szm888@psu.edu
% - In fcn_VD_stLateralDynamics
%   % * First write of the function
%
% As: fcn_VD_st2dofLateralDynamics
%
% 2026_11_10 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_st2dofLateralDynamics
%   % * Renamed function
%   % * Used to be called fcn_VD_stLateralDynamics
%   % * Changed inputs and outputs to be column vectors, not row vectors

% TO-DO:
% - 2026_11_10 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)



%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 5; % The largest Number of argument inputs to the function
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
		narginchk(MAX_NARGIN-1,MAX_NARGIN);

		% Check the X input to be sure it has 1 column, 3 rows
		fcn_DebugTools_checkInputsToFunctions(X, '2column_of_numbers',[1 1]);

		% Check the U input to be sure it has 1 column and 1 row, positive
		fcn_DebugTools_checkInputsToFunctions(U, 'positive_1column_of_numbers',[1 1]);

		% Check the lateralForces input to be sure it has 1 column, 3 rows
		fcn_DebugTools_checkInputsToFunctions(lateralForces, '2column_of_numbers',[1 1]);

		% Check the vehicleParameters input is a structure
		assert(isstruct(vehicleParameters));

	end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
	temp = varargin{end};
	if ~isempty(temp) % Did the user NOT give an empty figure number?
		figNum = temp; %#ok<NASGU>
		flag_do_plots = 1;
	end
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

% Pull out the parameters
L = vehicleParameters.DistancesWheelbase;
b = vehicleParameters.DistancesCGFromRearAxle;
a = L - b;
m = vehicleParameters.Mass; % mass in kg
Izz = vehicleParameters.Izz; % z-axis moment of inertia

% Grab the states
r = X(2);
dVdt = (1/m)*(lateralForces(1)+lateralForces(2))-r*U;
drdt = (1/Izz)*(a*lateralForces(1)-b*lateralForces(2));
dXdt = [dVdt drdt];

%% Plot the results (for debugging)?
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
if flag_do_plots
	% Plot results
	% figure(figNum);
	% hold on;
	% axis equal;
	% grid on;
	%
	%
	% % Plot the tires
	% tireRadius = L/10;
	%
	% tireDimensions = [tireRadius 0; -tireRadius 0];
	%
	% % Plot front tire
	% angleRad = 0; %slipAnglesInRadians(1,1);
	% rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	% offsetTire = a*[ones(2,1) zeros(2,1)] + tireDimensions*rotation;
	% plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',20);
	%
	% % Plot rear tire
	% angleRad = 0; % slipAnglesInRadians(1,2);
	% rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	% offsetTire = -b*[ones(2,1) zeros(2,1)] + tireDimensions*rotation;
	% plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',20);
	%
	% % Plot the vehicle as a stick figure
	% plot([a -b],[0 0],'-','LineWidth',5);
	%
	% tireForceDirectionsInRadians = [90 90]*pi/180;
	%
	% % Plot the normalizedTireForces
	% normalizedTireForces = lateralForces./[frontDownForce rearDownForce];
	%
	% slipVectors = [normalizedTireForces(1,1) 0; 0 0];
	% angleRad = tireForceDirectionsInRadians(1,1);
	% rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	% offsetTire = a*[ones(2,1) zeros(2,1)] + slipVectors*rotation;
	% plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',1);
	%
	% slipVectors = [normalizedTireForces(1,2) 0; 0 0];
	% angleRad = tireForceDirectionsInRadians(1,2);
	% rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	% offsetTire = -b*[ones(2,1) zeros(2,1)] + slipVectors*rotation;
	% plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',1);
	%
	% % quiver(a, 0,wheelVelocities(1,1),wheelVelocities(1,2),0);
	% % quiver(-b, 0,wheelVelocities(2,1),wheelVelocities(2,2),0);
end

if flag_do_debug
	fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
