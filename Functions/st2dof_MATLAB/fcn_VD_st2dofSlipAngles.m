function slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, varargin)

%% fcn_VD_st2dofSlipAngles
%   This function computes slip-angles for front and rear wheels. It
%   uses a single-track/bicycle model.
%
% FORMAT:
%
%      slipAnglesInRadians = fcn_VD_st2dofSlipAngles(U, V, r, steeringAngles, vehicleParameters, (flagUseLinear), (figNum))
%
% INPUTS:
%
%      U: x-axis velocity, in m/s
%
%      V: y-axis velocity, in m/s
%
%      r: yawRate of the vehilce in rad/sec
%
%      steeringAngles: a 1x2 vector denoting 
%      [frontSteering rearSteering] 
%      in units of [rad]. This is the road wheel angle.
%
%      vehicleParameters: a structure with the following subfields
%          vehicleParameters.DistancesWheelbase
%          vehicleParameters.DistancesCGFromRearAxle
%            
%      (OPTIONAL INPUTS)
%
%      flagUseLinear: a flag to set the calculation format for slip
% 
%          0: (default) uses inverse tangent calculation
% 
%          1: uses linear division calculation
%    
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      slipAnglesInRadians: An 1x2 vector of the slip angles of the front
%      and rear tires in rows 1 and 2 respectively, in radians
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_st2dofWheelVelocities
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_st2dofSlipAngles
%     for a full test suite.
%
% This function was written on 2026_10_22 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_slipAnglesBicycleModel
%
% 2026_10_22 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_slipAnglesBicycleModel
%   % * Used fcn_VD_kinematicBicycleModelRK4 as starter
%   % * Copied material over from fcn_VD_stSlipAngle
%
% As: fcn_VD_st2dofSlipAngles
%
% 2026_11_10 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_st2dofSlipAngles
%   % * Renamed function
%   % * Used to be called fcn_VD_slipAnglesBicycleModel
%
% 2026_11_11 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_st2dofSlipAngles
%   % * Updated header to match correct vehicleParameters structure
%   % * Changed inputs and outputs to 1x2 format, not 2x1
%   % * Added intermediate variable for debugging, velocityAnglesInRadians
%   % * Added flagUseLinear optional input to allow linear or nonlinear calculations

% TO-DO:
% - 2026_10_22 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 7; % The largest Number of argument inputs to the function
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

        % Validate that the U input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(U, '1column_of_numbers',[1 1]);

        % Validate that the V input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(V, '1column_of_numbers',[1 1]);

        % Validate that the r input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(r, '1column_of_numbers',[1 1]);

        % Validate that the steeringAngles input has 2 columns, 1 row
        fcn_DebugTools_checkInputsToFunctions(steeringAngles, '2column_of_numbers',[1 1]);

		% Check the vehicleParameters input is a structure
        assert(isstruct(vehicleParameters));

    end
end

%  Does the user want to specify flagUseLinear input?
flagUseLinear = 0; % Default is to NOT use linear
if 6 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        % Set the flagUseLinear value
        flagUseLinear = 1;
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; 
        flag_do_plots = 1;
    end
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
wheelVelocities         = fcn_VD_st2dofWheelVelocities( U, V, r, vehicleParameters, (-1));
if 1==flagUseLinear
	velocityAnglesInRadians = (wheelVelocities(:,2)./wheelVelocities(:,1))';
else
	velocityAnglesInRadians = (atan(wheelVelocities(:,2)./wheelVelocities(:,1)))';
end
slipAnglesInRadians     = velocityAnglesInRadians - steeringAngles; % Slip angle

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
	figure(figNum);
	hold on;
	axis equal;
	grid on;


	% Plot the tires
	L = vehicleParameters.DistancesWheelbase;
	b = vehicleParameters.DistancesCGFromRearAxle;
	a = L-b;

	tireRadius = L/10;
	
	tireDimensions = [tireRadius 0; -tireRadius 0];	

	% Plot front tire
	angleRad = steeringAngles(1);
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	offsetTire = a*[ones(2,1) zeros(2,1)] + tireDimensions*rotation;
	plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',20);

	% Plot rear tire
	angleRad = steeringAngles(2);
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	offsetTire = -b*[ones(2,1) zeros(2,1)] + tireDimensions*rotation;
	plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',20);

	% Plot the vehicle as a stick figure
	plot([a -b],[0 0],'-','LineWidth',5);
	

	% Plot the slip angles
	slipVectors = [tireRadius*1.5 0; 0 0];	
	angleRad = slipAnglesInRadians(1);
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	offsetTire = a*[ones(2,1) zeros(2,1)] + slipVectors*rotation;
	plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',1);

	slipVectors = [tireRadius*1.5 0; 0 0];	
	angleRad = slipAnglesInRadians(2);
	rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
	offsetTire = -b*[ones(2,1) zeros(2,1)] + slipVectors*rotation;
	plot(offsetTire(:,1), offsetTire(:,2),'-','LineWidth',1);

	% quiver(a, 0,wheelVelocities(1,1),wheelVelocities(1,2),0);
	% quiver(-b, 0,wheelVelocities(2,1),wheelVelocities(2,2),0);
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


