function rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, varargin)

%% fcn_VD_steeringApplyRateLimit
%   Fill in test steering cases for the front steering using typical
%   maneuvers
%
% FORMAT:
%
%      rateLimitedSteeringInputsVsTime = fcn_VD_steeringApplyRateLimit(steeringInputsVsTime, (rateLimitInDegPerSec), (figNum))
%
% INPUTS:
%
%      steeringInputsVsTime: A Nx2 vector of steering inputs to use for
%      applying the rate limit. The format is:
%      [timeInSeconds  steeringInputsInRadians]
%
%      (OPTIONAL INPUTS)
%
%      rateLimitInDegPerSec: a positive value of the rate limit to use, in
%      deg/sec. Default is 2 degrees per second of roadwheel angle. Of
%      note, handwheel angle actuators typically limit to about 300
%      degrees/second to 900 degrees/second of handwheel angle. the
%      steering gear ratio is usually about 15:1 to 18:1 from handwheel to
%      roadwheel angle. For highway driving, the roadwheel rate limit
%      should be about 1 to 3 deg/sec. For low speed (parking lot), this
%      should be about 20 to 35 deg/sec.
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%   rateLimitedSteeringInputsVsTime: A Nx2 vector of steering inputs that
%   are rate-limited. The format is:
%      [timeInSeconds  steeringInputsInRadians]
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_steeringApplyRateLimit
%     for a full test suite.
%
% This function was written on 2026_11_12
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_steeringApplyRateLimit
%
% 2026_11_12 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_steeringApplyRateLimit
%   % * First write of function,
%   % * Using fcn_VD_derivatives2DOFBicycleModel as starter

% TO-DO:
% - 2026_11_12 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
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
		narginchk(0,MAX_NARGIN);

		% % Check the X input to be sure it has 1 column, 3 rows
		% fcn_DebugTools_checkInputsToFunctions(X, '5column_of_numbers',[1 1]);
		% 
		% % Check the inputFrontRoadWheelAngleRadians input to be sure it has 1 column and 1 row
		% fcn_DebugTools_checkInputsToFunctions(inputFrontRoadWheelAngleRadians, '1column_of_numbers',[1 1]);
		% 
		% % Check the vehicleParameters input is a structure
		% assert(isstruct(vehicleParameters));
		% 
		% % Check the U input to be sure it has 1 column and 1 row, positive
		% fcn_DebugTools_checkInputsToFunctions(U, 'positive_1column_of_numbers',[1 1]);
	end
end


% Does the user want to specify the rateLimitInDegPerSec input?
rateLimitInDegPerSec = 2; % Default 
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)        
        rateLimitInDegPerSec = temp; % Set the value from user input
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

% Initialize output
rateLimitedSteeringInputsVsTime = steeringInputsVsTime; % [timeInSeconds  steeringInputsInRadians]

steeringTimes = steeringInputsVsTime(:,1);
Nsteering = length(steeringTimes);
for ith_time = 2:Nsteering
	thisDeltaT = steeringTimes(ith_time,1)-steeringTimes(ith_time-1,1);
	maxPossibleSteeringChangeInRadians = thisDeltaT*rateLimitInDegPerSec*pi/180;
	previousSteering = rateLimitedSteeringInputsVsTime(ith_time-1,2);
	maxSteering = previousSteering+maxPossibleSteeringChangeInRadians;
	minSteering = previousSteering-maxPossibleSteeringChangeInRadians;

	thisSteering = steeringInputsVsTime(ith_time,2);

	steeringToUse = thisSteering;
	if thisSteering>maxSteering
		steeringToUse = maxSteering;
	elseif thisSteering<minSteering
		steeringToUse = minSteering;
	end

	rateLimitedSteeringInputsVsTime(ith_time,2) = steeringToUse;
end

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

	% plot the derivative outputs as a vector
	figure(figNum);
	hold on;

	[xb, yb] = stairs(steeringInputsVsTime(:,1),steeringInputsVsTime(:,2)*180/pi);
	plot( xb, yb,'-','DisplayName','Input steering');

	[xb, yb] = stairs(rateLimitedSteeringInputsVsTime(:,1),rateLimitedSteeringInputsVsTime(:,2)*180/pi);
	plot( xb, yb,'-','DisplayName','Rate-limited steering (plotted as sample and hold)');

	plot(rateLimitedSteeringInputsVsTime(:,1),rateLimitedSteeringInputsVsTime(:,2)*180/pi,'-','DisplayName','Rate-limited steering (plotted with connections, to show rate)');


	legend('Interpreter','none');
	xlabel('Time [sec]','Interpreter','none');
	ylabel('Front steering input [deg]','Interpreter','none');
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

