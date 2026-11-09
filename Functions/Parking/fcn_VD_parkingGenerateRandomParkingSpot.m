function [goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance]  = ...
	fcn_VD_parkingGenerateRandomParkingSpot(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, varargin)

%% 
%   Generates a random parking spot to use as a target. Allows user to set
%   flags to specify difficulty.
%
% FORMAT:
%
%      [goalParkingPointInXYphi, approachXYPsi, curveXYPsi, exitXYPsi , approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingGenerateRandomParkingSpot(...
%    	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
%       startingPointInXYphi, (flagDifficulty), (inputRange), (roundingLevel), (figNum))
%
% INPUTS:
%
%      standardParkingTurnPath: an Nx3 vector of the [X, Y, yaw] positions
%      of the reference path used for turning. The reference path is
%      assumed to start at [0 0 0], with units of [m m rad].
%
%      tDistancesEachAngle: a Nx1 scalar of the distances where the  angle
%      of each standardParkingTurnPath would intersect the local x-axis of
%      the standardParkingTurnPath, as measured in the direction of the
%      x-axis
%
%      uDistancesEachAngle: a Nx1 scalar of the distances where the  angle
%      of each standardParkingTurnPath would intersect the local x-axis of
%      the standardParkingTurnPath, as measured in the direction of the
%      sampled point
%
%      startingPointInXYphi: a 1x3 vector denoting the [X, Y, yaw]
%      positions of the initial condition of the vehicle, before the
%      maneuver, in units of [m m rad]
%
%      (OPTIONAL INPUTS)
%
%      flagDifficulty: a setting to allow different difficulty levels
%
%         0: (default) only positive steering, entry, and exit distances
%
%         1: positive or negative steering, positive entry and exit distances
%
%         2: positive or negative steering, entry, and exit distances
%
%      inputRange: the  range to use in the format:
%        [approachDistanceMin approachDistanceMax...
%         deltaYawMin delatYawMax ...
%         exitDistanceMin exitDistanceMax] 
%      This is the range to randomly search. Default is [0 10 45*pi/180
%      90*pi/180 0 10]
%
%      roundingLevel: the number of decimal places to use when searching
%      for final values. The default is 2 decimal places. The final
%      solutions are based on using the rounded values.
%      
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      goalParkingPointInXYphi: a 1x3 vector denoting the [X, Y, yaw] positions of
%      the final condition of the vehicle, after the maneuver, in units of
%      [m m rad]
%
%      approachXYPsi: an Nx3 vector of the [X, Y, yaw] positions
%      of the approach path used for parking
%
%      curveXYPsi: an Nx3 vector of the [X, Y, yaw] positions
%      of the turning curved path used for parking
%
%      exitXYPsi: an Nx3 vector of the [X, Y, yaw] positions
%      of the exit path used for parking
%
%      approachPathDistance: a 1x1 scalar denoting how far (or back) the
%      vehicle must move before starting the turning maneuver
%
%      changeInArcAngleInRadians: a 1x1 scalar denoting how much the angle of
%      the vehicle must change
%
%      exitPathDistance: a 1x1 scalar denoting how far (or back) the
%      vehicle must move after the turning arc, to complete the
%      parking
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_parkingXYPhiToTU
%      fcn_Path_findProjectionHitOntoPath
%      fcn_VD_parkingPlotSolution
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingGenerateRandomParkingSpot
%     for a full test suite.
%
% This function was written on 2026_10_29 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% 2026_11_08 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingGenerateRandomParkingSpot
%   % * First write of function
%   % * Used fcn_VD_parkingCalcEntryArcExitPaths as starter

% TO-DO:
% - 2026_11_08 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 8; % The largest Number of argument inputs to the function
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
		narginchk(MAX_NARGIN-4,MAX_NARGIN);

        % Validate the inputs
        fcn_DebugTools_checkInputsToFunctions(standardParkingTurnPath, '3columns_of_numbers',[2 3]);
        Npoints = length(standardParkingTurnPath(:,1));
        fcn_DebugTools_checkInputsToFunctions(tDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
        fcn_DebugTools_checkInputsToFunctions(uDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
		fcn_DebugTools_checkInputsToFunctions(startingPointInXYphi, '3column_of_numbers',[1 1]);

	end
end


% Does user want to specify the flagDifficulty input?
flagDifficulty = 0; % Default
if 5 <= nargin
	temp = varargin{1};
	if ~isempty(temp)
		flagDifficulty = temp;
		if 0 == flag_max_speed
			if flag_check_inputs == 1
				% Check the flagDifficulty input, make sure it is a
				% positive integer, 1 row and 1 column
				fcn_DebugTools_checkInputsToFunctions(flagDifficulty,'nonzero_1columns_of_integers',[1 1]);
			end
		end
	end
end


% Does user want to specify the inputRange input?
inputRange = [0 10 45*pi/180 90*pi/180 0 10]; % Default
if 6 <= nargin
	temp = varargin{2};
	if ~isempty(temp)
		inputRange = temp;
		if 0 == flag_max_speed
			if flag_check_inputs == 1
				% Check the inputRange input, make sure it is characters
				fcn_DebugTools_checkInputsToFunctions(inputRange,'6columns_of_numbers', [1 1]);
			end
		end
	end
end


% Does user want to specify the roundingLevel input?
roundingLevel = 2; % Default
if 1 <= nargin
	temp = varargin{1};
	if ~isempty(temp)
		roundingLevel = temp;
		if 0 == flag_max_speed
			if flag_check_inputs == 1
				% Check the roundingLevel input, make sure it is a
				% positive integer, 1 row and 1 column
				fcn_DebugTools_checkInputsToFunctions(roundingLevel,'nonzero_1columns_of_integers',[1 1]);
			end
		end
	end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
figNum = [];
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; 
        flag_do_plots = 1;
    end
end

if isempty(figNum)
	figNum = gcf;
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

% Generate raw values
approachDistance = fcn_INTERNAL_generateRandomValueInRange(inputRange(1,1:2), roundingLevel);
changeInAngle = fcn_INTERNAL_generateRandomValueInRange(inputRange(1,3:4), roundingLevel);
exitDistance = fcn_INTERNAL_generateRandomValueInRange(inputRange(1,5:6), roundingLevel);

startOfCurve = zeros(1,3);
startingYawAngle = startingPointInXYphi(3);
deltaXY = approachDistance*[cos(startingYawAngle) sin(startingYawAngle)];
startOfCurve(1,1:2) = startingPointInXYphi(1,1:2)+deltaXY;
startOfCurve(1,3) = startingYawAngle;

% Generate a solution
[~, ~, exitXYPsi, ~] = fcn_VD_parkingPlotSolution(...
	standardParkingTurnPath, startOfCurve, ...
	approachDistance, changeInAngle, exitDistance, ...
	([]), (-1));

goalParkingPointInXYphi = exitXYPsi(end,:);

[approachXYPsi, curveXYPsi, exitXYPsi,...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
    standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalParkingPointInXYphi, (-1));



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
	figure(figNum);

	% Plot the results
	fcn_VD_parkingCalcEntryArcExitPaths(...
		standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
		startingPointInXYphi, goalParkingPointInXYphi, (figNum));

	% Approximately a Ford Taurus
	clear vehicleParameters
	vehicleParameters(1).a          = 0.9271; % Distance from front axle to CG, in meters
	vehicleParameters(1).b          = 1.5621;  % Distance from rear axle to CG, in meters

	fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startingPointInXYphi, [0; 0], true, (figNum));
	fcn_VD_parkingPlotParkingSpot(vehicleParameters, goalParkingPointInXYphi, (true), (figNum));

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

%% fcn_INTERNAL_generateRandomValueInRange
function roundedValue = fcn_INTERNAL_generateRandomValueInRange(range, rounding)
rawValue = (range(2)-range(1))*rand + range(1);
roundedValue = round(rawValue,rounding);
if 0==roundedValue
	roundedValue = 10^(-1*rounding);
end
end % Ends fcn_INTERNAL_generateRandomValueInRange