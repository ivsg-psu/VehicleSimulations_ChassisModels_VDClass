function [approachXYPsi, curveXYPsi, exitXYPsi, ...
    approachPathDistance, changeInArcAngleInRadians, exitPathDistance] = fcn_VD_parkingCalcEntryArcExitPaths(...
	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
    startingPointInXYphi, goalPointInXYphi, varargin)

%% fcn_VD_parkingCalcEntryArcExitPaths
%   Given an initial pose in [X Y phi] with units of [m m rad], and a final
%   pose of [X Y phi] in same units, finds the "t" and "u" coordinates of
%   the projection of the final pose onto the initial pose. 
%
% FORMAT:
%
%      [approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingCalcEntryArcExitPaths(...
%    	standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle, ...
%       startingPointInXYphi, goalPointInXYphi, (figNum))
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
%      goalPointInXYphi: a 1x3 vector denoting the [X, Y, yaw] positions of
%      the final condition of the vehicle, after the maneuver, in units of
%      [m m rad]
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
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
%     See the script: script_test_fcn_VD_parkingCalcEntryArcExitPaths
%     for a full test suite.
%
% This function was written on 2026_10_29 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingCalcEntryArcExitPaths
%
% 2026_10_29 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingCalcEntryArcExitPaths
%   % * Used fcn_VD_parking+XYPhiToTU as starter

% TO-DO:
% - 2026_10_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 6; % The largest Number of argument inputs to the function
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

        % Validate the inputs
        fcn_DebugTools_checkInputsToFunctions(standardParkingTurnPath, '3columns_of_numbers',[2 3]);
        Npoints = length(standardParkingTurnPath(:,1));
        fcn_DebugTools_checkInputsToFunctions(tDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
        fcn_DebugTools_checkInputsToFunctions(uDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
		fcn_DebugTools_checkInputsToFunctions(startingPointInXYphi, '3column_of_numbers',[1 1]);
		fcn_DebugTools_checkInputsToFunctions(goalPointInXYphi, '3column_of_numbers',[1 1]);

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

%% Step 1: Find the x-axis intercept of the given pose
[tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(startingPointInXYphi, goalPointInXYphi, (-1));
changeInArcAngleInRadians = goalPointInXYphi(1,3) - startingPointInXYphi(1,3);

%% Step 2: Make sure this is in the feasible range
% Look up t-value for this goalPhiInRadians value
minimumTfeasible = interp1(standardParkingTurnPath(:,3),tDistancesEachAngle,changeInArcAngleInRadians,'linear','extrap');

% Look up U-value for this goalPhiInRadians value
maximumUfeasible = interp1(standardParkingTurnPath(:,3),uDistancesEachAngle,changeInArcAngleInRadians,'linear','extrap');

if tDistance<minimumTfeasible
	error('t-value is not feasible');
end

if uDistance>maximumUfeasible
	error('u-value is not feasible');
end

%% Step 3: Calculate the "approach" distance, the "exit" distance, and the arc angle

approachPathDistance = tDistance-minimumTfeasible;
exitPathDistance = maximumUfeasible-uDistance;
endOfArcAngleInRadians = changeInArcAngleInRadians;

% Find where the turn starts
unitStartVectorXYPsi = [cos(startingPointInXYphi(1,3)) sin(startingPointInXYphi(1,3)) 0];
turnStartXYPsi = startingPointInXYphi + unitStartVectorXYPsi*approachPathDistance;

[approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
    standardParkingTurnPath, turnStartXYPsi, ...
    approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
    ([]), (-1));


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
    % Plot the results
    anglesToPlotTangents = [startingPointInXYphi(1,3); goalPointInXYphi(1,3)] - startingPointInXYphi(1,3);
    [approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
    	standardParkingTurnPath, turnStartXYPsi, ...
    	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
    	(anglesToPlotTangents), (figNum));
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

