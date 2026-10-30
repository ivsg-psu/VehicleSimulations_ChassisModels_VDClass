function [standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle] = ...
	fcn_VD_parkingFillStandardParkingTurnPath(varargin)

%% fcn_VD_parkingFillStandardParkingTurnPath
%   Given an initial pose in [X Y phi] with units of [m m rad], and a final
%   pose of [X Y phi] in same units, finds the "t" and "u" coordinates of
%   the projection of the final pose onto the initial pose. 
%
% FORMAT:
%
%      [standardParkingTurnPath, tDistancesEachAngle, uDistancesEachAngle] = ...
%      fcn_VD_parkingFillStandardParkingTurnPath((figNum))
%
% INPUTS:
%
%      (none)
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
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
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_parkingXYPhiToTU
%      fcn_Path_findProjectionHitOntoPath
%      fcn_VD_parkingPlotSolution
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingFillStandardParkingTurnPath
%     for a full test suite.
%
% This function was written on 2026_10_30 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingFillStandardParkingTurnPath
%
% 2026_10_30 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingFillStandardParkingTurnPath
%   % * Used fcn_VD_parkingCalcEntry+ArcExitPaths as starter

% TO-DO:
% - 2026_10_30 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 1; % The largest Number of argument inputs to the function
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

        % % Validate the inputs
        % fcn_DebugTools_checkInputsToFunctions(standardParkingTurnPath, '3columns_of_numbers',[2 3]);
        % Npoints = length(standardParkingTurnPath(:,1));
        % fcn_DebugTools_checkInputsToFunctions(tDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
        % fcn_DebugTools_checkInputsToFunctions(uDistancesEachAngle, '1columns_of_numbers',[Npoints Npoints]);
		% fcn_DebugTools_checkInputsToFunctions(startingPointInXYphi, '3column_of_numbers',[1 1]);
		% fcn_DebugTools_checkInputsToFunctions(goalPointInXYphi, '3column_of_numbers',[1 1]);

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
dataFileName = fullfile(pwd,'Data','parkingData.mat');
if 1==1 && exist(dataFileName,'file')
	load(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle'); 
else

	standard_parking_steering_amplitude_degrees = 20; % 2 degrees of steering amplitude for input sinewave
	standard_parking_steering_amplitude_radians = standard_parking_steering_amplitude_degrees*pi/180;

	% Set the simulation time/state arguments
	initialStates = [0 0 0]; % [X Y phi] in [m],[m],[rad]
	deltaT = 0.01; % Units are [sec]
	startTime = 0;
	arcLength = 170*pi/180*predictedTurnRadius;
	endTime = arcLength/vehicleParameters.U;
	timeInterval = [startTime endTime];  % Units are [sec]

	simulationTimes = (startTime:deltaT:endTime)';
	inputsVsTime = [simulationTimes standard_parking_steering_amplitude_radians*ones(length(simulationTimes(:,1)),1)]; % [times steering angles]

    % Call a function to generate a standard parking path
	
	figure(figNum);clf;
	[standardParkingTurnPath, ~, ~] = ...
		fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
		timeInterval, inputsVsTime, vehicleParameters, (figNum));

	% Call a function to convert the standard path into t- and u-distances
	[tDistancesEachAngle, uDistancesEachAngle] = fcn_VD_parkingSolveUTFromStandardTurn(...
		standardParkingTurnPath, (-1));

	save(dataFileName,'standardParkingTurnPath','tDistancesEachAngle','uDistancesEachAngle');

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
    % Plot the results
    fcn_VD_parkingSolveUTFromStandardTurn(standardParkingTurnPath, (figNum));
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

