function [approachXYPsi, curveXYPsi, exitXYPsi, hPlots] = fcn_VD_parkingPlotSolution(...
	referenceTurnPathXYPhi, initialXYPhi, ...
	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, varargin)

%% fcn_VD_parkingPlotSolution
%   This function plots the 2DOF bicycle model as a stick figure
%
% FORMAT:
%
%      [approachXYPsi, curveXYPsi, exitXYPsi] = fcn_VD_parkingPlotSolution(...
%    	referenceTurnPathXYPhi, initialXYPhi, ...
%   	approachPathDistance, endOfArcAngleInRadians, exitPathDistance, ...
%       (anglesToPlotTangents), (figNum))
%
% INPUTS:
%
%      referenceTurnPathXYPhi: an Nx3 vector of the [X, Y, yaw] positions
%      of the reference path used for turning. The reference path is
%      assumed to start at [0 0 0].
%
%      initialXYPhi: a 1x3 vector denoting the [X, Y, yaw] positions of the
%      initial condition of the vehicle, before the maneuver, in units of
%      [m m rad]
%
%      approachPathDistance: a 1x1 scalar denoting how far (or back) the
%      vehicle must move before starting the turning maneuver
%
%      endOfArcAngleInRadians: a 1x1 scalar denoting how much the angle of
%      the vehicle must change
%
%      exitPathDistance: a 1x1 scalar denoting how far (or back) the
%      vehicle must move after the turning arc, to complete the
%      parking
%
%      (OPTIONAL INPUTS)
%
%      anglesToPlotTangents: a Nx1 column of angles, in radians, to plot
%      tangents. If left empty, no tangents are plotted. Default is empty
%      ([]).
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
%      hPlots: a cell array of the handles of each plotted dataset
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingPlotSolution
%     for a full test suite.
%
% This function was written on 2026_10_28 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingPlotSolution
%
% 2026_10_28 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingPlotSolution
%   % * Used fcn_VD_plotStickFigure2DOFBicycleModel as starter
%
% 2026_10_30 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingPlotSolution
%   % * Modified to allow negative input angles

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
		narginchk(MAX_NARGIN-2,MAX_NARGIN);

		% Validate the reference path and maneuver parameters
		fcn_DebugTools_checkInputsToFunctions(referenceTurnPathXYPhi, '3columns_of_numbers',[2 3]);
		fcn_DebugTools_checkInputsToFunctions(initialXYPhi, '3column_of_numbers',[1 1]);
		assert(isscalar(approachPathDistance) && isscalar(endOfArcAngleInRadians) && isscalar(exitPathDistance));

	end
end

% Does user want to specify the anglesToPlotTangents input?
anglesToPlotTangents = []; % Default empty
if 6 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        anglesToPlotTangents = temp;
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

% Check the reference path for reversed situation
if endOfArcAngleInRadians>=0
	signCorrectedReferenceTurnPathXYPhi = referenceTurnPathXYPhi;
else
	signCorrectedReferenceTurnPathXYPhi = [referenceTurnPathXYPhi(:,1) -1*referenceTurnPathXYPhi(:,2) -1*referenceTurnPathXYPhi(:,3)];
end

% Calculate the absolute turn angle
absoluteAngleChange = abs(endOfArcAngleInRadians);



% Find the location in the reference turn path where the angle matches the
% end of arc angle
endOfArcIndex = find(abs(referenceTurnPathXYPhi(:,3))>absoluteAngleChange,1);
if isempty(endOfArcIndex)
	error('The reference path does not include angles with magnitudes up to: %.2f degrees. Must quit.',absoluteAngleChange*180/pi);
end
endOfArcIndex = endOfArcIndex-1;

curveXYPsiUnfixed = signCorrectedReferenceTurnPathXYPhi(1:endOfArcIndex,:);

% Check to see if the end of the reference arc is EXACTLY equal to the
% indexed value
if referenceTurnPathXYPhi(endOfArcIndex,3)<absoluteAngleChange
	% Interpolate the last X value in arc
	interpolatedPoint = fcn_INTERNAL_interpolateAtAngle(signCorrectedReferenceTurnPathXYPhi,endOfArcAngleInRadians);
	curveXYPsiUnRotated = [curveXYPsiUnfixed; interpolatedPoint];
else
	curveXYPsiUnRotated = curveXYPsiUnfixed;
end

% Rotate the unrotated points
rotatedArcXY = fcn_INTERNAL_XYrotateAndTranslate(curveXYPsiUnRotated(:,1:2), initialXYPhi(1,3), initialXYPhi(1,1:2));
rotatedArcPsi = curveXYPsiUnRotated(:,3)+initialXYPhi(1,3);
curveXYPsi = [rotatedArcXY rotatedArcPsi];



% Calculate approach and exit
firstAngle = curveXYPsi(1,3);
approachVector = [cos(firstAngle) sin(firstAngle)];
approachPoint = curveXYPsi(1,1:2)+approachPathDistance*(-1)*approachVector;
approachXYPsi = [approachPoint firstAngle; curveXYPsi(1,:)];

lastAngle = curveXYPsi(end,3);
exitVector = [cos(lastAngle) sin(lastAngle)];
exitPoint = curveXYPsi(end,1:2)+ exitPathDistance*(1)*exitVector;
exitXYPsi = [curveXYPsi(end,:); exitPoint lastAngle];

% % For debugging
% lastSegment = [curveXYPsi(end,1:2); curveXYPsi(end,1:2)+exitVector];
% plot(lastSegment(:,1),lastSegment(:,2),'-','Color',[0 0 0],'LineWidth',5,'DisplayName','Approach Path')

hPlots = cell(3+length(anglesToPlotTangents),1);

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
	grid on;
	axis equal;

	% Plot trajectories
	hPlots{1,1} = plot(approachXYPsi(:,1),approachXYPsi(:,2),'-','Color',[0 1 0],'LineWidth',5,'DisplayName','Approach Path');
	hPlots{2,1} = plot(curveXYPsi(:,1),curveXYPsi(:,2),'-','LineWidth',5,'DisplayName','Turn Path');
	hPlots{3,1} = plot(exitXYPsi(:,1),exitXYPsi(:,2),'-','Color',0.7*[0 1 0],'LineWidth',5,'DisplayName','Exit Path');

	% Plot tangents?
	if ~isempty(anglesToPlotTangents)
		allData = [approachXYPsi; curveXYPsi; exitXYPsi];
		rangeToPlot = max(max(allData(:,1:2))) - min(min(allData(:,1:2)));
		for ith_angle = 1:length(anglesToPlotTangents)
			thisAngle = anglesToPlotTangents(ith_angle);
			interpolatedPoint = fcn_INTERNAL_interpolateAtAngle(signCorrectedReferenceTurnPathXYPhi,thisAngle);
			unitVector = [cos(thisAngle) sin(thisAngle)];
			tangentPath = [unitVector*rangeToPlot; -unitVector*rangeToPlot] + repmat(interpolatedPoint(1,1:2),2,1);

			rotatedTangentPath = fcn_INTERNAL_XYrotateAndTranslate(tangentPath, initialXYPhi(1,3), initialXYPhi(1,1:2));

			hPlots{3+ith_angle,1} = plot(rotatedTangentPath(:,1),rotatedTangentPath(:,2),'-','Color',0.7*[1 1 0],'LineWidth',1,'DisplayName',sprintf('Tangent Path, %.0f deg in local xy frame',thisAngle*180/pi));

		end
	end
	

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


%% fcn_INTERNAL_interpolateAtAngle
function interpolatedPoint = fcn_INTERNAL_interpolateAtAngle(referenceTurnPathXYPhi,angleToInterpolate)
	% Interpolate the last X value in arc
	arcLastX = interp1(referenceTurnPathXYPhi(:,3),referenceTurnPathXYPhi(:,1),angleToInterpolate,'linear');
	arcLastY = interp1(referenceTurnPathXYPhi(:,3),referenceTurnPathXYPhi(:,2),angleToInterpolate,'linear');
	interpolatedPoint = [arcLastX arcLastY angleToInterpolate];
end % Ends fcn_INTERNAL_interpolateAtAngle

%% fcn_INTERNAL_XYrotateAndTranslate
function rotatedPointsXY = fcn_INTERNAL_XYrotateAndTranslate(XYPointsUnRotated, angleToRotateInRadians, translationXY)
Npoints = length(XYPointsUnRotated(:,1));
rotation = [cos(angleToRotateInRadians) sin(angleToRotateInRadians); -sin(angleToRotateInRadians) cos(angleToRotateInRadians)];
rotatedPointsXY = XYPointsUnRotated(:,1:2)*rotation + [translationXY(1,1)*ones(Npoints,1) translationXY(1,2)*ones(Npoints,1)] ;
end % Ends fcn_INTERNAL_XYrotateAndTranslate