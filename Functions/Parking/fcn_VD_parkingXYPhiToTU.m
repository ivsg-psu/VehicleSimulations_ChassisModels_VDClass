function [tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(...
	initialXYPhi, finalXYPhi, varargin)

%% fcn_VD_parkingXYPhiToTU
%   Given an initial pose in [X Y phi] with units of [m m rad], and a final
%   pose of [X Y phi] in same units, finds the "t" and "u" coordinates of
%   the projection of the final pose onto the initial pose. 
%
% FORMAT:
%
%      [tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(initialXYPhi, finalXYPhi, (figNum))
%
% INPUTS:
%
%      initialXYPhi: a 1x3 vector denoting the [X, Y, yaw] positions of the
%      initial condition of the vehicle, before the maneuver, in units of
%      [m m rad]
%
%      finalXYPhi: a 1x3 vector denoting the [X, Y, yaw] positions of the
%      final condition of the vehicle, after the maneuver, in units of
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
%      tDistance: a 1x1 scalar of the distance where the finalXYPhi would
%      intersect the initialXYPhi, as measured in the direction of
%      initialXYPhi
%
%      uDistance: a 1x1 scalar of the distance where the finalXYPhi would
%      intersect the initialXYPhi, as measured in the direction of
%      finalXYPhi
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findProjectionHitOntoPath
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingXYPhiToTU
%     for a full test suite.
%
% This function was written on 2026_10_28 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingXYPhiToTU
%
% 2026_10_28 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingXYPhiToTU
%   % * Used fcn_VD_parking+PlotSolution as starter

% TO-DO:
% - 2026_10_22 by Sean Brennan, sbrennan@psu.edu
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
		narginchk(MAX_NARGIN-1,MAX_NARGIN);

		% Validate the initialXYPhi and finalXYPhi
		fcn_DebugTools_checkInputsToFunctions(initialXYPhi, '3column_of_numbers',[1 1]);
		fcn_DebugTools_checkInputsToFunctions(finalXYPhi, '3column_of_numbers',[1 1]);

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

initialSegment = fcn_INTERNAL_convertPointIntoUnitSegment(initialXYPhi);
finalSegment   = fcn_INTERNAL_convertPointIntoUnitSegment(finalXYPhi);

% Find where initialSegment intercepts the finalSegment
% Use Path library to find where this segment intersects the x-axis
% FORMAT:
%      [distance, location, path_segment, t, u] = ...
%         fcn_Path_findProjectionHitOntoPath(path,...
%         sensor_vector_start,sensor_vector_end,...
%         (flag_search_type),(figNum))
flagSearchType = 4; % any intersection of the sensor and path lines

[~, location, ~, tDistance, uDistance] = ...
	fcn_Path_findProjectionHitOntoPath(initialSegment,finalSegment(1,:),finalSegment(2,:),...
	flagSearchType,-1);


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
	clear vehicleParameters
	vehicleParameters.U = 2;  % U is forward velocity of vehicle in longitudinal direction, [m/s] (rule of thumb: 5 mph ~= 2* m/s)
	vehicleParameters.L = 2.700; % wheelbase in meters of a 2013 Prius, standard class
	% How to get the a and b values? (used for plotting)
	% Use the weight split!
	% FOR PRIUS
	% * Total Curb Weight: ~3,042 lbs (1,380 kg)
	% * Estimated Front Weight: ~1,855 lbs (61%)
	% * Estimated Rear Weight: ~1,187 lbs (39%)
	vehicleParameters.a =  vehicleParameters.L * 1187/3042;
	vehicleParameters.b =  vehicleParameters.L * 1855/3042;

	fcn_INTERNAL_plotUandT(vehicleParameters, ...
		initialXYPhi, finalXYPhi, ...
		uDistance, tDistance, location, figNum);
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
%% fcn_INTERNAL_convertPointIntoUnitSegment
function thisSegment = fcn_INTERNAL_convertPointIntoUnitSegment(thisXYphipoint)
thisXYpoint    = thisXYphipoint(1,1:2);
thisAngle      = thisXYphipoint(1,3);

% Define a vector based on the goal point phi angle
unitVector = [cos(thisAngle) sin(thisAngle)];
thisSegment = [thisXYpoint; thisXYpoint+unitVector];
end % Ends fcn_INTERNAL_convertPointIntoUnitSegment

%% fcn_INTERNAL_plotUandT
function fcn_INTERNAL_plotUandT(vehicleParameters, startPointXYphi, endPointXYphi, uDistance, tDistance, location, figNum)

startPhiInRadians = startPointXYphi(1,3);
unitStartVector = [cos(startPhiInRadians) sin(startPhiInRadians)];

endPhiInRadians = endPointXYphi(1,3);
unitEndVector = [cos(endPhiInRadians) sin(endPhiInRadians)];

% Plot scenario
h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, startPointXYphi', [0; 0], true, (figNum));
set(h_start{1},'DisplayName','Start');

h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, endPointXYphi', [0; 0], true, (figNum));
set(h_end{1},'DisplayName','Goal');

% Plot x-axis distance, t
quiver(startPointXYphi(1,1),startPointXYphi(1,2),unitStartVector(1)*tDistance, unitStartVector(2)*tDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
textLocation = [startPointXYphi(1,1), startPointXYphi(1,2)+0.5] + unitEndVector*uDistance/2;
text(textLocation(1,1), textLocation(1,2),sprintf('t = %.2f meters',tDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');

% Plot orientation projection onto X, u
quiver(endPointXYphi(1,1),endPointXYphi(1,2),unitEndVector(1)*uDistance, unitEndVector(2)*uDistance,0,'LineWidth',3,'ShowArrowHead','on','MaxHeadSize',0.3);
textLocation = [endPointXYphi(1,1), endPointXYphi(1,2)] + unitEndVector*uDistance/2;
text(textLocation(1,1), textLocation(1,2),sprintf('u = %.2f meters',uDistance),'FontSize',20,'BackgroundColor',[1 1 1],'HorizontalAlignment','center');


% Plot hit
plot(location(1,1), location(1,2),'.','MarkerSize',50);
end % Ends fcn_INTERNAL_plotUandT

