function [tDistancesEachAngle, uDistancesEachAngle] = fcn_VD_parkingSolveUTFromStandardTurn(...
	standardParkingTurnPath, varargin)

%% fcn_VD_parkingSolveUTFromStandardTurn
%   Given a standard parking turn profile, starting from [0 0 0] in [X Y
%   phi] coordinates in [m m rad] respectively, finds the x-axis intercepts
%   for each location on the turn profile (t) and the distance from the profile
%   to that x-axis intercept (u). These are used in the parking algorithm
%   to solve for parking manuevers given arbitrary parking turn profiles
%
% FORMAT:
%
%      [tDistancesEachAngle, uDistancesEachAngle] = fcn_VD_parkingSolveUTFromStandardTurn(...
%    	standardParkingTurnPath, (figNum))
%
% INPUTS:
%
%      standardParkingTurnPath: an Nx3 vector of the [X, Y, yaw] positions
%      of the reference path used for turning. The reference path is
%      assumed to start at [0 0 0], with units of [m m rad].
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
%      tDistance: a Nx1 scalar of the distances where the  angle of each
%      standardParkingTurnPath would intersect the local x-axis of the
%      standardParkingTurnPath, as measured in the direction of
%      the x-axis
%
%      uDistance: a Nx1 scalar of the distances where the  angle of each
%      standardParkingTurnPath would intersect the local x-axis of the
%      standardParkingTurnPath, as measured in the direction of
%      the sampled point
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_Path_findProjectionHitOntoPath
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingSolveUTFromStandardTurn
%     for a full test suite.
%
% This function was written on 2026_10_28 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingSolveUTFromStandardTurn
%
% 2026_10_28 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingSolveUTFromStandardTurn
%   % * Used fcn_VD_parking+PlotSolution as starter

% TO-DO:
% - 2026_10_22 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
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

		% Validate the reference path and maneuver parameters
		fcn_DebugTools_checkInputsToFunctions(standardParkingTurnPath, '3columns_of_numbers',[2 3]);

	end
end

% % Does user want to specify the anglesToPlotTangents input?
% anglesToPlotTangents = []; % Default empty
% if 6 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         anglesToPlotTangents = temp;
%     end
% end

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

Nsteps = length(standardParkingTurnPath(:,1));

% Initialize outputs
tDistancesEachAngle = nan(Nsteps,1);
uDistancesEachAngle = nan(Nsteps,1);

standardParkingStartPoint = standardParkingTurnPath(1,1:3); % For plotting
for ith_step = 2:1:Nsteps
	thisXYphipoint = standardParkingTurnPath(ith_step,1:3);

	flag_doAnimationPlot = (1==0 && (ith_step==2 || (mod(ith_step,10)==0)));

	if flag_doAnimationPlot
		tempFigNum = 234;
		figure(tempFigNum);clf;
	else
		tempFigNum = -1;
	end

	[tDistance, uDistance] = fcn_VD_parkingXYPhiToTU(standardParkingStartPoint(1,:), thisXYphipoint, (tempFigNum));

	tDistancesEachAngle(ith_step,1) = tDistance;
	uDistancesEachAngle(ith_step,1) = uDistance;

	if flag_doAnimationPlot
		axis equal
		axis([-5 20 -5 20]);
		title(sprintf('Parking Angle: %.0f deg',thisXYphipoint(1,3)*180/pi))

		pause(0.01);
		exportgraphics(gcf,"parkingXIntercepts.gif","Append",ith_step>2);
	end
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


	% Plot results
	figure(figNum);

	subplot(1,3,1);
	% Call the function just to plot
	plot(standardParkingTurnPath(:,1),standardParkingTurnPath(:,2),'-','LineWidth',5,'DisplayName','Standard Parking Maneuver');

	h_start = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(1,:)', [0; 0], true, (figNum));
	set(h_start{1},'DisplayName','Start');

	h_end = fcn_VD_plotStickFigure2DOFBicycleModel(vehicleParameters, standardParkingTurnPath(end,:)', [0; 0], true, (figNum));
	set(h_end{1},'DisplayName','Goal');
	legend('Interpreter','none');

	subplot(1,3,2);
	plot(standardParkingTurnPath(:,3)*180/pi, tDistancesEachAngle,'.-');
	xlabel('Angle [deg]');
	ylabel('t-distances [m]');

	subplot(1,3,3);
	plot(standardParkingTurnPath(:,3)*180/pi, uDistancesEachAngle,'.-');
	hold on;
	plot(standardParkingTurnPath(:,3)*180/pi, -1*tDistancesEachAngle,'.-');
	xlabel('Angle [deg]');
	ylabel('u-distances [m]');


	sgtitle('Calculating x-axis projection intercepts (t) and distances (u) for the standard turn', 'Interpreter','none');

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

