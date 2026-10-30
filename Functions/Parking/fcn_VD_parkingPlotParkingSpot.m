function [h_fig, xyPlottingPoints] = fcn_VD_parkingPlotParkingSpot(vehicleParameters, vehicleGlobalStates, varargin)

%% fcn_VD_parkingPlotParkingSpot
%   This function plots the 2DOF bicycle model as a stick figure
%
% FORMAT:
%
%      [h_fig, xyPlottingPoints] = fcn_VD_parkingPlotParkingSpot(vehicleParameters, vehicleGlobalStates, (flagPlotFromRear), (figNum))
%
% INPUTS:
%
%      vehicleParameters: a structure with the following subfields
%          % Example below is approximately a Ford Taurus
%            vehicleParameters(1).U          = 20;     % x-axis velocity, in m/s
%            vehicleParameters(1).m          = 1031.9; % mass in kg
%            vehicleParameters(1).Iz         = 1850;   % z-axis rotational inertia in kg-m^2
%            vehicleParameters(1).a          = 0.9271; % distance from front axle to CG, in meters
%            vehicleParameters(1).b          = 1.5621; % distance from rear axle to CG, in meters
%            vehicleParameters(1).Caf        = 77500;  % Cornering stiffness of front axle, N/rad;
%            vehicleParameters(1).Car	     = 116250; % Cornering stiffness of rear axle, N/rad;
%
%      vehicleGlobalStates: a 1x3 vector denoting the [X, Y, yaw] positions
%      respectively, in units of [m m rad]
%
%      (OPTIONAL INPUTS)
%
%      flagPlotFromRear: set equal to true to plot vehicle from
%      measurements at rear axle, rather than at CG. Default is "false".
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      h_fig: a cell array of the figure handles for the plots
%
%      xyPlottingPoints: an Nx2 matrix of the [x y] points being plotted
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_parkingPlotParkingSpot
%     for a full test suite.
%
% This function was written on 2026_10_30 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_parkingPlotParkingSpot
%
% 2026_10_30 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingPlotParkingSpot
%   % * Used fcn_VD_plotSt+ickFigure2DOFBicycleModel as starter
%
% 2026_10_30b by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_parkingPlotParkingSpot
%   % * Used fcn_VD_plotSt+ickFigure2DOFBicycleModel as starter
%   % * Added xyPlottingPoints as an output
%   % * Fixed XYphi to be a row vector

% TO-DO:
% - 2026_10_30 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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

		% Check the vehicleParameters input is a structure
        assert(isstruct(vehicleParameters));

		% Validate that the vehicleGlobalStates input has 1 row 3 columns
        fcn_DebugTools_checkInputsToFunctions(vehicleGlobalStates, '3column_of_numbers',[1 1]);

    end
end

% Does user want to specify the flagPlotFromRear input?
flagPlotFromRear = false; % Default is plot from the CG
if 4 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flagPlotFromRear = temp;
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
% Nothing to do in main

h_fig = cell(1,1);

if ~isfield('L',vehicleParameters)
	vehicleParameters.L = vehicleParameters.a + vehicleParameters.b;
end

X = vehicleGlobalStates(1);
Y = vehicleGlobalStates(2);
yawAngle = vehicleGlobalStates(3);

% Calculate the vehicle longitudinal axis
angleRad = yawAngle;
frontGapDistance = 0.5; % Meters, this is the distance from the front axle to the stripe
rotation = [cos(angleRad) sin(angleRad); -sin(angleRad) cos(angleRad)];
if flagPlotFromRear
	offsetVehicle = [0 0; vehicleParameters.L 0]*rotation + [X*ones(2,1) Y*ones(2,1)] ;
	distToFrontAxle = vehicleParameters.L + frontGapDistance;
else
	offsetVehicle = [-vehicleParameters.b 0; vehicleParameters.a 0;]*rotation + [X*ones(2,1) Y*ones(2,1)] ;
	distToFrontAxle = vehicleParameters.a + frontGapDistance;

end

% A standard commercial parking space in the United States is 9 feet
% (108 inches) wide and 18 feet long. There's 3.281 feet in a meter
parkingWidth = 9/3.218;
parkingLength = 18/3.281;

unitXvector = (offsetVehicle(2,:)-offsetVehicle(1,:))./vehicleParameters.L;
unitYvector = unitXvector*[0 1; -1 0];


mainAxis = ones(2,1)*[X Y] + [distToFrontAxle-parkingLength; distToFrontAxle]*unitXvector;
leftStripe  = mainAxis + ones(2,1)*(parkingWidth/2)*unitYvector;
rightStripe = mainAxis + ones(2,1)*(-1*parkingWidth/2)*unitYvector;
topStripe = [rightStripe(2,:); leftStripe(2,:)];

xyPlottingPoints = [offsetVehicle; nan nan; rightStripe; nan nan; leftStripe; nan nan; topStripe];

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


	

	% Plot the vehicle as a stick figure
	h_fig = plot(xyPlottingPoints(:,1), xyPlottingPoints(:,2),'-','Color',0.1*[1 1 1], 'LineWidth',5,'DisplayName','Parking Spot');

	axis equal;
	grid on;

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


