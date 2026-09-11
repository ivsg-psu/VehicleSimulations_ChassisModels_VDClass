function dXdt = fcn_VD_body2GlobalCoordinatesDerivatives( X, xdot, varargin)

%% fcn_VD_body2GlobalCoordinatesDerivatives
%   This function calculates velocites in global coordinates.
%
% FORMAT:
%
%      dXdt = fcn_VD_body2GlobalCoordinatesDerivatives( X, xdot, (figNum))
%
% INPUTS:
%
%      X: A 3x1 vector of global pose in the form of 
%         [X; Y; Phi], which stand for:
% 
%         X: Global X position in meters
%
%         Y: Global Y position in meters
%
%         phi: Global yaw angle in radians, measured positive from X axis
%         to Y axis
%
%      xdot: a 3x1 vector of the body-fixed velocities in the form of
%         [U; V; r], which stand for:
%
%         U: Longitudinal velocity [m/s]
%
%         V: Lateral velocity [m/s]
%
%         r: yawrate [rad/s]
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%   dXdt: A 3x1 vector of velocities in global coordinates
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_derivativesKinematicPointMassModel
%     for a full test suite.
%
% This function was written on 2026_01_26 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_Body2GlobalCoordinates
% 
% 2021_05_16 by Satya Prasad, szm888@psu.edu
% - In fcn_VD_Body2GlobalCoordinates
%   % * First write of function
%
% As: fcn_VD_Body2GlobalCoordinates
%
% 2026_09_10 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_body2GlobalCoordinatesDerivatives
%   % * Renamed function to more "standard" form
%   % * Created function from fcn_VD_Body2GlobalCoordinates in "old" folder
%   % * Changed number of input arguments


% TO-DO:
% - 2026_01_26 by Sean Brennan, sbrennan@psu.edu
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

        % Check the X input to be sure it has 1 column, 3 rows
        fcn_DebugTools_checkInputsToFunctions(X, '1column_of_numbers',[3 3]);

        % Check the xdot input to be sure it has 1 column, 3 rows
        fcn_DebugTools_checkInputsToFunctions(xdot, '1column_of_numbers',[3 3]);

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

% Fill in variables
psi = X(3);

U = xdot(1);
V = xdot(2);
r = xdot(3);

% Planar rotation
dXdt   = U*cos(psi)-V*sin(psi);
dYdt   = U*sin(psi)+V*cos(psi);
dPhidt = r;

% Fill in derivative vector output
dXdt   = [dXdt; dYdt; dPhidt];


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

	quiver(X(1),X(2),dXdt(1), dXdt(2), 0);
	hold on;
	axis equal;
	axis padded;

	% Plot the rotation. this is done by putting a small, grey vector at
	% the end that points in the direction and magnitude of the rotation.
	rotationAngle = dPhidt*0.01; % Have to guess a time step
	changeVector = [dXdt(1) dXdt(2)];
	newChange = changeVector*[cos(rotationAngle) sin(rotationAngle); -sin(rotationAngle) cos(rotationAngle)];
	rotationStartPoint = [X(1),X(2)]+changeVector;
	rotationEndPoint = [X(1),X(2)]+newChange;
	differenceVector = rotationEndPoint - rotationStartPoint;
	quiver(rotationStartPoint(1),rotationStartPoint(2), differenceVector(1), differenceVector(2),0,'Color',0.8*[1 1 1]);


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
