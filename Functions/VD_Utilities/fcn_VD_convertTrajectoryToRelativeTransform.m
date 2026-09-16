function transformationMatrices = ...
    fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, varargin)

%% fcn_VD_convertTrajectoryToRelativeTransform
%   Calculates the sequence of homogenous transformation matrices that map
%   a state trajectory from an initial condition to the resulting
%   trajectory. The first row of the stateTrajectory input is treated as
%   the initial condition.
%
% FORMAT:
%
%      transformationMatrices = ...
%      fcn_VD_convertTrajectoryToRelativeTransform(stateTrajectory, (figNum))
%
% INPUTS:
%
%      stateTrajectory: A Nx3 vector of global poses in form of
%         [X Y Phi], which stand for:
% 
%         X: Global X position in meters
%
%         Y: Global Y position in meters
%
%         phi: Global yaw angle in radians, measured positive from X axis
%         to Y axis
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      transformationMatrices: An Nx1 cell array of 4x4 transformation
%      matrices in homogenous coordinates
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_convertTrajectoryToRelativeTransform
%     for a full test suite.
%
% This function was written on 2026_09_16 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_convertTrajectoryToRelativeTransform
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_convertTrajectoryToRelativeTransform
%   % * First write of function
%   % * Used fcn_VD_forwardReachabilityTreeRK4 as starter


% TO-DO:
% - 2026_09_16 by Sean Brennan, sbrennan@psu.edu
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

        % Validate that the stateTrajectory input has 3 columns, 2+ row
        fcn_DebugTools_checkInputsToFunctions(stateTrajectory, '3column_of_numbers',[2 3]);

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


%% Implements Transform calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the initial states
initialStates = stateTrajectory(1,:);

% How many time steps are involved?
N_trajectorySteps = length(stateTrajectory(:,1));

% Initialize the output cell array
transformationMatrices = cell(N_trajectorySteps,1);

% Loop through all the entries
for ith_step = 1:N_trajectorySteps
    deltaX = stateTrajectory(ith_step,1) - initialStates(1,1);
    deltaY = stateTrajectory(ith_step,2) - initialStates(1,2);
    deltaYaw = stateTrajectory(ith_step,3) - initialStates(1,3);
    rotations = [0 0 deltaYaw]; % radians 
    translations = [ deltaX deltaY 0]; % 1 2 3]; % radians
    transformationMatrix = fcn_VD_createTransformMatrix( rotations, translations, (-1));

    transformationMatrices{ith_step,1} = transformationMatrix;
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

    % Convert the transformationMatrices
    stackedMatrix = fcn_DebugTools_stackCellArrayIntoMatrix(transformationMatrices, (-1));
    stackedMatrixNoNans = stackedMatrix(~isnan(stackedMatrix(:,1)),:);

    % Fill in copies of initial position
    initialStatesHomogenous = [initialStates 1];
    
    % Calculate predicted positions
    predictedPositionsAllPoints = stackedMatrixNoNans*initialStatesHomogenous';
    predictedPositions_homogenousForm = (reshape(predictedPositionsAllPoints,4,[]))';
    predictedPositions = predictedPositions_homogenousForm(:,1:3);
    
    % plot the outputs
    h_plot = fcn_VD_plotTrajectory(stateTrajectory(:,1:2),(figNum));
	set(h_plot,'DisplayName','Original Input','LineWidth',5)
    h_plot = fcn_VD_plotTrajectory(predictedPositions(:,1:2),(figNum));
	set(h_plot,'DisplayName','Transformation','LineWidth',3)
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

