function [stateTrajectories, times, steeringAnglesUsed] = ...
    fcn_VD_forwardReachabilityTreeRK4(...
    initialStates, ...
    deltaT, timeInterval, steeringInterval, vehicleParameters, ...
    modelIDToUse, varargin)

%% fcn_VD_forwardReachabilityTreeRK4
%   Calculates the forward reachability tree of a vehicle using Runga Kutta
%   4th-order numerical solvers
%
% FORMAT:
%
%      [stateTrajectories, times, steeringAnglesUsed] =
%      fcn_VD_forwardReachabilityTreeRK4(initialStates, ...
%      deltaT, timeInterval, steeringInterval, vehicleParameters, ...
%      modelIDToUse, (figNum))
%
% INPUTS:
%
%      initialStates: A 1x3 vector of inital global pose in form of
%         [X Y Phi], which stand for:
% 
%         X: Global X position in meters
%
%         Y: Global Y position in meters
%
%         phi: Global yaw angle in radians, measured positive from X axis
%         to Y axis
%
%      deltaT: a 1x1 positive number denoting the time step to use, in
%      seconds
%
%      timeInterval: a 1x2 vector denoting [startTime endTime] in seconds
%
%      steeringInterval: an Mx1 vector denoting each steering angle to
%      evaluate in radians
%
%
%      vehicleParameters: a structure containing subfields of the following:
%
%          vehicleParameters.U: A 1x1 positive numeric value representing
%          the longitudinal velocity, in [m/s]
%
%      modelIDToUse: an 1x1 integer to indicate which model to use, from
%      one of the following:
%
%          0:  the kinematic point mass model is used by calling
%              fcn_VD_kinematicPointMassModelRK4
%
%          1:  the kinematic bicycle model is used by calling
%              fcn_VD_kinematicBicycleModelRK4

%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      stateTrajectories: An Nx3 vector of the state trajectories, with the
%      columns as [X Y Phi] in units of [m],[m],[rad]. If more than one
%      steering input is given, each trajectory is separated by NaN values.
%
%      times: An Nx1 vector of the simulation times, in seconds, with each
%      trajectory separated by NaN values
%
%      steeringAnglesUsed: An Nx1 vector of the steering values used in the
%      sim, in units of [rad], with each trajectory separated by NaN values
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_kinematicPointMassModelRK4
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_forwardReachabilityTreeRK4
%     for a full test suite.
%
% This function was written on 2026_09_16 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_forwardReachabilityTreeRK4
%
% 2026_09_16 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_forwardReachabilityTreeRK4
%   % * First write of function
%   % * Used fcn_VD_kinematicPointMassModelRK4 as starter


% TO-DO:
% - 2026_09_16 by Sean Brennan, sbrennan@psu.edu
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
        narginchk(MAX_NARGIN-1,MAX_NARGIN);

        % Validate that the initialStates input has 3 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(initialStates, '3column_of_numbers',[1 1]);

        % Validate that the deltaT input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(deltaT, '1column_of_numbers',[1 1]);

        % Validate that the timeInterval input has 2 columns, 1 row
        fcn_DebugTools_checkInputsToFunctions(timeInterval, '2column_of_numbers',[1 1]);

        % Validate that the steeringInterval input has 1 columns, 1+ rows
        fcn_DebugTools_checkInputsToFunctions(steeringInterval, '1column_of_numbers',[1 2]);

		% Check the parameters input is a structure
        assert(isstruct(vehicleParameters));

        % Validate that the modelIDToUse input has 1 column, 1 row of ints
        fcn_DebugTools_checkInputsToFunctions(modelIDToUse, '1column_of_integers',[1 1]);

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
%% RK4 in MATLAB Script
startTime = timeInterval(1);
endTime = timeInterval(2);
simulationTimes = (startTime:deltaT:endTime)';
N_timeSteps = length(simulationTimes); % This is the number of time steps we should have

if ~isequal(steeringInterval(1),steeringInterval(end))
    N_steeringInputs = length(steeringInterval);
else
    N_steeringInputs = 1;
end

% Fill names
switch modelIDToUse
    case 0
        simName = 'Kinematic Point Mass RK4';
    case 1
        simName = 'Kinematic Bicycle RK4';
    otherwise
        error('Unrecognized modelIDToUse input found. Function received entry of %.0d . Expecting values between 0 and 1.',modelIDToUse);
end

% Preallocate the output arrays
totalNumberOfNaNGapsRequired = N_steeringInputs-1;
stateTrajectories = nan(N_timeSteps*N_steeringInputs+totalNumberOfNaNGapsRequired,3);
times = nan(N_timeSteps*N_steeringInputs+totalNumberOfNaNGapsRequired,1);
steeringAnglesUsed = nan(N_timeSteps*N_steeringInputs+totalNumberOfNaNGapsRequired,1);

% Loop through steering inputs
offsetNaNEntries = 0; % How many NaN gaps have been inserted so far?

for ith_steeringInput = 1:N_steeringInputs
    thisRowOffset = N_timeSteps*(ith_steeringInput-1)+offsetNaNEntries;
    thisSteeringInput = steeringInterval(ith_steeringInput);

    % Initialize inputs for this situation
    inputsVsTime = [simulationTimes thisSteeringInput*ones(N_timeSteps,1)]; % [times steeringAngles]


    switch modelIDToUse
        case 0
            [thisSteeringStateTrajectories, thisSteeringTimes, thisSteeringAnglesUsed] = ...
                fcn_VD_kinematicPointMassModelRK4(initialStates, deltaT, ...
                timeInterval, inputsVsTime, vehicleParameters, (-1));
        case 1
            [thisSteeringStateTrajectories, thisSteeringTimes, thisSteeringAnglesUsed] = ...
                fcn_VD_kinematicBicycleModelRK4(initialStates, deltaT, ...
                timeInterval, inputsVsTime, vehicleParameters, (-1));
        otherwise
            error('Unrecognized modelIDToUse input found. Function received entry of %.0d . Expecting values between 0 and 1.',modelIDToUse);
    end


    % Save results
    rowRangeToFill = (thisRowOffset+1):(thisRowOffset+N_timeSteps);
    stateTrajectories(rowRangeToFill,:) = thisSteeringStateTrajectories;
    times(rowRangeToFill,:) = thisSteeringTimes;
    steeringAnglesUsed(rowRangeToFill,:) = thisSteeringAnglesUsed;

    % Increment the count of offsetNaNEntries?
    offsetNaNEntries = offsetNaNEntries+1;
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
    
    % plot the outputs
    h_plot = fcn_VD_plotTrajectory(stateTrajectories(:,1:2),(figNum));
	set(h_plot,'DisplayName',simName)
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

