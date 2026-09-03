function [stateTrajectory, t, inputHistory] = fcn_VD_FirstOrderEngineModelRK4(initialStates, deltaT, timeInterval, inputsVsTime, parameters, varargin)

%% fcn_VD_FirstOrderEngineModelRK4
%   Simulates the the first order engine model using Runga Kutta 4th-order
%
% FORMAT:
%
%      [stateTrajectory, t, inputHistory] =
%      fcn_VD_FirstOrderEngineModelRK4(initialStates, deltaT,
%      timeInterval, inputsVsTime, parameters, (figNum))
%
% INPUTS:
%
%      initialStates: A 1x1 vector of inital global pose in form of
%         [Omega], which stand for:
% 
%         Omega: engine rotational speed in rad/sec
%
%      deltaT: a 1x1 positive number denoting the time step to use, in
%      seconds
%
%      timeInterval: a 1x2 vector denoting [startTime endTime] in seconds
%
%      inputsVsTime: a Mx2 vector denoting 
%      [inputTime inputValues] 
%      where the list shows the time and input torque side-by-side
%      in units of [sec N-m] respectively. This is interpolated using
%      linear interpolation at the sampling times. For times outside the
%      given interval, zero values are used.
%
%      parameters: a structure containing subfields of the following
%          parameters.J: A 1x1 positive numeric value representing the engine
%          rotational inertia, in [kg-m^2]
%
%          parameters.B: A 1x1 positive numeric value representing the
%          viscous drag, in [(N-m)/(rad/sec)]
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      stateTrajectory: An Nx1 vector of the state trajectory, with the
%      columns as [Omega] in units of [rad/sec]
%
%      t: An Nx1 vector of the simulation times, in seconds
%
%      inputHistory: An Nx1 vector of the input values used in the sim,
%      in units of [N-m]
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_kinematicBicycleModel
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_FirstOrderEngineModelRK4
%     for a full test suite.
%
% This function was written on 2026_09_02 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_FirstOrderEngineModelRK4
%
% 2026_09_02 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_FirstOrderEngineModelRK4
%   % * First write of fcn_VD_FirstOrderEngineModelRK4 function
%   % * Used fcn_VD_kinematicBicycleModelRK4 as a starter

% TO-DO:
% - 2026_09_02 by Sean Brennan, sbrennan@psu.edu
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

        % Validate that the initialStates input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(initialStates, '1column_of_numbers',[1 1]);

        % Validate that the deltaT input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(deltaT, '1column_of_numbers',[1 1]);

        % Validate that the timeInterval input has 2 columns, 1 row
        fcn_DebugTools_checkInputsToFunctions(timeInterval, '2column_of_numbers',[1 1]);

        % Validate that the inputsVsTime input has 2 columns, 2+
		% rows
        fcn_DebugTools_checkInputsToFunctions(inputsVsTime, '2column_of_numbers',[2 3]);

		% Check the parameters input is a structure
        assert(isstruct(parameters));

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
N_variables = length(initialStates);

% Initialize variables
stateTrajectory = nan(N_timeSteps,N_variables);
t = nan(N_timeSteps,1);
inputHistory = interp1(inputsVsTime(:,1), inputsVsTime(:,2), simulationTimes,'linear',0);


% Set initial conditions
currentStates = initialStates;

% Fill in the parameters
J = parameters.J;
B = parameters.B;

for ith_time = 1:N_timeSteps
	thisTime       = simulationTimes(ith_time);

	% Fill in the results to save
	t(ith_time,1)  = thisTime; % Update time
	stateTrajectory(ith_time,:)   = currentStates;


	% Use Runga-Kutta to predict next position
	thisState = currentStates'; %#ok<NASGU>
	thisInput = inputHistory(ith_time,1);
	[~, y] = fcn_VD_RungeKutta(...
		@(t,y) fcn_VD_derivativesFirstOrderEngineModel(y,thisInput, J, B, -1), ...
		currentStates', thisTime, deltaT, -1);

	currentStates = y';
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
    % h_plot = fcn_VD_plotTrajectory(stateTrajectory(:,1:2),(figNum));
	% set(h_plot,'DisplayName','XY Trajectory (MATLAB RK4)')

    figure(figNum);
    hold on;
    subplot(2,1,1)
    h_plot = plot(t,inputHistory(:,1),'-');
	set(h_plot,'DisplayName','Input (MATLAB RK4)')

    subplot(2,1,2)
    h_plot = plot(t,stateTrajectory(:,1),'-');
	set(h_plot,'DisplayName','Engine Rotational Speed (MATLAB RK4)')
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

