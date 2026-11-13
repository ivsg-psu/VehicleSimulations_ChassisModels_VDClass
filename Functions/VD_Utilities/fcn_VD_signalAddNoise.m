function signalWithNoise = fcn_VD_signalAddNoise(signal, varargin)

%% fcn_VD_signalAddNoise
%   Adds noise to a signal
%
% FORMAT:
%
%      signalWithNoise = fcn_VD_signalAddNoise(signal, (noiseTypes), (noiseOptions), (figNum))
%
% INPUTS:
%
%      signal: an Nx1 input
%
%      (OPTIONAL INPUTS)
%
%      noiseTypes: an array of integers indicating which noise types to
%      add. Options include:
%
%           1: constant offset (default)
%
%           2: random normal noise
%
%           3: mean-returning random walk via auto-regression
%
%      noiseOptions: a cell array of length equal to the maximum noiseType,
%      allowing user to select options. Default is 1 for all values.
%
%           noiseOptions{1}(1): constant offset value
%
%           noiseOptions{2}(1): random normal noise standard deviation
%
%           noiseOptions{3}(1): rho, forgetting factor for random walk
%           (must be less than 1 for mean-returning) 
% 
%           noiseOptions{3}(2): standard deviation in integrated noise
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      signal: an Nx1 input

%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_signalAddNoise
%     for a full test suite.
%
% This function was written on 2026_11_12
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_signalAddNoise
%
% 2026_11_13 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_signalAddNoise
%   % * First write of function,
%   % * Using fcn_VD_derivatives2DOFBicycleModel as starter

% TO-DO:
% - 2026_11_12 by Sean Brennan, sbrennan@psu.edu
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
		narginchk(0,MAX_NARGIN);

		% Check the signal input to be sure it has 1 column, 2+ rows
		fcn_DebugTools_checkInputsToFunctions(signal, '1column_of_numbers',[2 3]);
		 
		% % Check the inputFrontRoadWheelAngleRadians input to be sure it has 1 column and 1 row
		% fcn_DebugTools_checkInputsToFunctions(inputFrontRoadWheelAngleRadians, '1column_of_numbers',[1 1]);
		% 
		% % Check the vehicleParameters input is a structure
		% assert(isstruct(vehicleParameters));
		% 
		% % Check the U input to be sure it has 1 column and 1 row, positive
		% fcn_DebugTools_checkInputsToFunctions(U, 'positive_1column_of_numbers',[1 1]);
	end
end


% Does the user want to specify the noiseTypes input?
noiseTypes = 1; % Default case
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)        
        noiseTypes = temp; % Set the value from user input
    end
end


% Does the user want to specify the noiseOptions input?
noiseOptions = cell(5,1); % Default case
if 3 <= nargin
    temp = varargin{2};
    if ~isempty(temp)        
        noiseOptions = temp; % Set the value from user input
    end
end


% % Does the user want to specify the simulationTimes input?
% simulationTimes = (0:0.01:4.5)';
% if 3 <= nargin
%     temp = varargin{3};
%     if ~isempty(temp)        
%         simulationTimes = temp; % Set the value from user input
%     end
% end

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

% Initialize output
signalWithNoise = signal;
Npoints = length(signal);

additionalNoise = 0*signal;

thisIndex = 1; % constant offset
if any(noiseTypes==thisIndex) 		
	if ~isempty(noiseOptions{thisIndex})
		noiseOffset = noiseOptions{thisIndex}(1); % Degrees of steering amplitude for input sinewave
	else
		noiseOffset = 1; % default
	end

	additionalNoise = additionalNoise + noiseOffset*ones(Npoints,1);
end

thisIndex = 2; % random normal
if any(noiseTypes==thisIndex) 		
	if ~isempty(noiseOptions{thisIndex})
		noiseAmplitude = noiseOptions{thisIndex}(1); % Degrees of steering amplitude for input sinewave
	else
		noiseAmplitude = 1; % default
	end

	additionalNoise = additionalNoise + noiseAmplitude*randn(Npoints,1);
end

thisIndex = 3; % mean-returning random walk via auto-regression
if any(noiseTypes==thisIndex) 		
	if ~isempty(noiseOptions{thisIndex})
		rho = noiseOptions{3}(1); % forgetting factor 
		sigma = noiseOptions{3}(2); % standard deviation in integrated noise
	else
		rho = 0.95;
		sigma = 1;
	end
	mu = 0; % Mean
	thisNoise = zeros(Npoints,1);
	thisNoise(1) = sigma*randn;
	for ith_point = 2:Npoints
		thisNoise(ith_point) = mu + rho*(thisNoise(ith_point-1)-mu) + sigma*randn;
	end


	additionalNoise = additionalNoise + thisNoise;
end




signalWithNoise =  signalWithNoise + additionalNoise; 


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
	hold on;

	plot(signalWithNoise,'-','DisplayName','With noise added');
	plot(signal,'-','DisplayName','Input','LineWidth',3);
	legend('Interpreter','none');
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

