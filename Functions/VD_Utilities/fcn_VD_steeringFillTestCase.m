function steeringInputsVsTime = fcn_VD_steeringFillTestCase(varargin)

%% fcn_VD_steeringFillTestCase
%   Fill in test steering cases for the front steering using typical
%   maneuvers
%
% FORMAT:
%
%      steeringInputsVsTime = fcn_VD_steeringFillTestCase( (testCaseID), (testCaseOptions), (simulationTimes), (figNum))
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%
%      testCaseID: An integer selecting the test case to use, from the
%      following:
%
%         0: (default) a constant steering input
%         1: a single wave of a sine wave input
%         2: a pseudo-random binary sequence
%         3: a continuously repeating sine wave input
%
%      testCaseOptions: a cell array allowing user to specify parameters
%      for each of the test cases, from the following:
%
%         0: (default) a constant steering input
%
%				testCaseOptions{1}: steeringAmplitudeDegrees, degrees of
%				steering amplitude for constant input (default is 5 deg)
%
%         1: a single wave of a sine wave input
%
%				testCaseOptions{1}: steeringAmplitudeDegrees, degrees of
%				steering amplitude for input sinewave (default is 5 deg)
%
%				testCaseOptions{2}: periodOfSineWave, the period of the
%				sine wave, units are seconds. A typical lane change is
%				about 3 to 4 seconds based on experimental highway
%				measurements. Default is 3.
%
%				testCaseOptions{3}: start time of the wave. Default is
%				after 1 second.
%
%         2: a pseudo-random binary sequence
%
%				testCaseOptions{1}: maxSteeringAmplitudeDegrees - the
%				maximum degrees of steering amplitude for input sinewave.
%				Default is 5 degrees.
%
%				testCaseOptions{2}: minSteeringAmplitudeDegrees - the
%				minimum degrees of steering amplitude for input sinewave.
%				Default is -5 degrees.
%
%				testCaseOptions{3}: minDwellTime - the minimum time within
%				a steering state, units are seconds. Default is 0.1
%				seconds.
%
%				testCaseOptions{4}: maxDwellTime - the maximum time within
%				a steering state, units are seconds. Default is 1
%				seconds.
%
%         3: (default) a sine wave input
%
%				testCaseOptions{1}: steeringAmplitudeDegrees, degrees of
%				steering amplitude for input sinewave (default is 5 deg)
%
%				testCaseOptions{2}: periodOfSineWave, the period of the
%				sine wave, units are seconds. A typical lane change is
%				about 3 to 4 seconds based on experimental highway
%				measurements. Default is 3.
%
%      simulationTimes: a Nx1 input allowing user to specify the time
%      vector to use. Default is:
%         simulationTimes = (0:0.01:4.5)';
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%   steeringInputsVsTime: A Nx2 vector of steering inputs to use for
%   testing. The format is:
%      [timeInSeconds  steeringInputsInRadians]
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_steeringFillTestCase
%     for a full test suite.
%
% This function was written on 2026_11_12
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_steeringFillTestCase
%
% 2026_11_12 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_steeringFillTestCase
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

		% % Check the X input to be sure it has 1 column, 3 rows
		% fcn_DebugTools_checkInputsToFunctions(X, '5column_of_numbers',[1 1]);
		% 
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


% Does the user want to specify the testCaseID input?
testCaseID = 0; % Default case
if 1 <= nargin
    temp = varargin{1};
    if ~isempty(temp)        
        testCaseID = temp; % Set the value from user input
    end
end


% Does the user want to specify the testCaseOptions input?
testCaseOptions = {}; % Default case
if 2 <= nargin
    temp = varargin{2};
    if ~isempty(temp)        
        testCaseOptions = temp; % Set the value from user input
    end
end


% Does the user want to specify the simulationTimes input?
simulationTimes = (0:0.01:4.5)';
if 3 <= nargin
    temp = varargin{3};
    if ~isempty(temp)        
        simulationTimes = temp; % Set the value from user input
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

% Initialize output
steeringInputsVsTime = [simulationTimes 0*simulationTimes]; % [times steering angles]

switch testCaseID
	case 0
		modeName = 'Constant steering input';
		if ~isempty(testCaseOptions)
			if length(testCaseOptions)~=1
				error('Incorrect number of test case options given. Expected 1 but found: %.0d',length(testCaseOptions));
			else
				steeringAmplitudeDegrees = testCaseOptions{1}; % Degrees of steering amplitude for input sinewave
			end
		else
			steeringAmplitudeDegrees = 5; % 2 degrees of steering amplitude for input sinewave
		end

		additionalSteeringInputsInRadians = steeringAmplitudeDegrees*pi/180*ones(length(simulationTimes),1);
	case 1
		modeName = 'Single wave of a sinewave steering input';
		if ~isempty(testCaseOptions)
			if length(testCaseOptions)~=3
				error('Incorrect number of test case options given. Expected 3 but found: %.0d',length(testCaseOptions));
			else
				steeringAmplitudeDegrees = testCaseOptions{1}; % Degrees of steering amplitude for input sinewave
				periodOfSineWave = testCaseOptions{2}; % Period, units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
				startTimeOfSineWave = testCaseOptions{3}; % When sine wave starts, units are seconds
			end
		else
			steeringAmplitudeDegrees = 5; % 2 degrees of steering amplitude for input sinewave
			periodOfSineWave = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
			startTimeOfSineWave = 1; % When sine wave starts, units are seconds
		end

		offsetTime = (simulationTimes - startTimeOfSineWave);
		indicesToFill = (offsetTime>0) & (offsetTime<periodOfSineWave);
		additionalSteeringInputsInRadians = steeringInputsVsTime(:,2)*0;
		additionalSteeringInputsInRadians(indicesToFill,1) = steeringAmplitudeDegrees*pi/180*sin((2*pi/periodOfSineWave)*offsetTime(indicesToFill));
	case 2
		modeName = 'Steering input is a pseudo-random binary sequence';
		if ~isempty(testCaseOptions)
			if length(testCaseOptions)~=4
				error('Incorrect number of test case options given. Expected 3 but found: %.0d',length(testCaseOptions));
			else
				maxSteeringAmplitudeDegrees = testCaseOptions{1}; % Degrees of steering amplitude for input sinewave
				minSteeringAmplitudeDegrees = testCaseOptions{2}; % Degrees of steering amplitude for input sinewave
				minDwellTime = testCaseOptions{3}; % Minimum time within a steering state, units are seconds
				maxDwellTime = testCaseOptions{4}; % When sine wave starts, units are seconds
			end
		else
			maxSteeringAmplitudeDegrees = 5; % Degrees of steering amplitude for input sinewave
			minSteeringAmplitudeDegrees = -5; % Degrees of steering amplitude for input sinewave
			minDwellTime = 0.1; % Minimum time within a steering state, units are seconds
			maxDwellTime = 1;   % Maximum time within a steering state, units are seconds
		end
				
		NsequencesNeededMax = length(simulationTimes);

		dwellTimeRange = maxDwellTime-minDwellTime;
		if minDwellTime<=0
			error('Min dwell time must be greater than zero');
		end
		if dwellTimeRange<0
			error('Max dwell time cannot be less than minimum dwell time');
		end

		dwellTimes = minDwellTime + rand(NsequencesNeededMax,1)*dwellTimeRange;
		dwellTimeChanges = cumsum(dwellTimes);
		dwellTimeChangesInTimeLimit = [0; dwellTimeChanges(dwellTimeChanges<max(simulationTimes)); max(simulationTimes)];
		Ndwells = length(dwellTimeChangesInTimeLimit);

		amplitudeRange = maxSteeringAmplitudeDegrees-minSteeringAmplitudeDegrees;
		if amplitudeRange<0
			error('maxSteeringAmplitudeDegrees must be greater than or equal to minSteeringAmplitudeDegrees');
		end
		
		steeringForEachDwellInRadians = (minSteeringAmplitudeDegrees + rand(Ndwells,1)*amplitudeRange)*pi/180;

		% Initialize values
		additionalSteeringInputsInRadians = steeringInputsVsTime(:,2)*0;
		
		% Loop through dwells, filling in steering for each
		for ith_dwell = 2:Ndwells
			dwellStartTime = dwellTimeChangesInTimeLimit(ith_dwell-1);
			dwellEndTime  =  dwellTimeChangesInTimeLimit(ith_dwell);
			theseIndices = (simulationTimes>=dwellStartTime) & (simulationTimes<dwellEndTime);
			additionalSteeringInputsInRadians(theseIndices,1) = steeringForEachDwellInRadians(ith_dwell);
		end
	case 3
		modeName = 'Continuous sinewave steering input';
		if ~isempty(testCaseOptions)
			if length(testCaseOptions)~=2
				error('Incorrect number of test case options given. Expected 2 but found: %.0d',length(testCaseOptions));
			else
				steeringAmplitudeDegrees = testCaseOptions{1}; % Degrees of steering amplitude for input sinewave
				periodOfSineWave = testCaseOptions{2}; % Period, units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
			end
		else
			steeringAmplitudeDegrees = 5; % 2 degrees of steering amplitude for input sinewave
			periodOfSineWave = 3; % Units are seconds. A typical lane change is about 3 to 4 seconds based on experimental highway measurements
		end

		additionalSteeringInputsInRadians = steeringAmplitudeDegrees*pi/180*sin((2*pi/periodOfSineWave)*simulationTimes);
	otherwise
		error('Unknown testCaseID encountered with value of: %.0d',testCaseID);
end
steeringInputsVsTime(:,2) =  steeringInputsVsTime(:,2) + additionalSteeringInputsInRadians; % [times steering angles]


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

	[xb, yb] = stairs(steeringInputsVsTime(:,1),steeringInputsVsTime(:,2)*180/pi);

	plot( xb, yb,'-','DisplayName',modeName);
	legend('Interpreter','none');
	xlabel('Time [sec]','Interpreter','none');
	ylabel('Front steering input [deg]','Interpreter','none');
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

