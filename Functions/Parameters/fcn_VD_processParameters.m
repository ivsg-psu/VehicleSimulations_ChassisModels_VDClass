function vehicleParametersProcessed = fcn_VD_processParameters(vehicleParameters, varargin)
%% fcn_VD_processParameters
% Performs checking on vehicle parameters to calculate, if necessary,
% missing terms. The calculations include:
%
%      
%
% FORMAT:
%
%      vehicleParametersProcessed = fcn_VD_processParameters(vehicleParameters, (figNum))
%
% INPUTS:
%
%      vehicleParameters: a structure containing the vehicle parameters.
%
%      (OPTIONAL INPUTS)
% 
%      figNum: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. 
%
% OUTPUTS:
%

%
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_processParameters
%     for a full test suite.
%
% This function was written on 2026_11_07 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% 2026_11_07 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_processParameters
%   % * First write of the function
%   % * Used fcn_VD_getStudentVehicle as starter


% TO-DO:
% 
% 2026_11_07 by Sean Brennan, sbrennan@psu.edu
% - fill in to-do items here.

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 2; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %     % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_DEBUGTOOLS_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_DEBUGTOOLS_FLAG_CHECK_INPUTS");
    MATLABFLAG_DEBUGTOOLS_FLAG_DO_DEBUG = getenv("MATLABFLAG_DEBUGTOOLS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_DEBUGTOOLS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_DEBUGTOOLS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_DEBUGTOOLS_FLAG_DO_DEBUG); 
        flag_check_inputs  = str2double(MATLABFLAG_DEBUGTOOLS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_figNum = 999978; %#ok<NASGU>
else
    debug_figNum = []; %#ok<NASGU>
end
%% check input arguments
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

if 0 == flag_max_speed
    if flag_check_inputs == 1
        % Are there the right number of inputs?
        narginchk(MAX_NARGIN-2,MAX_NARGIN);

    end
end

% % Does user want to specify the rangeToPickFrom input?
% vehicleParametersFile = fullfile(pwd,'Data','vehicleParameters.mat'); % Default
% if 1 <= nargin
% 	temp = varargin{1};
% 	if ~isempty(temp)
% 		vehicleParametersFile = temp;
% 		if 0 == flag_max_speed
% 			if flag_check_inputs == 1
% 				% Check the vehicleParametersFile input, make sure it is characters
% 				fcn_DebugTools_checkInputsToFunctions(vehicleParametersFile,'_of_char_strings');
% 			end
% 		end
% 	end
% end

% % Does user want to specify the numToPick input?
% numToPick = length(fromStructure); % Default is to select them all
% if 4 <= nargin
%     temp = varargin{2};
%     if ~isempty(temp)
%         numToPick = temp;
%     end
% end
% 
% % Does user want to specify the flagRandomizeSelections input?
% flagRandomizeSelections = false; % Default is to not randomize
% if 4 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         flagRandomizeSelections = temp;
%     end
% end

% Check to see if user specifies figNum?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp)
        figNum = temp; %#ok<NASGU>
        flag_do_plots = 1;
    end
end


%% Main code starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


vehicleParametersTemplate = fcn_INTERNAL_fillTemplate;

% Initialize variable
vehicleParameters = vehicleParametersTemplate;





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

    % Nothing to do here


end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
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

function selections = fcn_INTERNAL_queryField(selections, thisName, thisField,vehicleParameters,numQuestions)

textField = sprintf('%s_Text',thisField);
formatField = sprintf('%s_PrintFormat',thisField);
rangeField = sprintf('%s_AllowableRange',thisField);

defaultAnswer = '-missing-';
if ~isempty(vehicleParameters.(thisField))
	commandToEval = sprintf('sprintf(''%s'',vehicleParameters.%s)',vehicleParameters.(formatField),thisField);
	defaultAnswer = eval(commandToEval);
end

% Standard items
selections(numQuestions).MenuChar = num2str(numQuestions);
selections(numQuestions).Name = thisField;
selections(numQuestions).Text = sprintf(') %s: What is %s?\n',thisField, vehicleParameters.(textField));
selections(numQuestions).AnswerDefault = defaultAnswer;
selections(numQuestions).FunctionMore = sprintf('fcn_help_%s(%.0f, answers, selections);',thisName,numQuestions);
selections(numQuestions).FunctionPreMenu = {[]};
selections(numQuestions).FunctionMoreInputs = {30};
selections(numQuestions).FunctionSubmission = '[answers, numBadOptionInputs, flag_exitMain] = fcn_INTERNAL_enterData(answers, selections, selectedOptionCharacters, numBadOptionInputs)';
selections(numQuestions).FunctionSubmissionOptions = {'.'};
selections(numQuestions).isAllowableMenuOption = true;
selections(numQuestions).AnswerGradingCorrect = nan;
selections(numQuestions).AnswerGradingPoints = 1;
selections(numQuestions).AnswerGradingType = 'not empty';
selections(numQuestions).AnswerGradingOptions = {[]};
selections(numQuestions).AnswerAllowableRange = vehicleParameters.(rangeField);
selections(numQuestions).AnswerPrintFormat = vehicleParameters.(formatField);

if contains(vehicleParameters.(formatField),'d')
	% Integer
	selections(numQuestions).AnswerType = '1column_of_integers';
	selections(numQuestions).AnswerConversionFunction = 'str2double';
	selections(numQuestions).AnswerTypeOptions = [1 1];		
	selections(numQuestions).AnswerConversionFunction = 'str2double';

elseif contains(vehicleParameters.(formatField),'s')
	% String
	selections(numQuestions).AnswerType = '_of_char_strings';
	selections(numQuestions).AnswerConversionFunction = '';
	selections(numQuestions).AnswerTypeOptions = [];
	selections(numQuestions).AnswerConversionFunction = [];

elseif contains(vehicleParameters.(formatField),'f')
	% Float
	selections(numQuestions).AnswerType = '1column_of_numbers';
	selections(numQuestions).AnswerConversionFunction = 'str2double';
	selections(numQuestions).AnswerTypeOptions = [1 1];
	selections(numQuestions).AnswerConversionFunction = 'str2double';

end



end % Ends fcn_INTERNAL_queryField

%% fcn_INTERNAL_fillTemplate
function vehicleParameters = fcn_INTERNAL_fillTemplate
vehicleParameters = struct;

% Fill in descriptors
vehicleParameters.Year = [];
vehicleParameters.Year_Text = 'the model year of the vehicle (example: 2013)';
vehicleParameters.Year_PrintFormat = '%.0d';
vehicleParameters.Year_AllowableRange = [1800 2050];
vehicleParameters.Year_Units = [];

vehicleParameters.Make = [];
vehicleParameters.Make_Text = 'the make of the vehicle - the company name or brand (e.g., Toyota)';
vehicleParameters.Make_PrintFormat = '%s';
vehicleParameters.Make_AllowableRange = [];
vehicleParameters.Make_Units = [];

vehicleParameters.Model = [];
vehicleParameters.Model_Text = 'the model of the vehicle - the specific product line or name given to the vehicle by the manufacturer (e.g., Camry, Prius, or RAV4)';
vehicleParameters.Model_PrintFormat = '%s';
vehicleParameters.Model_AllowableRange = [];
vehicleParameters.Model_Units = [];

vehicleParameters.Trim = [];
vehicleParameters.Trim_Text = 'the trim line of the vehicle - the specific package level or feature set of the model (e.g., LE, XLE, or Limited)';
vehicleParameters.Trim_PrintFormat = '%s';
vehicleParameters.Trim_AllowableRange = [];
vehicleParameters.Trim_Units = [];

% Masses
vehicleParameters.Mass= [];
vehicleParameters.Mass_Text = 'the test mass of the vehicle - the mass being tested (should be between curb and gross masses), in kg';
vehicleParameters.Mass_PrintFormat = '%.2f';
vehicleParameters.Mass_AllowableRange = [100 9000];
vehicleParameters.Mass_Units = 'kg';

vehicleParameters.MassCurb= [];
vehicleParameters.MassCurb_Text = 'the base curb mass of the vehicle - the empty weight, in kg';
vehicleParameters.MassCurb_PrintFormat = '%.2f';
vehicleParameters.MassCurb_AllowableRange = [100 9000];
vehicleParameters.MassCurb_Units = 'kg';

vehicleParameters.MassGross= [];
vehicleParameters.MassGross_Text = 'the gross mass of the vehicle weight - the maximum fully loaded mass including passengers and cargo, in kg';
vehicleParameters.MassGross_PrintFormat = '%.2f';
vehicleParameters.MassGross_AllowableRange = [100 9000];
vehicleParameters.MassGross_Units = 'kg';

vehicleParameters.MassCurbFrontAxle = [];
vehicleParameters.MassCurbFrontAxle_Text = 'the front-axle curb mass of the vehicle, in kg';
vehicleParameters.MassCurbFrontAxle_PrintFormat = '%.2f';
vehicleParameters.MassCurbFrontAxle_AllowableRange = [100 9000];
vehicleParameters.MassCurbFrontAxle_Units = 'kg';

vehicleParameters.MassCurbRearAxle = [];
vehicleParameters.MassCurbRearAxle_Text = 'the rear-axle curb mass of the vehicle, in kg';
vehicleParameters.MassCurbRearAxle_PrintFormat = '%.2f';
vehicleParameters.MassCurbRearAxle_AllowableRange = [100 9000];
vehicleParameters.MassCurbRearAxle_Units = 'kg';

% Inertias
vehicleParameters.Ixx= [];
vehicleParameters.Ixx_Text = 'the test x-axis moment of inertia of the vehicle (roll), in kg-m^2';
vehicleParameters.Ixx_PrintFormat = '%.2f';
vehicleParameters.Ixx_AllowableRange = [1 90000];
vehicleParameters.Ixx_Units = 'kg-m^2';

vehicleParameters.IxxCurb= [];
vehicleParameters.IxxCurb_Text = 'the curb x-axis moment of inertia of the vehicle (roll)  - the minimum value, in kg-m^2';
vehicleParameters.IxxCurb_PrintFormat = '%.2f';
vehicleParameters.IxxCurb_AllowableRange = [1 90000];
vehicleParameters.IxxCurb_Units = 'kg-m^2';

vehicleParameters.IxxGross= [];
vehicleParameters.IxxGross_Text = 'the gross x-axis moment of inertia of the vehicle (roll) - the maximum value, in kg-m^2';
vehicleParameters.IxxGross_PrintFormat = '%.2f';
vehicleParameters.IxxGross_AllowableRange = [1 90000];
vehicleParameters.IxxGross_Units = 'kg-m^2';

vehicleParameters.Iyy= [];
vehicleParameters.Iyy_Text = 'the test y-axis moment of inertia of the vehicle (pitch), in kg-m^2';
vehicleParameters.Iyy_PrintFormat = '%.2f';
vehicleParameters.Iyy_AllowableRange = [1 90000];
vehicleParameters.Iyy_Units = 'kg-m^2';

vehicleParameters.IyyCurb= [];
vehicleParameters.IyyCurb_Text = 'the curb y-axis moment of inertia of the vehicle (pitch)  - the minimum value, in kg-m^2';
vehicleParameters.IyyCurb_PrintFormat = '%.2f';
vehicleParameters.IyyCurb_AllowableRange = [1 90000];
vehicleParameters.IyyCurb_Units = 'kg-m^2';

vehicleParameters.IyyGross= [];
vehicleParameters.IyyGross_Text = 'the gross y-axis moment of inertia of the vehicle (pitch) - the maximum value, in kg-m^2';
vehicleParameters.IyyGross_PrintFormat = '%.2f';
vehicleParameters.IyyGross_AllowableRange = [1 90000];
vehicleParameters.IyyGross_Units = 'kg-m^2';

vehicleParameters.Izz= [];
vehicleParameters.Izz_Text = 'the test x-axis moment of inertia of the vehicle (yaw), in kg-m^2';
vehicleParameters.Izz_PrintFormat = '%.2f';
vehicleParameters.Izz_AllowableRange = [1 90000];
vehicleParameters.Izz_Units = 'kg-m^2';

vehicleParameters.IzzCurb= [];
vehicleParameters.IzzCurb_Text = 'the curb x-axis moment of inertia of the vehicle (yaw)  - the minimum value, in kg-m^2';
vehicleParameters.IzzCurb_PrintFormat = '%.2f';
vehicleParameters.IzzCurb_AllowableRange = [1 90000];
vehicleParameters.IzzCurb_Units = 'kg-m^2';

vehicleParameters.IzzGross= [];
vehicleParameters.IzzGross_Text = 'the gross x-axis moment of inertia of the vehicle (yaw) - the maximum value, in kg-m^2';
vehicleParameters.IzzGross_PrintFormat = '%.2f';
vehicleParameters.IzzGross_AllowableRange = [1 90000];
vehicleParameters.IzzGross_Units = 'kg-m^2';

% Distances
vehicleParameters.DistancesWheelbase = [];
vehicleParameters.DistancesWheelbase_Text = 'the wheel base distance - the distance from the center of the front tires to the center of the rear tires - of the vehicle, in m';
vehicleParameters.DistancesWheelbase_PrintFormat = '%.3f';
vehicleParameters.DistancesWheelbase_AllowableRange = [0.01 90];
vehicleParameters.DistancesWheelbase_Units = 'm';

vehicleParameters.DistancesLength = [];
vehicleParameters.DistancesLength_Text = 'the vehicle length - the distance from the rearmost part of the rear bumper to frontmost part of the front bumper - of the vehicle, in m';
vehicleParameters.DistancesLength_PrintFormat = '%.3f';
vehicleParameters.DistancesLength_AllowableRange = [0.01 90];
vehicleParameters.DistancesLength_Units = 'm';

vehicleParameters.DistancesFrontOverhang = [];
vehicleParameters.DistancesFrontOverhang_Text = 'the front overhang distance - the distance from the center of the rear axle to the front-most location of the front bumper, in m';
vehicleParameters.DistancesFrontOverhang_PrintFormat = '%.3f';
vehicleParameters.DistancesFrontOverhang_AllowableRange = [0.01 90];
vehicleParameters.DistancesFrontOverhang_Units = 'm';

vehicleParameters.DistancesRearOverhang = [];
vehicleParameters.DistancesRearOverhang_Text = 'the rear overhang distance - the distance from the center of the rear axle to the rear-most location of the rear bumper, in m';
vehicleParameters.DistancesRearOverhang_PrintFormat = '%.3f';
vehicleParameters.DistancesRearOverhang_AllowableRange = [0.01 90];
vehicleParameters.DistancesRearOverhang_Units = 'm';

vehicleParameters.DistancesFrontTrack = [];
vehicleParameters.DistancesFrontTrack_Text = 'the front track distance - the distance from the center of the front right tires to the center of the front left tires - of the vehicle, in m';
vehicleParameters.DistancesFrontTrack_PrintFormat = '%.3f';
vehicleParameters.DistancesFrontTrack_AllowableRange = [0.01 90];
vehicleParameters.DistancesFrontTrack_Units = 'm';

vehicleParameters.DistancesRearTrack = [];
vehicleParameters.DistancesRearTrack_Text = 'the rear track distance - the distance from the center of the rear right tires to the center of the rear left tires - of the vehicle, in m';
vehicleParameters.DistancesRearTrack_PrintFormat = '%.3f';
vehicleParameters.DistancesRearTrack_AllowableRange = [0.01 90];
vehicleParameters.DistancesRearTrack_Units = 'm';

vehicleParameters.DistancesWidth = [];
vehicleParameters.DistancesWidth_Text = 'the vehicle width - the distance from the rightmost to leftmost component, usually excluding the mirrors - of the vehicle, in m';
vehicleParameters.DistancesWidth_PrintFormat = '%.3f';
vehicleParameters.DistancesWidth_AllowableRange = [0.01 90];
vehicleParameters.DistancesWidth_Units = 'm';

vehicleParameters.DistancesWidthWithMirrors = [];
vehicleParameters.DistancesWidthWithMirrors_Text = 'the vehicle width with mirrors - the distance from the rightmost to leftmost component, INCLUDING the mirrors - of the vehicle, in m';
vehicleParameters.DistancesWidthWithMirrors_PrintFormat = '%.3f';
vehicleParameters.DistancesWidthWithMirrors_AllowableRange = [0.01 90];
vehicleParameters.DistancesWidthWithMirrors_Units = 'm';

vehicleParameters.DistancesHeight = [];
vehicleParameters.DistancesHeight_Text = 'the vehicle height - the distance from the ground to the top-most part of the vehicle, excluding any antennas, in m';
vehicleParameters.DistancesHeight_PrintFormat = '%.3f';
vehicleParameters.DistancesHeight_AllowableRange = [0.01 90];
vehicleParameters.DistancesHeight_Units = 'm';

vehicleParameters.DistancesGroundClearance = [];
vehicleParameters.DistancesGroundClearance_Text = 'the ground clearance - the distance from the ground to the bottom-most part of the vehicle, excluding the tires, in m';
vehicleParameters.DistancesGroundClearance_PrintFormat = '%.3f';
vehicleParameters.DistancesGroundClearance_AllowableRange = [0.01 90];
vehicleParameters.DistancesGroundClearance_Units = 'm';

vehicleParameters.DistancesCGFromRearAxle = [];
vehicleParameters.DistancesCGFromRearAxle_Text = 'the CG distance forward from the center of the rear axle of the vehicle, in m';
vehicleParameters.DistancesCGFromRearAxle_PrintFormat = '%.3f';
vehicleParameters.DistancesCGFromRearAxle_AllowableRange = [0.01 90];
vehicleParameters.DistancesCGFromRearAxle_Units = 'm';

vehicleParameters.DistancesCGAboveGround = [];
vehicleParameters.DistancesCGAboveGround_Text = 'the CG distance of the vehicle above the ground, in m';
vehicleParameters.DistancesCGAboveGround_PrintFormat = '%.3f';
vehicleParameters.DistancesCGAboveGround_AllowableRange = [0.01 90];
vehicleParameters.DistancesCGAboveGround_Units = 'm';

% Tires
vehicleParameters.TireSidewallMarkings = [];
vehicleParameters.TireSidewallMarkings_Text = 'the tire sidewall markings for the vehicle (example: 195/65R15)';
vehicleParameters.TireSidewallMarkings_PrintFormat = '%s';
vehicleParameters.TireSidewallMarkings_AllowableRange = [];
vehicleParameters.TireSidewallMarkings_Units = [];

vehicleParameters.TireFrontCorneringCoefficient = [];
vehicleParameters.TireFrontCorneringCoefficient_Text = 'the front tire cornering coefficient - the cornering stiffness divided by the normal force, in units of 1/rad (typical values are 10 rad^-1)';
vehicleParameters.TireFrontCorneringCoefficient_PrintFormat = '%.3f';
vehicleParameters.TireFrontCorneringCoefficient_AllowableRange = [0.01 90];
vehicleParameters.TireFrontCorneringCoefficient_Units = '1/rad';

vehicleParameters.TireRearCorneringCoefficient = [];
vehicleParameters.TireRearCorneringCoefficient_Text = 'the rear tire cornering coefficient - the cornering stiffness divided by the normal force, in units of 1/rad (typical values are 14 rad^-1)';
vehicleParameters.TireRearCorneringCoefficient_PrintFormat = '%.3f';
vehicleParameters.TireRearCorneringCoefficient_AllowableRange = [0.01 90];
vehicleParameters.TireRearCorneringCoefficient_Units = '1/rad';

% Turn Radius
vehicleParameters.MinimumTurnRadius = [];
vehicleParameters.MinimumTurnRadius_Text = 'the minimum turn radius of the vehicle, in meters';
vehicleParameters.MinimumTurnRadius_PrintFormat = '%.2f';
vehicleParameters.MinimumTurnRadius_AllowableRange = [1 100];
vehicleParameters.MinimumTurnRadius_Units = 'm';

% Engine
vehicleParameters.EngineDisplacement = [];
vehicleParameters.EngineDisplacement_Text = 'the engine displacement, in liters';
vehicleParameters.EngineDisplacement_PrintFormat = '%.2f';
vehicleParameters.EngineDisplacement_AllowableRange = [1 100];
vehicleParameters.EngineDisplacement_Units = 'L';

vehicleParameters.EngineCylinders = [];
vehicleParameters.EngineCylinders_Text = 'the number of cylinders';
vehicleParameters.EngineCylinders_PrintFormat = '%.0d';
vehicleParameters.EngineCylinders_AllowableRange = [1 50];
vehicleParameters.EngineCylinders_Units = [];

vehicleParameters.EnginePeakHorsepower = [];
vehicleParameters.EnginePeakHorsepower_Text = 'the peak power of the engine, in horsepower';
vehicleParameters.EnginePeakHorsepower_PrintFormat = '%.2f';
vehicleParameters.EnginePeakHorsepower_AllowableRange = [1 5000];
vehicleParameters.EnginePeakHorsepower_Units = 'hp';

vehicleParameters.EngineSpeedAtPeakHorsepower = [];
vehicleParameters.EngineSpeedAtPeakHorsepower_Text = 'the engine speed at peak power, in RPM';
vehicleParameters.EngineSpeedAtPeakHorsepower_PrintFormat = '%.2f';
vehicleParameters.EngineSpeedAtPeakHorsepower_AllowableRange = [1 100000];
vehicleParameters.EngineSpeedAtPeakHorsepower_Units = 'rpm';

vehicleParameters.EngineTorqueAtPeakHorsepower = [];
vehicleParameters.EngineTorqueAtPeakHorsepower_Text = 'the engine torque at peak power, in ft-lb';
vehicleParameters.EngineTorqueAtPeakHorsepower_PrintFormat = '%.2f';
vehicleParameters.EngineTorqueAtPeakHorsepower_AllowableRange = [1 10000];
vehicleParameters.EngineTorqueAtPeakHorsepower_Units = 'ft-lb';

vehicleParameters.AccelerationZeroToSixtyTime = [];
vehicleParameters.AccelerationZeroToSixtyTime_Text = 'the 0 to 60 mph fastest time, in seconds';
vehicleParameters.AccelerationZeroToSixtyTime_PrintFormat = '%.2f';
vehicleParameters.AccelerationZeroToSixtyTime_AllowableRange = [0 1000];
vehicleParameters.AccelerationZeroToSixtyTime_Units = 'sec';
end % Ends fcn_INTERNAL_fillTemplate