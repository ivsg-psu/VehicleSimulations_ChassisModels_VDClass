function transformationMatrix = fcn_VD_fitTransformationSVD( sourcePoints, targetPoints, varargin)

%% fcn_VD_fitTransformationSVD 
% computes the best-fitting rigid transformation between two point sets using SVDx
%
% FORMAT:
%
%      transformationMatrix = fcn_VD_fitTransformationSVD(inputPoints, targetPoints, (weightingArray), (figNum))
%
% INPUTS:
%
%      sourcePoints: Nx3 array Source point cloud coordinates in LiDAR
%      frame.
%
%      targetPoints: Nx3 array Corresponding reference point coordinates
%      (e.g., in GPS or ENU frame).
%
%      (OPTIONAL INPUTS)
%
%      weightingArray: Nx1 array Weighting coefficients for each point pair. If
%      empty, uniform weights assumed.
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%   transformationMatrix: A 4x4 homogenous transformation matrix from
%   inputPoints to targetPoints
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_fitTransformationSVD
%     for a full test suite.
%
% This function was written on 2026_09_15 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu
% Function transcribed from Xinyu Cao, xfc5113@psu.edu, created 2024-01-24

% REVISION HISTORY:
%
% As: fcn_VD_fitTransformationSVD
%
% 2026_09_15 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_fitTransformationSVD
%   % * First write of function, 
%   % * using fcn_LiDARPoseEstimation_FitTransformationSVD as starter
%   % * See https://github.com/ivsg-psu/Publications_Journals_2025_JAVS_Cao_ExtrinsicCalibration

% TO-DO:
% - 2026_09_15 by Sean Brennan, sbrennan@psu.edu
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
        narginchk(MAX_NARGIN-2,MAX_NARGIN);

        % Check the inputPoints input
        fcn_DebugTools_checkInputsToFunctions(...
            sourcePoints, '3column_of_numbers');
        fcn_DebugTools_checkInputsToFunctions(...
            targetPoints, '3column_of_numbers');

    end
end

% Does user want to specify the weightingArray?
weightingArray = ones(size(sourcePoints, 1), 1);
if 3 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        weightingArray = temp;
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


%% SVD transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the weighted centroids of both point sets

source_center  = sum(weightingArray .* sourcePoints, 1) / sum(weightingArray);
target_center = sum(weightingArray .* targetPoints, 1) / sum(weightingArray);

% Center coordinates
source_centered  = sourcePoints - source_center;
target_centered = targetPoints - target_center;
W_diag = diag(weightingArray);
% Compute weighted cross-covariance
H = source_centered' * W_diag * target_centered;

% SVD for optimal rotation
[U, ~, V] = svd(H);
R_rotation = V * U';

% Ensure proper right-handed rotation
if det(R_rotation) < 0
    F = eye(3);
    F(3,3) = -1;
    R_rotation = V * F * U';
end

% Translation vector
translation_vector = target_center' - R_rotation * source_center';

% Construct transformation
transformationMatrix = tform(se3(R_rotation, translation_vector'));

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

	fprintf('Rotation Matrix:\n'); disp(R_rotation);
	fprintf('Translation Vector:\n'); disp(translation_vector');
	fprintf('Homogeneous Transformation:\n'); disp(transformationMatrix);

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

