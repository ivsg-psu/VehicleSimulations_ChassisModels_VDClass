function T_transformation = fcn_LiDARPoseEstimation_FitTransformationSVD(sourcePoints, targetPoints, varargin)
% fcn_LiDARPoseEstimation_FitTransformationSVD
% Computes the best-fitting rigid transformation between two point sets using SVD.
%
% FORMAT:
%   T_transformation = fcn_LiDARPoseEstimation_FitTransformationSVD(inputPoints, targetPoints, W_array)
%
% INPUTS:
%   inputPoints: Nx3 array
%       Source point cloud coordinates in LiDAR frame.
%
%   targetPoints: Nx3 array
%       Corresponding reference point coordinates (e.g., in GPS or ENU frame).
%
%   W_array: Nx1 array
%       Weighting coefficients for each point pair. If empty, uniform weights assumed.
%
% OUTPUTS:
%   T_transformation: SE(3) object
%       Homogeneous transformation from inputPoints to targetPoints.
%
% DEPENDENCIES:
%   - fcn_DebugTools_checkInputsToFunctions
%   - se3 class
%
% EXAMPLE USAGE:
%   See: script_test_fcn_LiDARPoseEstimation_FitTransformationSVD
%
% AUTHOR:
%   Xinyu Cao, xfc5113@psu.edu
%   Created: 2024-01-24

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

flag_do_debug = false;
flag_check_inputs = true;

if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the inputPoints input
        fcn_DebugTools_checkInputsToFunctions(...
            sourcePoints, '3column_of_numbers');
        fcn_DebugTools_checkInputsToFunctions(...
            targetPoints, '3column_of_numbers');
       

end

W_array = ones(size(sourcePoints, 1), 1);
if nargin >= 3
    W_array = varargin{1};
end


%% Solve for the circle
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

source_center  = sum(W_array .* sourcePoints, 1) / sum(W_array);
target_center = sum(W_array .* targetPoints, 1) / sum(W_array);

% Center coordinates
source_centered  = sourcePoints - source_center;
target_centered = targetPoints - target_center;
W_diag = diag(W_array);
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
T_transformation = se3(R_rotation, translation_vector');

%% Debug Output
if flag_do_debug
    fprintf('Rotation Matrix:\n'); disp(R_rotation);
    fprintf('Translation Vector:\n'); disp(translation_vector');
    fprintf('Homogeneous Transformation:\n'); disp(T_transformation.tform);
end

end
