function transformationMatrix = fcn_VD_createTransformMatrix( rotations, translations, varargin)

%% fcn_VD_createTransformMatrix 
% creates a homogeneous transformation matrix
%
% FORMAT:
%
%      transformationMatrix = fcn_VD_createTransformMatrix( rotations, translations, (figNum))
%
% INPUTS:
%
%      rotations: A 1x3 vector of rotation values in the form 
%         [roll pitch yaw], in radians, where the terms stand for:
% 
%         roll: rotation about the x-axis
%
%         pitch: rotation about the y-axis
%
%         yaw: rotation about the z-axis
%
%      translations: a 1x3 vector containing translation distances in [x y
%      z] direction
%
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%   transformationMatrix: A 4x4 homogenous transformation matrix
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_createTransformMatrix
%     for a full test suite.
%
% This function was written on 2026_09_14 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% As: fcn_VD_createTransformMatrix
%
% 2026_09_14 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_createTransformMatrix
%   % * First write of function, 
%   % * using fcn_VD_derivativesKinematicPointMassModel as starter

% TO-DO:
% - 2026_09_14 by Sean Brennan, sbrennan@psu.edu
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

        % Check the rotations input to be sure it has 3 columns, 1 rows
        fcn_DebugTools_checkInputsToFunctions(rotations, '3column_of_numbers',[1 1]);

        % Check the translations input to be sure it has 3 columns and 1 row
        fcn_DebugTools_checkInputsToFunctions(translations, '3column_of_numbers',[1 1]);

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
% Fill in values for rotations
roll = rotations(1,1); % in rad
pitch = rotations(1,2); % in rad
yaw = rotations(1,3); % in rad

% Perform trig calculations
cosR = cos(roll);
sinR = sin(roll);
cosP = cos(pitch);
sinP = sin(pitch);
cosY = cos(yaw);
sinY = sin(yaw);

Rx = [
    1 0 0;
    0 cosR -sinR;
    0 sinR  cosR;
    ];

Ry = [
    cosP 0 sinP;
    0    1    0;
    -sinP 0 cosP;
    ];

Rz = [
    cosY -sinY 0;
    sinY  cosY 0;
    0     0    1;
    ];
Rotations = Rz*Ry*Rx;

transformationMatrix = [...
    Rotations translations';
    zeros(1,3) 1;
    ];

% For debugging
% Check that this is correct?
if 1==1
    Ttemp = tform(se3([yaw, pitch, roll],"eul",'ZYX',translations));

    assert(isequal(round(transformationMatrix,4),round(Ttemp,4)));

    % Create the translation matrix
    translation_matrix = makehgtform('translate',translations);
    % Create the z-rotate matrix
    rotation_matrix_z = makehgtform('zrotate',yaw);
    % Create the y-rotate matrix
    rotation_matrix_y = makehgtform('yrotate',pitch);
    % Create the x-rotate matrix
    rotation_matrix_x = makehgtform('xrotate',roll);
    % Compute the transformation matrix
    Transformation_Matrix = translation_matrix*rotation_matrix_z*...
        rotation_matrix_y*rotation_matrix_x;

    assert(isequal(round(transformationMatrix,4),round(Transformation_Matrix,4)));

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
    
    % plot the outputs as a transformation on a unit cube

    % Clear/prepare figure
    hold on;
    grid on;
    figure(figNum);

    % Define vertices for a unit-length 3D cube centered at the origin
    % Length = 1, so coordinates span from -0.5 to 0.5
    [X, Y, Z] = meshgrid([-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5]);

    % Plot the cube vertices as points
    scatter3(X(:), Y(:), Z(:), 60, 'filled', 'MarkerFaceColor', 'b');

    % Draw the edges of the cube for better visualization
    edges = [
        -0.5 -0.5 -0.5;  0.5 -0.5 -0.5;
        0.5 -0.5 -0.5;  0.5  0.5 -0.5;
        0.5  0.5 -0.5; -0.5  0.5 -0.5;
        -0.5  0.5 -0.5; -0.5 -0.5 -0.5;
        -0.5 -0.5  0.5;  0.5 -0.5  0.5;
        0.5 -0.5  0.5;  0.5  0.5  0.5;
        0.5  0.5  0.5; -0.5  0.5  0.5;
        -0.5  0.5  0.5; -0.5 -0.5  0.5;
        -0.5 -0.5 -0.5; -0.5 -0.5  0.5;
        0.5 -0.5 -0.5;  0.5 -0.5  0.5;
        0.5  0.5 -0.5;  0.5  0.5  0.5;
        -0.5  0.5 -0.5; -0.5  0.5  0.5
        ];
    for i = 1:2:size(edges, 1)
        plot3(edges(i:i+1, 1), edges(i:i+1, 2), edges(i:i+1, 3), 'k-', 'LineWidth', 1);
    end

    % Define starting points (centers) of the positive-facing faces
    % Positive X face center: (0.5, 0, 0)
    % Positive Y face center: (0, 0.5, 0)
    % Positive Z face center: (0, 0, 0.5)
    face_centers = [
        0.5, 0.0, 0.0; % +X face
        0.0, 0.5, 0.0; % +Y face
        0.0, 0.0, 0.5  % +Z face
        ];

    % Define the unit vectors pointing outwards from these faces
    vectors = [
        1, 0, 0; % Points in +X direction
        0, 1, 0; % Points in +Y direction
        0, 0, 1  % Points in +Z direction
        ];
	endPointsOfVectors = face_centers+vectors;

    % Plot the unit vectors using quiver3
    % 'AutoScale', 'off' ensures the vectors keep their literal unit length of 1
    quiver3(face_centers(:,1), face_centers(:,2), face_centers(:,3), ...
        vectors(:,1), vectors(:,2), vectors(:,3), ...
        0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Formatting the plot
    xlabel('X Axis');
    ylabel('Y Axis');
    zlabel('Z Axis');
    title('Unit Cube with Positive-Facing Normal Vectors');
    axis equal;
    view(3); % Set to default 3D view orientation


	% Add the rotation
	vertices = [X(:) Y(:) Z(:) ones(length(X(:)),1)];
	rotated_vertices = (transformationMatrix*vertices')';

	homogenous_edges = [edges ones(length(edges(:,1)),1)];
	rotated_edges = (transformationMatrix*homogenous_edges')';

	homogenous_face_centers = [face_centers ones(length(face_centers(:,1)),1)];
	rotated_face_centers = (transformationMatrix*homogenous_face_centers')';

	homogenous_endPointsOfVectors = [endPointsOfVectors ones(length(endPointsOfVectors(:,1)),1)];
	rotated_endPointsOfVectors = (transformationMatrix*homogenous_endPointsOfVectors')';
	rotated_vectors = rotated_endPointsOfVectors - rotated_face_centers;

	% Plot the cube vertices as points
    scatter3(rotated_vertices(:,1), rotated_vertices(:,2), rotated_vertices(:,3), 60, 'filled', 'MarkerFaceColor', 'g');

	for i = 1:2:size(edges, 1)
		plot3(rotated_edges(i:i+1, 1), rotated_edges(i:i+1, 2), rotated_edges(i:i+1, 3), '-', 'LineWidth', 1,'Color',0.3*[1 1 1]);
	end

    % Plot the unit vectors using quiver3
    % 'AutoScale', 'off' ensures the vectors keep their literal unit length of 1
    quiver3(rotated_face_centers(:,1), rotated_face_centers(:,2), rotated_face_centers(:,3), ...
        rotated_vectors(:,1), rotated_vectors(:,2), rotated_vectors(:,3), ...
        0, 'Color', 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

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

