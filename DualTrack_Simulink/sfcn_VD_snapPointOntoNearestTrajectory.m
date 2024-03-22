%%%%%%%%%%%%%%%%% sfcn_VD_snapPointOntoNearestTrajectory %%%%%%%%%%%%%%%%%%
% Purpose:
%   Finds the pose on the interpolated trajectory that is nearest to the 
%   input vehicle position
% 
%   Note: Change the sample time depending on the simulation.
%
% INPUTS:
%   INPUT 1: 'ref_trajectory' is the complete vehicle_trajectory
%   INPUT 2: 'vehicle_cg_ahead' is the cg of the vehicle in EN
% 
% OUTPUTS:
%   OUTPUT 1: Trajectory of the vehicle nearest to INPUT 2
% 
% PARAMETERS:
%   Parameter 1: Structure mapping column names of INPUT 1 with column
%   numbers
%   Parameter 2: Parameter limiting the search range for nearest neighbour
%   Parameter 3: Size of INPUT 1
% 
% DEPENDENCIES:
%   fcn_Path_snapPointOntoNearestTraversal
%
% Author: Satya Prasad (szm888@psu.edu)
% Created: 2021/07/20
% ======== to do list ============
% 1. Method to compute 'length_of_search_window' need to be updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_VD_snapPointOntoNearestTrajectory(block)
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Register the parameters
block.NumDialogPrms     = 3;
block.DialogPrmsTunable = {'Nontunable', 'Nontunable', 'Nontunable'};

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = block.DialogPrm(3).Data;
block.InputPort(2).Dimensions = 2;

% Override output port properties
block.OutputPort(1).Dimensions = block.DialogPrm(3).Data(2);

% Set block sample time to be inherited
block.SampleTimes = [0.01, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitConditions);
block.RegBlockMethod('Outputs', @Output); % Required

function DoPostPropSetup(block)
%% Setup Dwork
block.NumDworks = 1;
block.Dwork(1).Name = 'flag'; % variable to reset pertinent variables for every new simulation
block.Dwork(1).Dimensions = 1;
block.Dwork(1).DatatypeID = 0;
block.Dwork(1).Complexity = 'Real';

function InitConditions(block)
%% Initialize Dwork
block.Dwork(1).Data = 1; % acts as a flag 

function Output(block)
%% Output function with search window length depending on reference trajectory
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
persistent ref_trajectory % reference trajectory with unique station coordinates
persistent length_of_ref_trajectory % length of reference trajectory
persistent start_index_for_local_search % beginning index for the search window
persistent length_of_search_window % length of the search window

%% Read INPUTS
if block.Dwork(1).Data
    % Happens only at the beginning of a new simulation
    ref_trajectory = block.InputPort(1).Data; % reference trajectory
end
vehicle_cg = [block.InputPort(2).Data(1), block.InputPort(2).Data(2)]; % CG of the vehicle in EN

%% Read parameters
fields_trajectory = block.DialogPrm(1).Data; % MATLAB structure mapping column names to column numbers

if block.Dwork(1).Data
    % Happens only at the beginning of a new simulation
    %% Initialization
    length_of_ref_trajectory = size(ref_trajectory,1); % length of reference trajectory
    
    % variable to limit the search range
    length_of_search_window = max(2,ceil(block.DialogPrm(2).Data/min(diff(ref_trajectory(:,fields_trajectory.station)))));
    
    % Create a traversal using reference trajectory
    reference_traversal.X       = ref_trajectory(:,fields_trajectory.east);
    reference_traversal.Y       = ref_trajectory(:,fields_trajectory.north);
    reference_traversal.Station = ref_trajectory(:,fields_trajectory.station);
    reference_traversal.Yaw     = ref_trajectory(:,fields_trajectory.yaw);
    
    [closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,~,~] = ...
        fcn_Path_snapPointOntoNearestTraversal(vehicle_cg,reference_traversal); % call to the snap function
    
    output_index = first_path_point_index; % index that decides ouput properties other than pose, velocity, and station
    start_index_for_local_search = first_path_point_index; % shift start_index depending on the nearest neighbour
    vehicle_trajectory = ref_trajectory(output_index,:); % Initialize the Output-1
    % update the output with estimated pose and station
    vehicle_trajectory(fields_trajectory.east)    = closest_path_point(1);
    vehicle_trajectory(fields_trajectory.north)   = closest_path_point(2);
    vehicle_trajectory(fields_trajectory.yaw)     = path_point_yaw;
    vehicle_trajectory(fields_trajectory.station) = s_coordinate;
    
    block.Dwork(1).Data = 0;
    
else
    % if start_index_for_local_search reaches maximum, then set it to one less than the maximum
    if start_index_for_local_search == length_of_ref_trajectory
        start_index_for_local_search = length_of_ref_trajectory-1;
    end
    % set the upper bound to limit the search range depending on size of 
    % INPUT-1/Parameter-1 and length_of_search_window
    end_index_for_local_search = min(start_index_for_local_search+length_of_search_window,length_of_ref_trajectory);
    
    % Arrange the trajectory in traversal structure
    reference_traversal.X       = ref_trajectory(start_index_for_local_search:end_index_for_local_search,...
        fields_trajectory.east);
    reference_traversal.Y       = ref_trajectory(start_index_for_local_search:end_index_for_local_search,...
        fields_trajectory.north);
    reference_traversal.Station = ref_trajectory(start_index_for_local_search:end_index_for_local_search,...
        fields_trajectory.station);
    reference_traversal.Yaw     = ref_trajectory(start_index_for_local_search:end_index_for_local_search,...
        fields_trajectory.yaw);
    
    [closest_path_point,s_coordinate,path_point_yaw,first_path_point_index,~,~] = ...
        fcn_Path_snapPointOntoNearestTraversal(vehicle_cg, reference_traversal); % call to the snap function
    
    output_index = start_index_for_local_search+first_path_point_index-1; % index that decides ouput properties other than pose, velocity, and station
    start_index_for_local_search = output_index; % shift start_index depending on the nearest neighbour
    vehicle_trajectory = ref_trajectory(output_index,:); % Initialize the Output-1
    % update the output with estimated pose and station
    vehicle_trajectory(fields_trajectory.east)    = closest_path_point(1);
    vehicle_trajectory(fields_trajectory.north)   = closest_path_point(2);
    vehicle_trajectory(fields_trajectory.yaw)     = path_point_yaw;
    vehicle_trajectory(fields_trajectory.station) = s_coordinate;
end

block.OutputPort(1).Data = vehicle_trajectory;