%%%%%%%%%%%%%%%%%% S-Function sfcn_udpDB2V_receiveWrapper %%%%%%%%%%%%%%%%%
% Purpose:
%   The purpose of this function is to call 'fcn_udpDB2V_receiveFriction' 
%   to execute in the background.
% 
%   Note: Change the sample time depending on the simulation.
% 
% Inputs:
%   Nothing, it receives data over UDP.
% 
% Outputs:
%   Output 1: Time at which query corresponding to the UDP packet is
%   sent.
%   Output 2: Average friction
% 
% Author:  Satya Prasad
% Created: 2022/05/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_udpDB2V_receiveWrapper(block)
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 0;
block.NumOutputPorts = 2;

% Setup functional port to default
block.SetPreCompPortInfoToDefaults;

% Override output port properties
block.OutputPort(1).Dimensions = 1;
block.OutputPort(2).Dimensions = 1;

% Set block sample time to be inherited
block.SampleTimes = [0.01, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output); % Required

function Output(block)
% persistent variables are local to the function in which they are declared,
% yet their values are retained in memory between calls to the function
persistent p F
if isempty(p)
    % get current parallel pool
    % If no parallel pool exists, gcp starts a new parallel pool and 
    % returns a pool object for that
    p = gcp;
    F = parfeval(p,@fcn_udpDB2V_receiveFriction,2);
end

if 0~=size(F.OutputArguments,1)
    % write received data to output variable
    query_time   = F.OutputArguments{1};
    avg_friction = F.OutputArguments{2};
    F = parfeval(p,@fcn_udpDB2V_receiveFriction,2);
else
    query_time   = nan(1);
    avg_friction = nan(1);
end

block.OutputPort(1).Data = query_time;
block.OutputPort(2).Data = avg_friction;