%%%%%%%%%%%%%%%%%% S-Function sfcn_udpV2DB_queryFriction %%%%%%%%%%%%%%%%%%
% Purpose:
%   The purpose of this function is to send data from Simulink to 
%   Python/ROS over UDP.
% 
%   Note: Change the sample time depending on the simulation.
% 
% Inputs:
%   Input 1: Flag to decide whether to query or not.
%   Input 2: Simulation time.
%   Input 3: Road ID
%   Input 4: Station at the beginning of 'local'
%   Input 5: Station at the end of 'local'
% 
% Outputs:
%   Nothing. It sends data over UDP.
% 
% Author:  Satya Prasad
% Created: 2022/06/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sfcn_udpV2DB_queryFriction(block)
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%% Function: setup ===================================================
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
% Register number of input and output ports
block.NumInputPorts  = 5;
block.NumOutputPorts = 0;

% Setup functional port properties to dynamically inherited
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions = 1; % Flag to decide whether to query or not
block.InputPort(2).Dimensions = 1; % Time
block.InputPort(3).Dimensions = 1; % Road_id
block.InputPort(4).Dimensions = 1; % Start station
block.InputPort(5).Dimensions = 1; % End station

% Set block sample time to be inherited
block.SampleTimes = [0.01, 0];

% Set the block simStateCompliance to default (i.e., same as a built-in block)
block.SimStateCompliance = 'DefaultSimState';

% Register methods
block.RegBlockMethod('Outputs', @Output); % Required

function Output(block)
query_flag = block.InputPort(1).Data;

if query_flag
    %% Inputs
    query_time    = block.InputPort(2).Data;
    road_id       = block.InputPort(3).Data;
    start_station = block.InputPort(4).Data;
    end_station   = block.InputPort(5).Data;
    
    %% Data to send
    data_to_send = [num2str(query_time) ',' ...
                    num2str(road_id) ',' ...
                    num2str(start_station) ',' ...
                    num2str(end_station)];
    data_encoded = unicode2native(data_to_send, 'UTF-8');
    
    %% Create socket
    remote_host = '10.0.0.160'; % Receiver IP-address
    remote_port = 5002; % Receiver port
    send_socket = udpport(remote_host,remote_port); % create UDP object with 'closed' status
    send_socket.EnablePortSharing = 'on';
    
    %% Open-Send-Close
    fopen(send_socket); % connect/open the socket
    fwrite(send_socket,data_encoded); % send data over UDP
    fclose(send_socket); % disconnect/close the socket
end % NOTE: END IF statement 'query_flag'