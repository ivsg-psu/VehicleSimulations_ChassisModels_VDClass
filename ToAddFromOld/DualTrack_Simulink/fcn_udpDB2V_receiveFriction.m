function [query_time,avg_friction] = fcn_udpDB2V_receiveFriction()
%% fcn_udpDB2V_receiveFriction
% Purpose:
%   The purpose of this function it to receive data over UDP.
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

%% Create receive socket
local_host = '10.0.0.160'; % Receiver IP-address
local_port = 5001; % Receiving port
timeout    = 10; % [seconds]
receive_socket = udp('','LocalHost',local_host,'LocalPort',local_port,...
                     'Timeout',timeout); % create UDP object with 'closed' status
receive_socket.EnablePortSharing = 'on';
receive_socket.InputBufferSize   = 1024;

%% Open-Receive-Close
fopen(receive_socket); % open the socket

try
    data_encoded = fread(receive_socket); % read data
    if ~isempty(data_encoded)
        data_decoded  = native2unicode(data_encoded,'UTF-8'); %#ok<N2UNI> % converts bytes to unicode representation
        data_received = str2double(split(data_decoded',','));
        query_time    = data_received(1);
        avg_friction  = data_received(2);
    end % NOTE: END IF statement '~isempty(data_encoded)'
catch ME
    rethrow(ME);
end

fclose(receive_socket); % close the socket
end