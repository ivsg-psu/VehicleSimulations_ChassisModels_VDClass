function fcn_VD_plotCompareSteeringAngle(ref_time,ref_steering_angle,ref_sim_name,...
    time,steering_angle,sim_name,varargin)
%% fcn_VD_plotCompareSteeringAngle
% Purpose:
%   To plot steering angle(s) against time
%
% Inputs:
%   ref_time, time: A Nx1 vector of time [sec]
%   ref_steering_angle, steering_angle: A Nx4 matrix of steering angle(s) [rad]
%   ref_sim_name, sim_name: Simulation name as string
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_07_05
% 

%% Check input arguments
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
% Are there the right number of inputs?
if 6>nargin || 7<nargin
    error('Incorrect number of input arguments')
end

%% Plots the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 7==nargin
    fig_num = varargin{1};
else
    fig = figure;
    fig_num = fig.Number;
end
max_value = max([ref_steering_angle, steering_angle], [], 'all');
min_value = min([ref_steering_angle, steering_angle], [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 0.5;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotCompareSteeringAngle');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

subplot(2,2,1)
plot(ref_time, ref_steering_angle(:,1), 'b', 'Linewidth', 1)
hold on
plot(time, steering_angle(:,1), 'r', 'Linewidth', 1)
grid on
legend(['FL ' ref_sim_name], ['FL' sim_name], 'Location', 'best')
ylabel('Steering Angle [rad]')
ylim([min_value-offset max_value+offset])

subplot(2,2,2)
plot(ref_time, ref_steering_angle(:,2), 'b', 'Linewidth', 1)
hold on
plot(time, steering_angle(:,2), 'r', 'Linewidth', 1)
grid on
legend(['FR ' ref_sim_name], ['FR ' sim_name], 'Location', 'best')
ylim([min_value-offset max_value+offset])

subplot(2,2,3)
plot(ref_time, ref_steering_angle(:,3), 'b', 'Linewidth', 1)
hold on
plot(time, steering_angle(:,3), 'r', 'Linewidth', 1)
grid on
legend(['RL ' ref_sim_name], ['RL ' sim_name], 'Location', 'best')
ylabel('Steering Angle [rad]')
xlabel('Time [s]')
ylim([min_value-offset max_value+offset])

subplot(2,2,4)
plot(ref_time, ref_steering_angle(:,4), 'b', 'Linewidth', 1)
hold on
plot(time, steering_angle(:,4), 'r', 'Linewidth', 1)
grid on
legend(['RR ' ref_sim_name], ['RR ' sim_name], 'Location', 'best')
xlabel('Time [s]')
ylim([min_value-offset max_value+offset])

sgtitle('Steering Angle')
end