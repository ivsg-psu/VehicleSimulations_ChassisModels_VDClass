function fcn_VD_plotStationFrictionUtilization(wheel_station,friction_utilization,varargin)
%% fcn_VD_plotStationFrictionUtilization
% Purpose:
%   To plot friction utilization against wheel station
%
% Inputs:
%   wheel_station: A Nx4 matrix of wheel station [m]
%   steering_angle: A Nx4 matrix of friction utilization [No Units]
% 
% Returned Results:
%   A plot
% 
% Author: Satya Prasad
% Created: 2021_011_10
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
if 2>nargin || 3<nargin
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
if 3 == nargin
    fig_num = varargin{1};
else
    fig = figure;
    fig_num = fig.Number;
end

max_value = max(friction_utilization, [], 'all');
min_value = min(friction_utilization, [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 0.5;
end

max_value_x = max(wheel_station, [], 'all');
min_value_x = min(wheel_station, [], 'all');
offset_x    = 0.1*(max_value_x-min_value_x);
if 0 == offset_x
    offset_x = 0.5;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotStationFrictionUtilization');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

subplot(2,2,1)
plot(wheel_station(:,1), friction_utilization(:,1), 'b', 'Linewidth', 1)
grid on
legend('front left wheel', 'Location', 'best')
ylabel('Friction Utilization')
ylim([min_value-offset max_value+offset])
xlim([min_value_x-offset_x max_value_x+offset_x])

subplot(2,2,2)
plot(wheel_station(:,2), friction_utilization(:,2), 'b', 'Linewidth', 1)
grid on
legend('front right wheel', 'Location', 'best')
ylim([min_value-offset max_value+offset])
xlim([min_value_x-offset_x max_value_x+offset_x])

subplot(2,2,3)
plot(wheel_station(:,3), friction_utilization(:,3), 'b', 'Linewidth', 1)
grid on
legend('rear left wheel', 'Location', 'best')
ylabel('Friction Utilization')
xlabel('Station [m]')
ylim([min_value-offset max_value+offset])
xlim([min_value_x-offset_x max_value_x+offset_x])

subplot(2,2,4)
plot(wheel_station(:,3), friction_utilization(:,4), 'b', 'Linewidth', 1)
legend('rear right wheel', 'Location', 'best')
grid on
xlabel('Station [m]')
ylim([min_value-offset max_value+offset])
xlim([min_value_x-offset_x max_value_x+offset_x])

sgtitle('Friction Utilization')
end