function fcn_VD_plotTimeSteeringAngle(time, steering_angle, varargin)
%% fcn_VD_plotTimeSteeringAngle
% Purpose:
%   To plot steering angle(s) against time
%
% Inputs:
%   time: A Nx1 vector of time [sec]
%   steering_angle: A Nx4 matrix of steering angle(s) [rad]
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_07_03
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
max_value = max(steering_angle, [], 'all');
min_value = min(steering_angle, [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 0.5;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotTimeSteeringAngle');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf
hold on
plot(time,steering_angle(:,1),'b','Linewidth',2.4)
plot(time,steering_angle(:,2),'r','Linewidth',1.2)
plot(time,steering_angle(:,3),'m--','Linewidth',2.4)
plot(time,steering_angle(:,4),'g--','Linewidth',1.2)
grid on
legend('FL','FR','RL','RR','Interpreter','latex','Fontsize',13,...
       'NumColumns',2,'Location','best')
set(gca,'Fontsize',13)
ylabel('Steering Angle $[rad]$','Interpreter','latex','Fontsize',18)
xlabel('Time $[s]$','Interpreter','latex','Fontsize',18)
ylim([min_value-offset max_value+offset])
xlim([time(1) time(end)])
end