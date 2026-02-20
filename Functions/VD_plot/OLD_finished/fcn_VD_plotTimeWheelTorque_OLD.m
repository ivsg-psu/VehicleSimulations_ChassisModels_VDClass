function fcn_VD_plotTimeWheelTorque(time, wheel_torque, varargin)
%% fcn_VD_plotTimeWheelTorque
% Purpose:
%   To plot the wheel torques against time
%
% Inputs:
%   time: A Nx1 vector of time [sec]
%   wheel_torque: A Nx4 matrix of wheel_torque(s) [Nm]
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_07_03

% Revision history
% 2024_03_21 - S. Brennan
% -- fixed a bug when plotting multiple wheel torques

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

if 3 == nargin
    fig_num = varargin{1};
else
    fig = figure;
    fig_num = fig.Number;
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

max_value = max(wheel_torque, [], 'all');
min_value = min(wheel_torque, [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 25;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotTimeWheelTorque');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf
hold on

plot(time,wheel_torque(:,1),'b','Linewidth',2.4)
plot(time,wheel_torque(:,2),'r','Linewidth',1.2)
plot(time,wheel_torque(:,3),'m--','Linewidth',2.4)
plot(time,wheel_torque(:,4),'g--','Linewidth',1.2)
grid on
legend('FL','FR','RL','RR','Interpreter','latex','Fontsize',13,...
       'NumColumns',2,'Location','best')
set(gca,'Fontsize',13)
ylabel('Wheel Torque $[Nm]$','Interpreter','latex','Fontsize',18)
xlabel('Time $[s]$','Interpreter','latex','Fontsize',18)
ylim([min_value-offset max_value+offset])
xlim([time(1) time(end)])
end