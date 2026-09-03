function fcn_VD_plotTimeWheelSpeed(time, wheel_speed, varargin)
%% fcn_VD_plotTimeWheelSpeed
% Purpose:
%   To plot wheel speed(s) against time
%
% Inputs:
%   time: A Nx1 vector of time [sec]
%   wheel_speed: A Nx2 matrix of wheel speed(s) [rad/s]
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_08_13

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
max_value = max(wheel_speed, [], 'all');
min_value = min(wheel_speed, [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 0.5;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotTimeWheelSpeed');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf
plot(time,wheel_speed(:,1),'b','Linewidth',1.2)
hold on
plot(time,wheel_speed(:,2),'g--','Linewidth',1.2)
grid on
legend('Front','Rear','Location','best','Interpreter','Latex','Fontsize',13)
set(gca,'Fontsize',13)
ylabel('Wheel Speed [rad/s]','Interpreter','Latex','Fontsize',18)
xlabel('Time [s]','Interpreter','Latex','Fontsize',18)
ylim([min_value-offset max_value+offset])
end