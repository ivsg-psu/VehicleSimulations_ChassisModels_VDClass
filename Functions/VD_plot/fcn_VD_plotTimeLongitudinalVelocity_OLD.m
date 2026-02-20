function fcn_VD_plotTimeLongitudinalVelocity(time, longitudinal_velocity, varargin)
%% fcn_VD_plotTimeLongitudinalVelocity
% Purpose:
%   To plot the longitudinal velocity against time
%
% Inputs:
%   time: A Nx1 vector of time [sec]
%   longitudinal_velocity: A Nx1 vector of longitudinal velocity [m/s]
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_07_03

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

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotTimeLongitudinalVelocity');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf
plot(time, longitudinal_velocity,'g','Linewidth',1.2)
grid on
set(gca,'Fontsize',13)
xlabel('Time $[s]$','Interpreter','Latex','Fontsize',18)
ylabel('Longitudinal Velocity $[m/s]$','Interpreter','Latex','Fontsize',18)
end