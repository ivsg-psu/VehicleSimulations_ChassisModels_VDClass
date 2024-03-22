function fcn_VD_plotCompareLongitudinalVelocity(ref_time,ref_longitudinal_velocity,...
    ref_sim_name,time,longitudinal_velocity,sim_name,varargin)
%% fcn_VD_plotCompareLongitudinalVelocity
% Purpose:
%   To plot the longitudinal velocity against time
%
% Inputs:
%   ref_time, time: A Nx1 vector of time [sec]
%   ref_longitudinal_velocity, longitudinal_velocity: A Nx1 vector of 
%   longitudinal velocity [m/s]
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
if 7 == nargin
    fig_num = varargin{1};
else
    fig = figure;
    fig_num = fig.Number;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotCompareLongitudinalVelocity');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

plot(ref_time, ref_longitudinal_velocity, 'b', 'Linewidth', 1.2)
hold on
plot(time, longitudinal_velocity, 'r--', 'Linewidth', 1.2)
legend(ref_sim_name, sim_name, 'Location', 'best',...
       'Interpreter','Latex','Fontsize',13)
grid on
set(gca,'Fontsize',13)
xlabel('Time $[s]$','Interpreter','Latex','Fontsize',18)
ylabel('Longitudinal Velocity $[m/s]$','Interpreter','Latex','Fontsize',18)
end