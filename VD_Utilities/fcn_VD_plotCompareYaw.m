function fcn_VD_plotCompareYaw(ref_time,ref_yaw,ref_sim_name,...
    time,yaw,sim_name,varargin)
%% fcn_VD_plotCompareYaw
% Purpose:
%   To plot the yaw against time
%
% Inputs:
%   ref_time, time: A Nx1 vector of time [sec]
%   ref_yaw, yaw: A Nx1 vector of yaw [rad]
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
set(h_fig, 'Name', 'fcn_VD_plotCompareYaw');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

plot(ref_time, ref_yaw, 'b', 'Linewidth', 1.2)
hold on
plot(time, yaw, 'r--', 'Linewidth', 1.2)
legend(ref_sim_name, sim_name, 'Location', 'best',...
       'Interpreter','Latex','Fontsize',13)
grid on
set(gca,'Fontsize',13)
xlabel('Time $[s]$','Interpreter','Latex','Fontsize',18)
ylabel('Yaw $[rad]$','Interpreter','Latex','Fontsize',18)
end