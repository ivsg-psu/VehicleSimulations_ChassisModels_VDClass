function fcn_VD_plotCompareTrajectory(ref_trajectory,ref_sim_name,...
    trajectory,sim_name,varargin)
%% fcn_VD_plotCompareTrajectory
% Purpose:
%   To plot the trajectory
%
% Inputs:
%   ref_trajectory, trajectory: A Nx2 vector of trajectory
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
if 4>nargin || 5<nargin
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
if 5 == nargin
    fig_num = varargin{1};
else
    fig = figure;
    fig_num = fig.Number;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotCompareTrajectory');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

plot(ref_trajectory(:,1), ref_trajectory(:,2), 'b', 'Linewidth', 1.2)
hold on
plot(trajectory(:,1), trajectory(:,2), 'r--', 'Linewidth', 1.2)
legend(ref_sim_name, sim_name, 'Location', 'best',...
       'Interpreter','Latex','Fontsize',13)
axis equal
grid on
set(gca,'Fontsize',13)
xlabel('East $[m]$','Interpreter','Latex','Fontsize',18)
ylabel('North $[m]$','Interpreter','Latex','Fontsize',18)
end