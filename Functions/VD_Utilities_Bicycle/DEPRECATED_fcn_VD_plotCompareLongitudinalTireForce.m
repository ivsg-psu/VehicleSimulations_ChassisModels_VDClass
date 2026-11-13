function fcn_VD_plotCompareLongitudinalTireForce(ref_time,ref_lon_tire_force,...
    ref_sim_name,time,lon_tire_force,sim_name,varargin)
%% fcn_VD_plotCompareLongitudinalTireForce
% Purpose:
%   To plot longitudinal tire force(s) against time
%
% Inputs:
%   ref_time, time: A Nx1 vector of time [sec]
%   ref_lon_tire_force, lon_tire_force: A Nx2 matrix of longitudinal tire force(s) [N]
%   ref_sim_name, sim_name: Simulation name as string
%
% Returned Results:
%   A plot
%
% Author: Satya Prasad
% Created: 2021_08_13
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
max_value = max([ref_lon_tire_force, lon_tire_force], [], 'all');
min_value = min([ref_lon_tire_force, lon_tire_force], [], 'all');
offset    = 0.1*(max_value-min_value);
if 0 == offset
    offset = 25;
end

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotCompareLongitudinalTireForce');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf

subplot(2,1,1)
plot(ref_time, ref_lon_tire_force(:,1), 'b', 'Linewidth', 1.2)
hold on
plot(time, lon_tire_force(:,1), 'r--', 'Linewidth', 1.2)
grid on
legend(['Front ' ref_sim_name], ['Front ' sim_name], 'Location', 'best',...
       'Interpreter','Latex','Fontsize',13)
set(gca,'Fontsize',13)
ylabel('Longitudinal Tire Force $[N]$','Interpreter','Latex','Fontsize',11)
ylim([min_value-offset max_value+offset])

subplot(2,1,2)
plot(ref_time, ref_lon_tire_force(:,2), 'b', 'Linewidth', 1.2)
hold on
plot(time, lon_tire_force(:,2), 'r--', 'Linewidth', 1.2)
grid on
legend(['Rear ' ref_sim_name], ['Rear ' sim_name], 'Location', 'best',...
       'Interpreter','Latex','Fontsize',13)
set(gca,'Fontsize',13)
ylabel('Longitudinal Tire Force $[N]$','Interpreter','Latex','Fontsize',11)
xlabel('Time $[s]$','Interpreter','Latex','Fontsize',13)
ylim([min_value-offset max_value+offset])
end