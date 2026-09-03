function fcn_VD_plotTimeFriction(time, friction_coefficient, varargin)
%% fcn_VD_plotTimeFriction
% Purpose:
%   To plot the friction coeeficient against time
%
% Inputs:
%   time: A Nx1 vector of time [sec]
%   friction_coefficient: A Nx2 matrix of friction_coefficient(s) [No units]
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

h_fig = figure(fig_num);
set(h_fig, 'Name', 'fcn_VD_plotTimeFriction');
width = 600; height = 400; right = 100; bottom = 400;
set(gcf, 'position', [right, bottom, width, height])
clf
plot(time,friction_coefficient(:,1),'b.')
hold on
plot(time,friction_coefficient(:,2),'g.')
grid on
legend('Front','Rear','Location','best','Interpreter','Latex','Fontsize',13)
ylabel('Friction Coefficient [No Units]','Interpreter','Latex','Fontsize',13)
xlabel('Time $[s]$','Interpreter','Latex','Fontsize',13)
ylim([-0.1 1.5])
set(gca,'Fontsize',13)
end