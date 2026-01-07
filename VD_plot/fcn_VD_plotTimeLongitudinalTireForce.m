function fcn_VD_plotTimeLongitudinalTireForce(time, longTireForces, varargin)
%% fcn_VD_plotTimeLongitudinalTireForce
%   fcn_VD_plotTimeLongitudinalTireForce plots the longitudinal tire
%   force(s) versus time where the columns are long forces of each tire.
%
% FORMAT:
%
%       fcn_VD_plotTimeLongitudinalTireForce(time, longTireForces, (figNum))
%
% INPUTS:
%
%     time: an Nx1 vector representing time [sec]
%
%     longTireForces: an NxM vector representing the longitudinal tire
%     force(s) of M tires [unitless], in the following order: FL, FR, RL,
%     RR if 4 columns, or F and R if 2 columns
%
%     (OPTIONAL INPUTS)
%
%     figNum: a figure number to plot results.
%
% OUTPUTS:
%
%     (none)
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_plotTimeLongitudinalTireForce
%     for a full test suite.
%
% This function was written on 2021_08_13 by Satya Prasad and is maintained
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% 2021_07_03 by Satya Prasad, szm888@psu.edu
% - Wrote the code originally
% 
% 2026_01_07 by Sean Brennan, sbrennan@psu.edu
% - Fixed header formatting to standard form
% - Renamed "slip_+angle" variable to be longTireForces

% TO-DO:
%
% 2026_01_07 by Sean Brennan, sbrennan@psu.edu
% - Need to fill in the tire ordering

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
flag_max_speed = 0; % The default. This runs code with all error checking
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_VD_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_VD_FLAG_CHECK_INPUTS");
    MATLABFLAG_VD_FLAG_DO_DEBUG = getenv("MATLABFLAG_VD_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VD_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VD_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VD_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VD_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug % If debugging is on, print on entry/exit to the function
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_figNum = 999978; %#ok<NASGU>
else
    debug_figNum = []; %#ok<NASGU>
end

%% check input arguments?
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
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,MAX_NARGIN);

        % Check the time input to be sure it has 1 column and 1 row, positive
        fcn_DebugTools_checkInputsToFunctions(time, 'positive_1column_of_numbers',[1 2]);
        Ndata = size(time,1);

        % Check the longTireForces input to be sure it has 1 column and 1 row, positive
        fcn_DebugTools_checkInputsToFunctions(longTireForces, '1orMorecolumn_of_numbers',[Ndata Ndata]);

    end
end


% 
%   Set the start values
% [flag_start_is_a_point_type, start_zone_definition] = fcn_Laps_checkZoneType(start_zone_definition, 'start_definition', -1);
% 
% 
%   The following area checks for variable argument inputs (varargin)
% 
%   Does the user want to specify the end_definition?
%   Set defaults first:
% end_zone_definition = start_zone_definition; % Default case
% flag_end_is_a_point_type = flag_start_is_a_point_type; % Inheret the start case
%   Check for user input
% if 3 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         % Set the end values
%         [flag_end_is_a_point_type, end_zone_definition] = fcn_Laps_checkZoneType(temp, 'end_definition', -1);
%     end
% end
% 
%   Does the user want to specify excursion_definition?
% flag_use_excursion_definition = 0; % Default case
% flag_excursion_is_a_point_type = 1; % Default case
% if 4 <= nargin
%     temp = varargin{2};
%     if ~isempty(temp)
%         % Set the excursion values
%         [flag_excursion_is_a_point_type, excursion_definition] = fcn_Laps_checkZoneType(temp, 'excursion_definition',-1);
%         flag_use_excursion_definition = 1;
%     end
% end



% Does user want to show the plots?
flag_do_plots = 1; % Default is to NOT show plots
figNum = [];
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; 
        flag_do_plots = 1;
    end
end

if isempty(figNum)
    fig = figure;
    figNum = fig.Number;
end

%% Differential equation for 5-DOF model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntires = size(longTireForces,2);

%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots

    max_value = max(longTireForces, [], 'all');
    min_value = min(longTireForces, [], 'all');
    offset    = 0.1*(max_value-min_value);
    if 0 == offset
        offset = 0.5;
    end

    h_fig = figure(figNum);
    set(h_fig, 'Name', 'fcn_VD_plotTimeLongitudinalTireForce');
    % width = 600; height = 400; right = 100; bottom = 400;
    % set(gcf, 'position', [right, bottom, width, height])
    clf
    hold on

    if 4==Ntires
        colorStrings = {'b','r','m--','g--'};
        lineWidths = [2.4; 1.2; 2.4; 1.2];
        labelStrings = {'FL','FR','RL','RR'};
    elseif 2== Ntires
        colorStrings = {'b','r'};
        lineWidths = [2.4; 2.4];
        labelStrings = {'F','R'};
    else
        error('Plotting not set up yet to plot %.0f tires. Throwing an error.',Ntires);
    end

    for ith_tire = 1:Ntires
        plot(time,longTireForces(:,ith_tire),colorStrings{ith_tire},'Linewidth',lineWidths(ith_tire), 'DisplayName', labelStrings{ith_tire});
    end

    grid on
    legend('Interpreter','latex','Fontsize',13,...
        'NumColumns',1,'Location','best')

    set(gca,'Fontsize',13)
    ylabel('Long. Tire Force $[N]$','Interpreter','latex','Fontsize',18)
    xlabel('Time $[s]$','Interpreter','latex','Fontsize',18)
    ylim([min_value-offset max_value+offset])
    xlim([time(1) time(end)])


    
end

if flag_do_debug
    fprintf(1, 'ENDING function: %s, in file: %s\n\n', st(1).name, st(1).file);
end

end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§



