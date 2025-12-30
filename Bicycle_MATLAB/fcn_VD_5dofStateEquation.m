function dydt = fcn_VD_5dofStateEquation(~,y,acceleration, varargin)
%% fcn_VD_5dofStateEquation
%   fcn_VD_5dofStateEquation calculates the derivatives of longitudinal
%   velocity, lateral velocity, yaw rate, and angular velocities of wheels
%   given acceleration and position inputs. The primary purpose of this
%   script is to convert the accelerations (forces) from the Newtonian
%   frame to the body-fixed, non-Newtonian frame, for the body-centered
%   vector components.
%
% FORMAT:
%
%      dydt = fcn_VD_5dofStateEquation(~, y, acceleration, (figNum))
%      dydt ~ [dUdt; dVdt; drdt; domegadt];
%
% INPUTS:
%
%      (unused first input)
%      y: A 5x1 vector of velocities in body coordinates containing
%           [longitudinal velocity; 
%            lateral velocity; 
%            yaw rate; 
%            wheel angular_velocity_Front; 
%            wheel_angular_velocity_Rear]
%
%      acceleration: A 5x1 vector of Newtonian accelerations that contains:
%           [longitudinal acceleration; 
%            lateral acceleration; 
%            yaw acceleration;
%            wheel angular acceleration_Front; 
%            wheel_angular_acceleration_Rear] 
%
%      (OPTIONAL INPUTS)
%
%      figNum: a figure number to plot results.
%
% OUTPUTS:
%
%      dydt: A 5x1 vector of body-fixed accelerations containing:
%           [longitudinal acceleration; 
%            lateral acceleration; 
%            yaw acceleration;
%            wheel angular acceleration_Front; 
%            wheel_angular_acceleration_Rear] 
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_5dofStateEquation
%     for a full test suite.
%
% This function was written on 2021_08_13 by Satya Prasad and is maintained
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% 2021_08_13 by Satya Prasad, szm888@psu.edu
% - Wrote the code originally
% 
% 2025_12_30 by Sean Brennan, sbrennan@psu.edu
% - Fixed header formatting to standard form
% - Fixed input checking to standard form

% TO-DO:
%
% 2025_11_21 by Sean Brennan, sbrennan@psu.edu
% - (fill in items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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
        narginchk(3,MAX_NARGIN);

        % Check the y input to be sure it has 1 column, 5 rows
        fcn_DebugTools_checkInputsToFunctions(y, '1column_of_numbers',[5 5]);

        % Check the acceleration input to be sure it has 1 column, 5 rows
        fcn_DebugTools_checkInputsToFunctions(acceleration, '1column_of_numbers',[5 5]);

    end
end


% 
% % Set the start values
% [flag_start_is_a_point_type, start_zone_definition] = fcn_Laps_checkZoneType(start_zone_definition, 'start_definition', -1);
% 
% 
% % The following area checks for variable argument inputs (varargin)
% 
% % Does the user want to specify the end_definition?
% % Set defaults first:
% end_zone_definition = start_zone_definition; % Default case
% flag_end_is_a_point_type = flag_start_is_a_point_type; % Inheret the start case
% % Check for user input
% if 3 <= nargin
%     temp = varargin{1};
%     if ~isempty(temp)
%         % Set the end values
%         [flag_end_is_a_point_type, end_zone_definition] = fcn_Laps_checkZoneType(temp, 'end_definition', -1);
%     end
% end
% 
% % Does the user want to specify excursion_definition?
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
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; %#ok<NASGU>
        flag_do_plots = 1;
    end
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
dUdt = acceleration(1)+y(3)*y(2); % Vxdot
dVdt = acceleration(2)-y(3)*y(1); % Vydot
drdt = acceleration(3); % yaw acceleration

domegadt = acceleration(4:5); % angular acceleration of wheel

dydt = [dUdt; dVdt; drdt; domegadt];

%% Any debugging?
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
    
    % % plot the final XY result
    % figure(figNum);
    % clf;
    % 
    % % Everything put together
    % subplot(1,2,1);
    % hold on;
    % grid on
    % title('Results of breaking data into laps');
    % 
    % 
    % 
    % % Plot the indices per lap
    % all_ones = ones(length(input_path(:,1)),1);
    % 
    % % fill in data
    % start_of_lap_x = [];
    % start_of_lap_y = [];
    % lap_x = [];
    % lap_y = [];
    % end_of_lap_x = [];
    % end_of_lap_y = [];
    % for ith_lap = 1:Nlaps
    %     start_of_lap_x = [start_of_lap_x; cell_array_of_entry_indices{ith_lap}; NaN]; %#ok<AGROW>
    %     start_of_lap_y = [start_of_lap_y; all_ones(cell_array_of_entry_indices{ith_lap})*ith_lap; NaN]; %#ok<AGROW>;
    %     lap_x = [lap_x; cell_array_of_lap_indices{ith_lap}; NaN]; %#ok<AGROW>
    %     lap_y = [lap_y; all_ones(cell_array_of_lap_indices{ith_lap})*ith_lap; NaN]; %#ok<AGROW>;
    %     end_of_lap_x = [end_of_lap_x; cell_array_of_exit_indices{ith_lap}; NaN]; %#ok<AGROW>
    %     end_of_lap_y = [end_of_lap_y; all_ones(cell_array_of_exit_indices{ith_lap})*ith_lap; NaN]; %#ok<AGROW>;
    % end
    % 
    % % Plot results
    % plot(start_of_lap_x,start_of_lap_y,'g-','Linewidth',3,'DisplayName','Prelap');
    % plot(lap_x,lap_y,'b-','Linewidth',3,'DisplayName','Lap');
    % plot(end_of_lap_x,end_of_lap_y,'r-','Linewidth',3,'DisplayName','Postlap');
    % 
    % h_legend = legend;
    % set(h_legend,'AutoUpdate','off');
    % 
    % xlabel('Indices');
    % ylabel('Lap number');
    % axis([0 length(input_path(:,1)) 0 Nlaps+0.5]);
    % 
    % 
    % subplot(1,2,2);
    % % Plot the XY coordinates of the traversals
    % hold on;
    % grid on
    % title('Results of breaking data into laps');
    % axis equal
    % 
    % cellArrayOfPathsToPlot = cell(Nlaps+1,1);
    % cellArrayOfPathsToPlot{1,1}     = input_path;
    % for ith_lap = 1:Nlaps
    %     temp_indices = cell_array_of_lap_indices{ith_lap};
    %     if length(temp_indices)>1
    %         dummy_path = input_path(temp_indices,:);
    %     else
    %         dummy_path = [];
    %     end
    %     cellArrayOfPathsToPlot{ith_lap+1,1} = dummy_path;
    % end
    % h = fcn_Laps_plotLapsXY(cellArrayOfPathsToPlot,figNum);
    % 
    % % Make input be thin line
    % set(h(1),'Color',[0 0 0],'Marker','none','Linewidth', 0.75);
    % 
    % % Make all the laps have thick lines
    % for ith_plot = 2:(length(h))
    %     set(h(ith_plot),'Marker','none','Linewidth', 5);
    % end
    % 
    % % Add legend
    % legend_text = {};
    % legend_text = [legend_text, 'Input path'];
    % for ith_lap = 1:Nlaps
    %     legend_text = [legend_text, sprintf('Lap %d',ith_lap)]; %#ok<AGROW>
    % end
    % 
    % h_legend = legend(legend_text);
    % set(h_legend,'AutoUpdate','off');
    % 
    % 
    % 
    % %     % Plot the start, excursion, and end conditions
    % %     % Start point in green
    % %     if flag_start_is_a_point_type==1
    % %         Xcenter = start_zone_definition(1,1);
    % %         Ycenter = start_zone_definition(1,2);
    % %         radius  = start_zone_definition(1,3);
    % %         INTERNAL_plot_circle(Xcenter, Ycenter, radius, [0 .7 0], 4);
    % %     end
    % %
    % %     % End point in red
    % %     if flag_end_is_a_point_type==1
    % %         Xcenter = end_definition(1,1);
    % %         Ycenter = end_definition(1,2);
    % %         radius  = end_definition(1,3);
    % %         INTERNAL_plot_circle(Xcenter, Ycenter, radius, [0.7 0 0], 2);
    % %     end
    % %     legend_text = [legend_text, 'Start condition'];
    % %     legend_text = [legend_text, 'End condition'];
    % %     h_legend = legend(legend_text);
    % %     set(h_legend,'AutoUpdate','off');
    % 
    % % Plot start zone
    % h_start_zone = fcn_Laps_plotZoneDefinition(start_zone_definition,'g-',figNum);
    % 
    % % Plot end zone
    % h_end_zone = fcn_Laps_plotZoneDefinition(end_zone_definition,'r-',figNum);

    
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

