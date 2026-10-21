function [stateTrajectory, t, steeringUsed] = fcn_VD_2DOFBicycleModelRK4(initialStates, deltaT, timeInterval, inputsVsTime, vehicleParameters, varargin)

%% fcn_VD_2DOFBicycleModelRK4
%   Simulates the bicycle kinematic model using Runga Kutta 4th-order
%
% FORMAT:
%
%      [stateTrajectory, t, steeringUsed] =
%      fcn_VD_2DOFBicycleModelRK4(initialStates, deltaT,
%      timeInterval, inputsVsTime, vehicleParameters, (figNum))
%
% INPUTS:
%
%      initialStates: A 1x3 vector of inital global pose in form of
%         [X Y Phi], which stand for:
% 
%         X: Global X position in meters
%
%         Y: Global Y position in meters
%
%         phi: Global yaw angle in radians, measured positive from X axis
%         to Y axis
%
%      deltaT: a 1x1 positive number denoting the time step to use, in
%      seconds
%
%      timeInterval: a 1x2 vector denoting [startTime endTime] in seconds
%
%      inputsVsTime: a Mx2 vector denoting 
%      [steeringTime steeringValues] 
%      in units of [sec rad] respectively. This is interpolated using
%      linear interpolation at the sampling times. For times outside the
%      given interval, zero values are used.
%
%      vehicleParameters: a structure with the following subfields
%          % Example below is approximately a Ford Taurus
%            vehicleParameters(1).U          = 20;     % x-axis velocity, in m/s
%            vehicleParameters(1).m          = 1031.9; % mass in kg
%            vehicleParameters(1).Iz         = 1850;   % z-axis rotational inertia in kg-m^2
%            vehicleParameters(1).a          = 0.9271; % distance from front axle to CG, in meters
%            vehicleParameters(1).b          = 1.5621; % distance from rear axle to CG, in meters
%            vehicleParameters(1).Caf        = 77500;  % Cornering stiffness of front axle, N/rad;
%            vehicleParameters(1).Car	     = 116250; % Cornering stiffness of rear axle, N/rad;
%            
%      (OPTIONAL INPUTS)
%
%      figNum: a FID number to print results. If set to -1, skips any
%      input checking or debugging, no prints will be generated, and sets
%      up code to maximize speed.
%
% OUTPUTS:
%
%      stateTrajectory: An Nx3 vector of the state trajectory, with the
%      columns as [X Y Phi] in units of [m],[m],[rad]
%
%      t: An Nx1 vector of the simulation times, in seconds
%
%      steeringUsed: An Nx1 vector of the steering values used in the sim,
%      in units of [rad]
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%      fcn_VD_derivatives2DOFBicycleModel
%
% EXAMPLES:
%
%     See the script: script_test_fcn_VD_2DOFBicycleModelRK4
%     for a full test suite.
%
% This function was written on 2026_10_20 
% by Sean Brennan. Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
% 
% As: fcn_VD_2DOFBicycleModelRK4
%
% 2026_10_20 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_2DOFBicycleModelRK4
%   % * Used fcn_VD_kinematicBicycleModelRK4 as starter
%   % * Copied material over from fcn_VD_bicycle2dofModel

% TO-DO:
% - 2026_10_20 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 6; % The largest Number of argument inputs to the function
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
        narginchk(MAX_NARGIN-1,MAX_NARGIN);

        % Validate that the initialStates input has 3 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(initialStates, '3column_of_numbers',[1 1]);

        % Validate that the deltaT input has 1 column, 1 row
        fcn_DebugTools_checkInputsToFunctions(deltaT, '1column_of_numbers',[1 1]);

        % Validate that the timeInterval input has 2 columns, 1 row
        fcn_DebugTools_checkInputsToFunctions(timeInterval, '2column_of_numbers',[1 1]);

        % Validate that the inputsVsTime input has 2 columns, 2+
		% rows
        fcn_DebugTools_checkInputsToFunctions(inputsVsTime, '2column_of_numbers',[2 3]);

		% Check the vehicleParameters input is a structure
        assert(isstruct(vehicleParameters));

    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp; 
        flag_do_plots = 1;
    end
end


%% Implements Bicycle Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RK4 in MATLAB Script
startTime = timeInterval(1);
endTime = timeInterval(2);
simulationTimes = (startTime:deltaT:endTime)';
N_timeSteps = length(simulationTimes); % This is the number of time steps we should have

% The number of vehicles is the length of the structure array
N_vehicles = length(vehicleParameters);


%% Do simulations
% This section runs the codes to predict core behavior within SIMULINK and MATLAB,
% looping through each vehicle.
% For this simulation set, the core behavior is
% V: lateral speed
% r: rotational rate of vehicle around z-axis

if N_vehicles>1
	stateTrajectory = cell(N_vehicles,1);
else
	stateTrajectory = nan(N_timeSteps,3);
end

% Loop through all the vehicles, simulating the trajectory of each within
% the for loop
for ith_vehicle = 1:N_vehicles
    % Print to the console which vehicle we are working on
    fprintf(1,'Working on vehicle: %d\n', ith_vehicle);

	% Initialize variables
	thisVehicleParameters = parameters(ith_vehicle);
	thisStateTrajectory   = nan(N_timeSteps,3);
	thisT = nan(N_timeSteps,1);
	thisSteeringUsed = interp1(inputsVsTime(:,1), inputsVsTime(:,2), simulationTimes,'linear',0);

    %% RK4 in MATLAB Script
    % Fill in the parameters needed by the simulation for it to run

    % Set initial conditions
    currentStates = initialStates;

    for ith_time = 1:N_timeSteps
		thisTime       = simulationTimes(ith_time);

		% Fill in the results to save
		t(ith_time,1)  = thisTime; % Update time
		thisStateTrajectory(ith_time,:)   = currentStates;

		% Use Runga-Kutta to predict next position
		y = currentStates';
		steering_amplitude_degrees = thisSteeringUsed(ith_time,1);
        [~, y] = fcn_VD_RungeKutta(...
                @(t,y) fcn_VD_derivatives2DOFBicycleModel(y, steering_amplitude_degrees, thisVehicleParameters), ...
                currentStates', thisTime, deltaT, -1);
       
		currentStates = y';
    end
end

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
    
    % plot the outputs
	h_plot = fcn_VD_plotTrajectory(stateTrajectory(:,1:2),(figNum));
	set(h_plot,'DisplayName','XY Trajectory (MATLAB RK4)')

	% % Start the plotting. Both the plots for each situation should agree with
	% % each other (otherwise the equations are not consistent)
	%
	% index_to_label = (2==t);
	% % Loop through each of the vehicles
	% for ith_vehicle=1:N_vehicles
	%     % Plot the yaw rate
	%     h1 = figure(99);
	%     hold on;
	%     set(h1,'Name','Yawrate')
	%     plot(t(:,ith_vehicle),all_r(:,ith_vehicle),'b');
	%     grid on;
	%     text(2, all_r(index_to_label,ith_vehicle), ['vehicle' num2str(ith_vehicle)]);
	%     xlabel('Time (sec)');
	%     ylabel('Yawrate (rad/sec)');
	%     title('Yawrate');
	%
	%     % Plot the lateral velocity
	%     h2 = figure(88);
	%     hold on;
	%     set(h2,'Name','LatVel')
	%     plot(t(:,ith_vehicle),all_V(:,ith_vehicle),'b');
	%     grid on;
	%     text(2, all_V(index_to_label,ith_vehicle), ['vehicle' num2str(ith_vehicle)]);
	%     xlabel('Time (sec)');
	%     ylabel('Lateral Velocity (m/sec)');
	%     title('Lateral Velocity');
	%
	%     % The XY Plots
	%     h3 = figure(77);
	%     hold on;
	%     set(h2,'Name','XYposition')
	%     plot(all_X(:,ith_vehicle),all_Y(:,ith_vehicle),'b');
	%     grid on;
	%     text(all_X(index_to_label,ith_vehicle), all_Y(index_to_label,ith_vehicle), ...
	%         ['vehicle' num2str(ith_vehicle)]);
	%     xlabel('X position [m]');
	%     ylabel('Y position [m]');
	%     title('Position');
	% end
	%
	% % The Steering Input
	% h4 = figure(66);
	% hold on;
	% set(h4,'Name','Steering Input')
	% plot(t(:,ith_vehicle),df,'b');
	% grid on;
	% xlabel('Time (sec)');
	% ylabel('Steering input (rad)');
	% title('Steering Input');
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


