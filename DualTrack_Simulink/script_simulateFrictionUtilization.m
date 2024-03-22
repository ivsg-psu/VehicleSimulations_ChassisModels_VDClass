%%%%%%%%%%% Script script_run_mdl_VD_dt7dofModelCosimulationST.m %%%%%%%%%%
% Author:  Satya, Liming, Dr. Beal, Dr.Brennan
% Created: 2021/10/23

% To do list:
% 6) Concurrently, and building off of the work done by PI Beal [26], 
%   this effort will extend the characterization from direct friction 
%   sensing from state-of-the-art wheel force transducers to instead use 
%   indirect friction estimates obtained from load cells in the steering 
%   tie rods and/or torque measurements from the steering servomotors.

% Processing Flow:
%   Step 1: Prepare the workspace
%   Step 2: Add path to dependencies and define inputs
%   Step 3: Query for ENU reference
%           refLLA = fcn_ENUrefPointQuery();
%   Step 4: Get valid Unix Time and bounds on System Entrance Time
%           unixTimeAndAimsunTimeRange = fcn_findValidTimeRangeUnixAndAimsun();
%   Step 5: Get the valid combination of Vehilce ID and Global Time
%           [vehiclePopulation_ID_GlobalTime,SectionId_VehID_GlobalTime] = ...
%                   fcn_findValidVehicleIdandGlobalTime();
%   Step 6: Get centerline data (section/lane)
%   Step 7: Load a vehicle's trajectory
%           queryVehicleTrajectory = fcn_queryVehicleTrajectory();
%   Step 8: Add offset and swerving to the queryVehicleTrajectory
%           swervedVehicleTrajectory = fcn_swervedVehicleTrajectory();
%   Step 9: Run the simulation
%   Step 10: Plot the results
%   Step 11: Insert data into database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 1: Prepare the workspace
clear all   %#ok<CLALL> % Clears workspace to remove any old variables
close all   % Close all open figures
clc         % Clear console space

%% Step 2: Add path to dependencies and define inputs
addpath('C:\Users\prasa\Documents\GitHub\IVSG\Vehicle Simulations\VehicleSimulations_ProjectExamples_NSFForgetfulDatabases\Utilities'); % all the functions and wrapper class
addpath('C:\Users\prasa\Documents\GitHub\IVSG\Vehicle Simulations\VehicleSimulations_ProjectExamples_NSFForgetfulDatabases\DataFiles'); % all the .mat data files

dir.datafiles = ['.' filesep 'DataFiles' filesep]; % all the .mat data files

load('sc2pm'); load('pm2sc');
friction_ref_table = sc_to_pm_table;
% 'names' is a structure containing name(s) of database, tables/relations
% -------------------------------------------------------------------------
% name of the road (Look for this in the table 'road')
names.road = 'I99_StateCollege_To_PortMatilda';
% names.road = 'Onekm_curved_highway';
simulink_trip_id = 1;

% name of the table containing trip information
names.tableTrips = 'trips';
names.trip = 'aimsun simluation 2021-02-24'; % Look for this in the table 'trips'
% names.trip = 'aimsun simluation 2020-05-13';
% -------------------------------------------------------------------------
% name of the database
names.database         = 'nsf_roadtraffic_friction_v2';
% name of the table that contains the relation between road name and road id
names.tableRoad        = 'road';
% name of the table that contains the relation between road id and
% section/segment ids (Aimsun section is equivalent to Database segment)
names.tableRoadSegment = 'road_road_segment';
% name of the table containing processed traffic data traversing between 
% State College and Port Matilda
names.tableTraffic     = 'road_traffic_procsessed_sctopm';
% names.tableTraffic     = 'road_traffic_procsessed';

spheroid = referenceEllipsoid('wgs84'); % referenceEllipsoid

% flag triggers
% set to 'true' to query new data from DB. Otherwise Onekm_curved_highway default file will be loaded.
flag.dbQuery  = true;
% flag to print to the console
flag.verbose  = true;
% set to 'true' to print trajectory information to command window
flag.doDebug  = false;
% flag for plot
flag.plot     = false;

% swerve parameters for adding two sine waves
flag.swerve        = false;  % 'true' to add distortion and 'false' to not
flag.lateralOffset = false;  % 'true' to add random offset of the trajectory and 'false' to not
% amplitude(s) of sine distortion in [meters]
% perpendicular to center of the aimsun vehicle trajectory
swerve.amplitude  = [0.1, 0.4]; %[meters] [0.3,0.1]
% time period(s) of sine distortion in [seconds]
swerve.timePeriod = [37, 4]; %[seconds][60,5]
% maximum fixed lateral offset to a vehicle trajectory
swerve.MaxLateralOffset = 0.1; % [meters]

% field names of vehicle trajectory before being used in the simulation
fieldsTrajectory.vehId      = 1;
fieldsTrajectory.vehType    = 2;
fieldsTrajectory.vehLength  = 3;
fieldsTrajectory.vehWidth   = 4;
fieldsTrajectory.sectiondId = 5;
fieldsTrajectory.junctionId = 6;
fieldsTrajectory.laneNum    = 7;
fieldsTrajectory.direction  = 8;
fieldsTrajectory.lat        = 9;
fieldsTrajectory.lon        = 10;
fieldsTrajectory.alt        = 11;
fieldsTrajectory.east       = 12;
fieldsTrajectory.north      = 13;
fieldsTrajectory.up         = 14;
fieldsTrajectory.yaw        = 15;
fieldsTrajectory.speed      = 16;
fieldsTrajectory.friction_station = 17;
fieldsTrajectory.station    = 18;
fieldsTrajectory.aimsunTime = 19;
fieldsTrajectory.globalTime = 20;

% Define vehicle properties
% Define a MATLAB structure that specifies the physical values for a vehicle.
% For convenience, we ask that you call this stucture 'vehicle'.
vehicle.m   = 1600; % mass (kg)
vehicle.Izz = 2500; % mass moment of inertia (kg m^2)
vehicle.Iw  = 1.2; % mass moment of inertia of a wheel (kg m^2)
vehicle.Re  = 0.32; % effective radius of a wheel (m)
vehicle.h_cg = 0.42; % height of the cg (m)
vehicle.Ca  = [95000; 95000; 110000; 110000]; % wheel cornering stiffnesses
vehicle.Cx  = [65000; 65000; 65000; 65000]; % longitudinal stiffnesses

vehicle.contact_patch_length = 0.15; % [meters]
vehicle.friction_ratio = 1; % [No units]

vdParam.longitudinalTransfer = 1;
vdParam.lateralTransfer = 1;

controller.look_ahead_distance = 20; % look-ahead distance [meters]
controller.steering_Pgain = 0.1; % P gain for steering control
controller.velocity_Pgain = 200; % P gain for steering control

mdlParam.flagFriction = false; % true->Add noise; false->No noise

%% Friction parameters
mdlParam.stationCoeff = 0.0;
mdlParam.timeCoeff = 0.0;
mdlParam.frictionIntercept = 0.9;

%% Get the current UTC time and GPS time
gps_utc_time = 18; % [seconds]
dt = datetime('now','TimeZone','America/New_York');
matsim_unix_time = posixtime(dt); % [seconds]
matsim_gps_time  = matsim_unix_time+gps_utc_time; % [seconds]

%% Step 3: Query for ENU reference
% reference LLA for ENU to LLA transformation and viceversa
if flag.dbQuery
    % query reference point from dababase
    refLLA = fcn_ENUrefPointQuery(names.road);
else
    % NEMO Point, Centre of pacific ocean, is used as reference for virtual
    % simulations
    refLLA = [-48.876667, -123.393333, 0, 0];
end % NOTE: Ends flag.dbQuery (Step 3)
mdlParam.refLLA  = refLLA(1:3);
enu_reference_id = refLLA(4);

%% Step 4: Get valid Unix Time and bounds on System Entrance Time
% Unix Time is the number of seconds that have elapsed since 1970-01-01 00:00:00.
% It is computed from datetime using the function 'posixtime'. It is in [seconds]

% System Entrance Time is the time at which a vehicle enters the network.
% It is not wall time but the time measured by AIMSUN or time elapsed in
% the simulation. It is in [seconds]

% Global Time is Unix Time plus System Entrance Time

% unixTimeAndAimsunTimeRange: It's a Nx5 table. First and second columns
% contian trip id and trip description. Third column contains Unix Time.
% Fourth and fifth columns contains lower and upper bounds on System
% Entrance Time respectively. Output contains information of all trips if
% number of inputs is one. Output contains information of all trips on
% 'names.road' if number of inputs is two.
if flag.dbQuery
    unixTimeAndAimsunTimeRange = ...
        fcn_findValidTimeRangeUnixAndAimsun(names,names.road);
    
    indexTripOfInterest = 1;
    % UNIX Time at which simulation is started, [seconds]
    time.Unix = unixTimeAndAimsunTimeRange(indexTripOfInterest,'UnixTime');
    % lower bound on the System Entrance Time, [seconds]
    time.AimsunLB = unixTimeAndAimsunTimeRange(indexTripOfInterest,'aimsunLB');
    % upper bound on the System Entrance Time, [seconds]
    time.AimsunUB = unixTimeAndAimsunTimeRange(indexTripOfInterest,'aimsunUB');
else
    indexTripOfInterest = 0;
    time.Unix = 0;
    time.AimsunLB = 0;
    time.AimsunUB = 0; 
end % NOTE: Ends flag.dbQuery (Step 4)

%% Step 5: Get all the valid Vehicle ID, Global Time combinations lying
% within specified bounds of global time defined by Unix Time and System
% Entrance Time
% vehiclePopulation_ID_GlobalTime: Nx2 vector, consisting of unique
% (vehicle_id, global_time) pairs. N is number of unique trajectories.

% SectionId_VehID_GlobalTime: Mx3 vector, consisting of unique
% (section_id, vehicle_id, global_time) pairs. M is number of unique
% encounters of all trajectories with sections in a road.

% listOfSections: Lx1 vector,consisting of unique section_ids
if flag.dbQuery
    [~,SectionId_VehID_GlobalTime,listOfSections] = ...
        fcn_findValidVehicleIdandGlobalTime(time.Unix,time.AimsunLB,...
        time.AimsunUB,names,names.road);
else
    % load the default trajectory
    SectionId_VehID_GlobalTime = [201 593 1589418271];
    listOfSections = 201;
end % NOTE: Ends flag.dbQuery (Step 5)

%% Go over all the trajectories in each section
% NOTE: There are only two sections in I99 from SC to PM
for indexSectionOfInterest = 1:length(listOfSections)
    if 1 == mod(indexSectionOfInterest,2)
        mdlParam.sideOfRoad = 1;
    else
        mdlParam.sideOfRoad = -1;
    end
    %% Step 6: Get centerline data
    % find road section from SectionId_VehID_GlobalTime given a vehicle_id
%     if flag.dbQuery
%         % query lane center data
%         [~,lanesCenter_table] = ...
%             fcn_laneCenterQuerybySection(listOfSections(indexSectionOfInterest));
%         % query road segment reference data
%         [~,sectionRef_table] = ...
%             fcn_sectionRefQuerybySection(listOfSections(indexSectionOfInterest));
%     else
%         % load lane center data
%         load('lanesCenter_table_section201.mat','lanesCenter_table');
%         % load road segment reference data
%         load('sectionRef_table_section201.mat','sectionRef_table');
%     end % NOTE: Ends flag.dbQuery (Step 6)
%     
%     % prepare data for ST Road Geometry DB Query block
%     road_geometry_radius = 1./sectionRef_table.curvature;
%     road_geometry = [sectionRef_table.east, sectionRef_table.north, ...
%         sectionRef_table.up, sectionRef_table.grade, sectionRef_table.bank, ...
%         road_geometry_radius, sectionRef_table.yaw];
    if 1 == indexSectionOfInterest
        sectionRef_table = sc_to_pm_table;
    else
        sectionRef_table = pm_to_sc_table;
    end
    road_properties.grade = 0; road_properties.bank_angle = 0; % road properties
    
    % find vehiclePopulation_ID_GlobalTime by section
    vehiclePopulation_ID_GlobalTime = ...
        SectionId_VehID_GlobalTime((SectionId_VehID_GlobalTime(:,1)==...
        listOfSections(indexSectionOfInterest)),2:3);
    refGlobalTime = min(SectionId_VehID_GlobalTime(:,3));
    for indexTrajectoryOfInterest = 1:size(vehiclePopulation_ID_GlobalTime,1)
        %% Step 7: Load a vehicle's trajectory
        % queryVehicleTrajectory: Nx20 matrix, containing all the attributes
        % defined by trajectory_attributes(within the function) sorted in
        % the order of aimsun_time
        mdlParam.vehGlobalTime = ...
            vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,2)-...
            refGlobalTime;
        if flag.dbQuery
            % query vehicle trajectory from the database
            queryVehicleTrajectory = fcn_queryVehicleTrajectory(...
                vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,1), ...
                vehiclePopulation_ID_GlobalTime(indexTrajectoryOfInterest,2), ...
                names, listOfSections(indexSectionOfInterest));
        else
            % load queryVehicleTrajectory data, Just in case to dubug it offline
            load('vehicle_trajectory_593_1589418271.mat', 'queryVehicleTrajectory');
        end % NOTE: Ends flag.dbQuery (Step 7)
        
        if flag.verbose
            % print the duration of vehicle trajectory
            fprintf(1,'\nDuration of vehicle - %d trajectory is %.2f seconds \n', ...
                queryVehicleTrajectory(1,1), ...
                max(queryVehicleTrajectory(:,fieldsTrajectory.aimsunTime))-...
                min(queryVehicleTrajectory(:,fieldsTrajectory.aimsunTime)));
            % print the start and end point of vehicle trajectory
            fprintf(1,'Trajectory of vehicle - %d begins at (%.2f, %.2f, %.2f) and ends at (%.2f, %.2f, %.2f) \n', ...
                queryVehicleTrajectory(1,1), ...
                queryVehicleTrajectory(1,fieldsTrajectory.east), ...
                queryVehicleTrajectory(1,fieldsTrajectory.north), ...
                queryVehicleTrajectory(1,fieldsTrajectory.up), ...
                queryVehicleTrajectory(end,fieldsTrajectory.east), ...
                queryVehicleTrajectory(end,fieldsTrajectory.north), ...
                queryVehicleTrajectory(end,fieldsTrajectory.up));
        end % NOTE: Ends flag.verbose (Step 7)
        
        %% Step 8: Add offset and swerving to the queryVehicleTrajectory
        % fcn_swervedVehicleTrajectory takes in queryVehicleTrajectory and 
        % adds two sine waves to it based on swerve.amplitude, 
        % swerve.time_period, and fixed lateral offset to the vehicle trajectory
        if flag.lateralOffset
            % add random lateral offset
            swerve.lateralOffset = swerve.MaxLateralOffset*(2*rand-1);
        else
            swerve.lateralOffset = 0;
        end
        swervedVehicleTrajectory = ...
            fcn_VD_swervedVehicleTrajectoryST(queryVehicleTrajectory,sectionRef_table,...
            friction_ref_table,swerve.amplitude,swerve.timePeriod,swerve.lateralOffset,refLLA,...
            fieldsTrajectory,flag.swerve); % lateral offset works even flag.swerve is false
        
        % check loaded trajectory
        if flag.doDebug
            fcn_VD_plotTrajectory(swervedVehicleTrajectory(:,[12,13]),12345); % Plot output trajectory
            fcn_VD_plotTimeYaw(swervedVehicleTrajectory(:,19),...
                swervedVehicleTrajectory(:,15),12346); % Plot yaw
        end
        
        %% Step 9: Run the simulation
        inputTrajectory = swervedVehicleTrajectory;
        mdlParam.max_station = max(inputTrajectory(:,fieldsTrajectory.station));
        
        % set parameters in Trajectory query block
        mdlParam.trajectorySize = size(inputTrajectory);    % size of input trajectory
        mdlParam.flagFullWndow  = false; % parameter that defines the search window for trajectory query
        % parameter to decide the size of the search window by using the maximum
        % distance that can be traversed by a vehicle moving with velocity 40m/s in 0.1s
        mdlParam.searchDistance = 40*0.1; % [meters]
        
        % duration of the simulation
        sim_duration_from_aimsun_time = ...
            max(inputTrajectory(:,fieldsTrajectory.aimsunTime)) - ...
            min(inputTrajectory(:,fieldsTrajectory.aimsunTime));
        
        % initial conditions
        initial.east  = inputTrajectory(1,fieldsTrajectory.east); % [meters]
        initial.north = inputTrajectory(1,fieldsTrajectory.north); % [meters]
        initial.heading = inputTrajectory(1,fieldsTrajectory.yaw); % [rad]
        initial.longitudinalSpeed = inputTrajectory(1,fieldsTrajectory.speed); % [m/s]
        initial.wheelSpeeds = initial.longitudinalSpeed/vehicle.Re*ones(1,4); % [rad/s]
        set_friction = 0.9;
        initial.longitudinalSpeed = initial.longitudinalSpeed*sqrt(set_friction/0.9);
        fprintf(1,'\nSimulation is running ...\n');
        profile off
        profile on -timer 'performance'
        simlulation_start = tic;
        sim('mdl_VD_simulateFrictionUtilization.slx', round(5*sim_duration_from_aimsun_time,2));
        time_elapsed = toc(simlulation_start);
        profile viewer
        
        if flag.verbose
            fprintf(1,'Simulation is done.\n\tWall time to run the simulation: %.5f seconds. \n',time_elapsed);
            fprintf(1,'\tDuration of vehicle trajectory in Aimsun: %.5f seconds. \n',sim_duration_from_aimsun_time);
        end % NOTE: Ends flag.verbose (Step 9)
        
        % Estimate friction utilization
        longitudinal_tireForce = TireForces(:,[5,6,7,8]);
        lateral_tireForce      = TireForces(:,[1,2,3,4]);
        total_tireForce        = sqrt((longitudinal_tireForce.^2)+(lateral_tireForce.^2));
        frictionalForce        = (mean_true_friction.*Fz);
        friction_utilization   = total_tireForce./frictionalForce;
        
        h_fig = figure(09);
        set(h_fig, 'Name', 'fcn_VD_plotTimeSteeringAngle');
        width = 600; height = 400; right = 100; bottom = 400;
        set(gcf, 'position', [right, bottom, width, height])
        max_value = max(friction_utilization, [], 'all');
        min_value = min(friction_utilization, [], 'all');
        offset    = 0.1*(max_value-min_value);
        max_value_x = max(wheel_station, [], 'all');
        min_value_x = min(wheel_station, [], 'all');
        offset_x    = 0.1*(max_value_x-min_value_x);
        clf
        hold on
        plot(wheel_station(:,1), friction_utilization(:,1), 'r.', 'Markersize', 3)
        plot(wheel_station(:,2), friction_utilization(:,2), 'b.', 'Markersize', 3)
        plot(wheel_station(:,3), friction_utilization(:,3), 'g.', 'Markersize', 3)
        plot(wheel_station(:,3), friction_utilization(:,4), 'c.', 'Markersize', 3)
        grid on
        ylabel('Friction Utilization')
        xlabel('Station [m]')
        legend('FL','FR','RL','RR','Location','best')
        ylim([min_value-offset max_value+offset])
        xlim([min_value_x-offset_x max_value_x+offset_x])
        title('Friction Utilization')
        fcn_VD_plotStationFrictionUtilization(wheel_station,friction_utilization,10); % Plot friction utilization

        fcn_VD_plotTimeLongitudinalVelocity(time,States(:,1),21); % Plot longitudinal velocity
        fcn_VD_plotTimeLateralVelocity(time,States(:,2),22); % Plot lateral velocity
        fcn_VD_plotTimeYawRate(time,States(:,3),23); % Plot yaw rate
        fcn_VD_plotTimeWheelSpeed(time,States(:,(4:7)),24); % Plot wheel speed
        
        fcn_VD_plotTrajectory(pose(:,[1,2]),25); % Plot output trajectory
        fcn_VD_plotTimeYaw(time,pose(:,3),26); % Plot yaw
        %% Step 10: Plot the results
        if flag.plot
            fcn_VD_plotTimeSteeringAngle(time,delta,11); % Plot steering angle
            fcn_VD_plotTimeWheelTorque(time,wheel_torque,12); % Plot wheel torque
            
            fcn_VD_plotTimeSlipAngle(time,alpha,13); % Plot slip angles
            fcn_VD_plotTimeWheelSlip(time,kappa,14); % Plot wheel slip
            fcn_VD_plotTimeNormalForce(time,Fz,15); % Plot normla force
            
            fcn_VD_plotTimeLongitudinalTireForce(time,Fx,16); % Plot longitudinal tire force
            fcn_VD_plotTimeLateralTireForce(time,Fy,17); % Plot lateral tire force
            fcn_VD_plotTimeAligningMoment(time,Mz,18); % Plot aligning moment
            
            fcn_VD_plotTimeLongitudinalAcceleration(time,States(:,8),19); % Plot longitudinal acceleration
            fcn_VD_plotTimeLateralAcceleration(time,States(:,9),20); % Plot lateral acceleration
            
            fcn_VD_plotTimeLongitudinalVelocity(time,States(:,1),21); % Plot longitudinal velocity
            fcn_VD_plotTimeLateralVelocity(time,States(:,2),22); % Plot lateral velocity
            fcn_VD_plotTimeYawRate(time,States(:,3),23); % Plot yaw rate
            fcn_VD_plotTimeWheelSpeed(time,States(:,(4:7)),24); % Plot wheel speed
            fcn_VD_plotTrajectory(pose(:,[1,2]),25); % Plot output trajectory
            fcn_VD_plotTimeYaw(time,pose(:,3),26); % Plot yaw
            
            fcn_VD_plotTimeLateralOffsetError(time,lateral_offset_error,28);
            
            fcn_VD_plotTimeFriction(time,true_friction,29); % Plot true friction
            fcn_VD_plotTimeFriction(time,preview_true_friction,30); % Plot preview friction
        end % NOTE: Ends flag.plot
    end % NOTE: Ends for 'indexTrajectoryOfInterest'
    
end % NOTE: Ends for loop 'indexSectionOfInterest'

%% Play music at the end of code
load chirp
sound(y,Fs)
load handel
sound(y,Fs)