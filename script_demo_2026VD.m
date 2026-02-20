% Script to demonstrate the Vehicle Dynamics (VD) library (special release
% for 2026 VD class)
%
% This is the explanation of the code that can be found by running
%
%       script_demo_2026VD.m
%
% This is a script to demonstrate the functions within the "VD"
% library. This code repo is typically located at:
%
% https://github.com/ivsg-psu/Classes_VehicleDynamics_2026VDLibrary
%
% This script was written by S. Brennan using code from Satya
% Questions or comments? sbrennan@psu.edu

% REVISION HISTORY:
%
% 2024_03_21 by Sean Brennan, sbrennan@psu.edu
% - Created master demo script by copy/pasting from ParseXODR
%
% 2024_09_26 by Sean Brennan, sbrennan@psu.edu
% - Updated function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
%
% 2025_12_30 by Sean Brennan, sbrennan@psu.edu
% - Updated formatting of Bicycle_MATLAB subdirectory functions
%   % * Edited these to standard form 
%
% 2026_01_31 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_kinematicPointMassModelRK4
%   % * Renamed function to indicate that it is for derivatives only
%   % * Improved header comments
%   % * Fixed input checking to use DebugTools
%   % * Set plot handle DisplayName for RK4 MATLAB plot.
% - In fcn_VD_kinematicPointMassModelSimulink
%   % * Renamed function to indicate that it is for derivatives only
%   % * Improved header comments
%   % * Fixed input checking to use DebugTools
% - In script_test_fcn_VD_kinematicPointMassModelSimulink
%   % * Wrote the code originally, using breakDataIntoLaps as starter
% - In script_test_fcn_VD_kinematicPointMassModelRK4
%   % * Wrote the code originally, using breakDataIntoLaps as starter
%
% 2026_02_01 by Sean Brennan, sbrennan@psu.edu
% - In fcn_VD_kinematicPointMassModelSimulink
%   % * Added a Simulink-based interface and multiple release-specific
%   %   % model files, plus a test script, and a small MATLAB plot fix.
%   % * Added fcn_VD_kinematicPointMassModelSimulink.m: Simulink-based
%   %   % simulation wrapper that runs the appropriate release SLX, extracts
%   %   % outputs, and plots results.
%   % * Added many versioned SLX files (R2018b through R2025b) under
%   %   % KinematicPointMass_Simulink and a top-level
%   %   % mdl_VD_KinematicPointMassModel_2024b.slx.
% - Moved comparison scripts out of main codes and into subfolder
%   % * Added this to path
%
% 2026_02_10 by Sean Brennan, sbrennan@psu.edu
% - In script_demo_VD
%   % * Added check if repo is ready for release
%
% (As script_demo_2026VD)
%
%
% 2026_02_20 by Sean Brennan, sbrennan@psu.edu
% - In script_demo_2026VD
%   % * Ported code over from VDClass library


% TO-DO:
% - 2025_12_29 by Sean Brennan, sbrennan@psu.edu
%   % (add items here)


%% Make sure we are running out of root directory
st = dbstack;
thisFile = which(st(1).file);
[filepath,name,ext] = fileparts(thisFile);
cd(filepath);

%%% START OF STANDARD INSTALLER CODE %%%%%%%%%

%% Clear paths and folders, if needed
if 1==1
    clear flag_Laps_Folders_Initialized
end

if 1==0
    fcn_INTERNAL_clearUtilitiesFromPathAndFolders;
end

if 1==0
    % Resets all paths to factory default
    restoredefaultpath;
end

%% Install dependencies
% Define a universal resource locator (URL) pointing to the repos of
% dependencies to install. Note that DebugTools is always installed
% automatically, first, even if not listed:
clear dependencyURLs dependencySubfolders
ith_repo = 0;

ith_repo = ith_repo+1;
dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary';
dependencySubfolders{ith_repo} = {'Functions','Data'};

ith_repo = ith_repo+1;
dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary';
dependencySubfolders{ith_repo} = {'Functions','Data'};

ith_repo = ith_repo+1;
dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/FieldDataCollection_GPSRelatedCodes_GPSClass';
dependencySubfolders{ith_repo} = {'Functions'};


% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_PathTools_GetUserInputPath';
% dependencySubfolders{ith_repo} = {''};
%
% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/FieldDataCollection_VisualizingFieldData_PlotRoad';
% dependencySubfolders{ith_repo} = {'Functions','Data'};
%
% ith_repo = ith_repo+1;
% dependencyURLs{ith_repo} = 'https://github.com/ivsg-psu/PathPlanning_MapTools_MapGenClassLibrary';
% dependencySubfolders{ith_repo} = {'Functions','testFixtures','GridMapGen'};



%% Do we need to set up the work space?
if ~exist('flag_Laps_Folders_Initialized','var')

    % Clear prior global variable flags
    clear global FLAG_*

    % Navigate to the Installer directory
    currentFolder = pwd;
    cd('Installer');
    % Create a function handle
    func_handle = @fcn_DebugTools_autoInstallRepos;

    % Return to the original directory
    cd(currentFolder);

    % Call the function to do the install
    func_handle(dependencyURLs, dependencySubfolders, (0), (-1));

    % Add this function's folders to the path
    this_project_folders = {...
        'Functions', ...
        'Data',...
        'VD_Utilities', ...
		'KinematicPointMass_MATLAB',...
        'KinematicPointMass_Simulink',...
		'KinematicBicycle_MATLAB',...
		'KinematicBicycle_Simulink',... % 'Bicycle_MATLAB', ...
        'VD_Utilities_Bicycle',...  % 'Bicycle_Simulink',...     % 'DualTrack_MATLAB',...   % 'DualTrack_Simulink',...
        'VD_plot',...
		'ComparisonScripts',...
        'Datafiles'};
    fcn_DebugTools_addSubdirectoriesToPath(pwd,this_project_folders)

    flag_Laps_Folders_Initialized = 1;
end

%%% END OF STANDARD INSTALLER CODE %%%%%%%%%

%% Set environment flags for input checking in Laps library
% These are values to set if we want to check inputs or do debugging
setenv('MATLABFLAG_VD_FLAG_CHECK_INPUTS','1');
setenv('MATLABFLAG_VD_FLAG_DO_DEBUG','0');

%% Set environment flags that define the ENU origin
% This sets the "center" of the ENU coordinate system for all plotting
% functions
% Location for Test Track base station
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LATITUDE','40.86368573');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_LONGITUDE','-77.83592832');
setenv('MATLABFLAG_PLOTROAD_REFERENCE_ALTITUDE','344.189');


%% Set environment flags for plotting
% These are values to set if we are forcing image alignment via Lat and Lon
% shifting, when doing geoplot. This is added because the geoplot images
% are very, very slightly off at the test track, which is confusing when
% plotting data
setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LAT','-0.0000008');
setenv('MATLABFLAG_PLOTROAD_ALIGNMATLABLLAPLOTTINGIMAGES_LON','0.0000054');

%% Check if repo is ready for release
if 1==0
	figNum = 999999;
	repoShortName = '_VD_';
	fcn_DebugTools_testRepoForRelease(repoShortName, (figNum));
end

%% Start of Demo Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _            __   _____                          _____          _
%  / ____| |           | |          / _| |  __ \                        / ____|        | |
% | (___ | |_ __ _ _ __| |_    ___ | |_  | |  | | ___ _ __ ___   ___   | |     ___   __| | ___
%  \___ \| __/ _` | '__| __|  / _ \|  _| | |  | |/ _ \ '_ ` _ \ / _ \  | |    / _ \ / _` |/ _ \
%  ____) | || (_| | |  | |_  | (_) | |   | |__| |  __/ | | | | | (_) | | |___| (_) | (_| |  __/
% |_____/ \__\__,_|_|   \__|  \___/|_|   |_____/ \___|_| |_| |_|\___/   \_____\___/ \__,_|\___|
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Start%20of%20Demo%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the demo code for the Vehicle Dynamics (VD) class library!')





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


%% function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
function fcn_INTERNAL_clearUtilitiesFromPathAndFolders
% Clear out the variables
clear global flag* FLAG*
clear flag*
clear path

% Clear out any path directories under Utilities
path_dirs = regexp(path,'[;]','split');
utilities_dir = fullfile(pwd,filesep,'Utilities');
for ith_dir = 1:length(path_dirs)
    utility_flag = strfind(path_dirs{ith_dir},utilities_dir);
    if ~isempty(utility_flag)
        rmpath(path_dirs{ith_dir});
    end
end

% Delete the Utilities folder, to be extra clean!
if  exist(utilities_dir,'dir')
    [status,message,message_ID] = rmdir(utilities_dir,'s');
    if 0==status
        error('Unable remove directory: %s \nReason message: %s \nand message_ID: %s\n',utilities_dir, message,message_ID);
    end
end

end % Ends fcn_INTERNAL_clearUtilitiesFromPathAndFolders

