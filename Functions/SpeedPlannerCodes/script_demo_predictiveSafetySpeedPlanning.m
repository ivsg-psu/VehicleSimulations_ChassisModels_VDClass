% script_demo_predictiveSafetySpeedPlanning
%
% This script demonstrates predictive safety speed planning using a future
% vehicle path and a road boundary.
%
% The script detects whether the predicted vehicle path intersects a
% boundary. If a collision is predicted, the first collision station is
% passed to the speed planner as a zero-terminal-speed constraint.
%
% The resulting backward pass determines the station where braking must
% begin so that the vehicle reaches zero speed at the collision location.
%
% The script also compares the effect of different friction coefficients
% on braking distance.
%
% Run script_demo_2026VDLibrary before running this script to initialize
% the Vehicle Dynamics library and its dependencies.
%
% REVISION HISTORY:
%
% 2026_07_10 by Jaime Rodriguez
% - Created predictive safety speed-planning demonstration.
% - Added path-intersection detection.
% - Added collision-based speed-profile calculation.
% - Added braking-start visualization.
% - Added friction-coefficient comparison.
%
% TO-DO:
%
% - Replace the straight path with measured or predicted vehicle paths.
% - Add several obstacles or road boundaries.
% - Integrate the method with online path updates.
% - Validate the method experimentally on the vehicle platform.

%% Prep the workspaceclose all;
clear;
clc;

%% Start of Demo Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____ _             _            __   _____
%  / ____| |           | |          / _| |  __ \
% | (___ | |_ __ _ _ __| |_    ___ | |_  | |  | | ___ _ __ ___   ___
%  \___ \| __/ _` | '__| __|  / _ \|  _| | |  | |/ _ \ '_ ` _ \ / _ \
%  ____) | || (_| | |  | |_  | (_) | |   | |__| |  __/ | | | | | (_) |
% |_____/ \__\__,_|_|   \__|  \___/|_|   |_____/ \___|_| |_| |_|\___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Predictive safety speed-planning demonstration')

%% Part 1 - Path intersection

predictedPath = [
     0  0
     5  0
    10  0
    15  0
];

boundaryPath = [
    8 -3
    8  3
];

[intersectionPoints, ...
 sCoordinatesInPredictedPath, ...
 sCoordinatesInBoundary] = ...
    fcn_Path_findIntersectionsBetweenPaths( ...
        predictedPath, ...
        boundaryPath, ...
        10001);

%% Part 2 - Collision station

if isempty(intersectionPoints)
    error('No intersection was detected.');
end

sCollision = min(sCoordinatesInPredictedPath);

fprintf('Collision detected at s = %.3f m\n', sCollision);

%% Part 3 - Build a simple road ending at the collision point

clear roadSpec

roadSpec(1).type   = 'straight';
roadSpec(1).length = sCollision;
roadSpec(1).grade  = 0;
roadSpec(1).angle  = [];
roadSpec(1).radius = [];

clear frictionSpec

frictionSpec(1).mu_default = 0.6;
frictionSpec(1).s_start    = [];
frictionSpec(1).s_end      = [];
frictionSpec(1).mu         = [];

ds = 0.1;

[road, mu] = ...
    fcn_VD_buildSampleRoad( ...
        roadSpec, ...
        frictionSpec, ...
        ds, ...
        []);

fprintf('Road final station: %.3f m\n', road.s(end));

%% Part 4 - Common Planner settings

clear plannerSettings

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

tolerance = 1e-3;

%% Part 5 - Compute speed profile

figNumSpeedPlanner = 20001;

[UxDesired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile( ...
        road, ...
        mu, ...
        plannerSettings, ...
        figNumSpeedPlanner);

fprintf('Initial planned speed: %.3f m/s\n', UxDesired(1));
fprintf('Terminal planned speed: %.3f m/s\n', UxDesired(end));

%% Find braking start station

UForward = speedProfileData.U_forward;

brakingIndex = find( ...
    UxDesired < UForward - tolerance, ...
    1, ...
    'first');

if isempty(brakingIndex)
    sBrakingStart = NaN;
    brakingDistance = NaN;
    fprintf('No braking phase was detected.\n');
else
    sBrakingStart = road.s(brakingIndex);
    brakingDistance = sCollision - sBrakingStart;

    fprintf('Braking starts at s = %.3f m\n', sBrakingStart);
    fprintf('Braking distance: %.3f m\n', brakingDistance);
end

figure(figNumSpeedPlanner);
hold on;

if ~isnan(sBrakingStart)
    xline(sBrakingStart, '--', ...
        sprintf('Braking start: %.2f m', sBrakingStart));
end

xline(sCollision, ':', ...
    sprintf('Collision: %.2f m', sCollision));

%% Test different friction values

muValues = [0.2 0.4 0.6 0.8 1.0];

sBrakingStartAll = nan(size(muValues));
brakingDistanceAll = nan(size(muValues));

figure(30001);
clf;
hold on;
grid on;

for i = 1:length(muValues)

    %% Friction definition

    clear frictionSpec

    frictionSpec(1).mu_default = muValues(i);
    frictionSpec(1).s_start    = [];
    frictionSpec(1).s_end      = [];
    frictionSpec(1).mu         = [];

    %% Build road

    [road, mu] = ...
        fcn_VD_buildSampleRoad( ...
            roadSpec, ...
            frictionSpec, ...
            ds, ...
            []);

    %% Compute speed profile

    [UxDesired, speedProfileData] = ...
        fcn_VD_computeSpeedProfile( ...
            road, ...
            mu, ...
            plannerSettings, ...
            []);

    %% Find braking start

    UForward = speedProfileData.U_forward;

    brakingIndex = find( ...
        UxDesired < UForward - tolerance, ...
        1, ...
        'first');

    if ~isempty(brakingIndex)

        sBrakingStartAll(i) = road.s(brakingIndex);

        brakingDistanceAll(i) = ...
            sCollision - sBrakingStartAll(i);

    end

    %% Plot speed profile

    plot( ...
        road.s, ...
        UxDesired, ...
        'LineWidth', 1.5, ...
        'DisplayName', ...
        sprintf('\\mu = %.1f', muValues(i)));

end

%% Final plot settings

xline(sCollision, ':', ...
    sprintf('Collision: %.2f m', sCollision));

xlabel('Station s [m]');
ylabel('Desired speed [m/s]');
title('Effect of friction on braking profile');
legend('Location','best');

