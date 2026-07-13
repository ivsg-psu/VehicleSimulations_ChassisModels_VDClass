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
% - Integrate the method with online path updates.
% - Validate the method experimentally on the vehicle platform.

%% Prep the workspace
close all;
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

%% Part 1 - Define demonstration scenario

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

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNumSpeedPlanner = 20001;

%% Part 2 - Compute predictive safety speed profile

[UxDesired, road, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNumSpeedPlanner);

if isempty(intersectionPoints)
    error('No intersection was detected in the demonstration scenario.');
end


%% Part 3 - Display results

fprintf('Collision detected at s = %.3f m\n', sCollision);
fprintf('Road final station: %.3f m\n', road.s(end));
fprintf('Initial planned speed: %.3f m/s\n', UxDesired(1));
fprintf('Terminal planned speed: %.3f m/s\n', UxDesired(end));
fprintf('Braking starts at s = %.3f m\n', sBrakingStart);
fprintf('Braking distance: %.3f m\n', brakingDistance);

figure(figNumSpeedPlanner);
hold on;

if ~isnan(sBrakingStart)
    xline( ...
        sBrakingStart, ...
        '--', ...
        sprintf('Braking start: %.2f m', sBrakingStart));
end

xline( ...
    sCollision, ...
    ':', ...
    sprintf('Collision: %.2f m', sCollision));

%% Part 4 - Plot path geometry

geometryFigNum = 20002;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot( ...
    predictedPath(:,1), ...
    predictedPath(:,2), ...
    '-o', ...
    'LineWidth', 2, ...
    'DisplayName', 'Predicted path');

plot( ...
    boundaryPath(:,1), ...
    boundaryPath(:,2), ...
    '-s', ...
    'LineWidth', 2, ...
    'DisplayName', 'Boundary');

plot( ...
    intersectionPoints(:,1), ...
    intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersection');

xlabel('X [m]');
ylabel('Y [m]');
title('Predicted path and detected collision');
legend('Location','best');

%% Part 5 - Test different friction values

muValues = [0.2 0.4 0.6 0.8 1.0];

sBrakingStartAll = nan(size(muValues));
brakingDistanceAll = nan(size(muValues));

comparisonFigNum = 30001;

figure(comparisonFigNum);
clf;
hold on;
grid on;

for iMu = 1:length(muValues)

    currentMu = muValues(iMu);

    [currentUxDesired, currentRoad, ~, ...
     ~, currentSBrakingStart, ...
     currentBrakingDistance, ~] = ...
        fcn_VD_computePredictiveSafetySpeedProfile( ...
            predictedPath, ...
            boundaryPath, ...
            currentMu, ...
            plannerSettings, ...
            ds, ...
            []);

    sBrakingStartAll(iMu) = currentSBrakingStart;
    brakingDistanceAll(iMu) = currentBrakingDistance;

    plot( ...
        currentRoad.s, ...
        currentUxDesired, ...
        'LineWidth', 1.5, ...
        'DisplayName', ...
        sprintf('\\mu = %.1f', currentMu));

end

xline( ...
    sCollision, ...
    ':', ...
    sprintf('Collision: %.2f m', sCollision));

xlabel('Station s [m]');
ylabel('Desired speed [m/s]');
title('Effect of friction on braking profile');
legend('Location','best');

fprintf('\nFriction comparison:\n');
fprintf('mu\tBraking start [m]\tBraking distance [m]\n');

for iMu = 1:length(muValues)

    fprintf( ...
        '%.1f\t%.3f\t\t\t%.3f\n', ...
        muValues(iMu), ...
        sBrakingStartAll(iMu), ...
        brakingDistanceAll(iMu));

end
