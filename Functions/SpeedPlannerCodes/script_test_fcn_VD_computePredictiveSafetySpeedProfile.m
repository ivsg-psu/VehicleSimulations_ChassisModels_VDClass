%% script_test_fcn_VD_computePredictiveSafetySpeedProfile
%
% This script tests the function
%
%       fcn_VD_computePredictiveSafetySpeedProfile
%
% The following cases are evaluated:
%
% 1. A predicted path with one intersection.
% 2. A predicted path with no intersections.
% 3. A predicted path with two intersections.
%
% The tests verify the detected collision station, terminal speed,
% braking-start station, and behavior when no collision is present.
%
% REVISION HISTORY:
%
% 2026_07_10 by Jaime Rodriguez
% - Created test script.
% - Added one-intersection test.
% - Added no-intersection test.
% - Added two-intersection test.
% - Added geometric visualization of the predicted path and boundary.
%
% TO-DO:
%
% - Add tests for curved predicted paths.
% - Add tests for different friction coefficients.
% - Add tests for several independent obstacles.
%% Prep the workspace
close all;
clear;
clc;

%% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

%% BASIC test 1 - One intersection
% Inputs

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

figNum = 20001;

% Call function

[UxDesired, ~, ~, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

% Print results

fprintf('Collision station: %.3f m\n', sCollision);
fprintf('Braking start: %.3f m\n', sBrakingStart);
fprintf('Braking distance: %.3f m\n', brakingDistance);
fprintf('Terminal speed: %.3f m/s\n', UxDesired(end));

% Mark braking start and collision station

figure(figNum);
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

%% BASIC test 2 - No intersection
% Inputs

predictedPath = [
     0  0
     5  0
    10  0
    15  0
];

boundaryPath = [
    8  1
    8  5
];

figNum = 20002;

% Call function

[UxDesired, ~, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

% Print results

assert(isempty(intersectionPoints), ...
    'No intersections should be detected.');

assert(isnan(sCollision), ...
    'Collision station should be NaN.');

assert(isnan(sBrakingStart), ...
    'Braking start should be NaN.');

assert(isnan(brakingDistance), ...
    'Braking distance should be NaN.');

assert(~isempty(UxDesired), ...
    'A normal speed profile should still be calculated.');

assert(abs(UxDesired(end) - plannerSettings.U_posted) < 1e-6, ...
    'The vehicle should continue at the posted speed.');

% Mark braking start and collision station

figure(figNum);
hold on;

if ~isnan(sBrakingStart)
    xline( ...
        sBrakingStart, ...
        '--', ...
        sprintf('Braking start: %.2f m', sBrakingStart));
end

if ~isnan(sCollision)
    xline( ...
        sCollision, ...
        ':', ...
        sprintf('Collision: %.2f m', sCollision));
end

fprintf('\nNo-intersection test:\n');
fprintf('No collision detected.\n');
fprintf('The vehicle continues along the full predicted path.\n');
fprintf('Terminal speed: %.3f m/s\n', UxDesired(end));

%% BASIC test 3 - Two intersections
% Inputs

predictedPath = [
     0  0
     5  0
    10  0
    15  0
];

boundaryPath = [
     4 -3
     4  3
    11  3
    11 -3
];

figNum = 20003;

% Call function

[UxDesired, road, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 speedProfileData] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

% Find which detected intersection is the first collision
[~, firstCollisionIndex] = ...
    min(abs(intersectionPoints(:,1) - sCollision));

% Plot braking profile result

figure(figNum);
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
    sprintf('First collision: %.2f m', sCollision), ...
    'LineWidth', 2);

% BASIC test 3.1 - Plot path and both intersections

geometryFigNum = 30003;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot(predictedPath(:,1), predictedPath(:,2), ...
    '-o', ...
    'DisplayName', 'Predicted path');

plot(boundaryPath(:,1), boundaryPath(:,2), ...
    '-s', ...
    'DisplayName', 'Boundary');

plot(intersectionPoints(:,1), intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersections');

plot( ...
    intersectionPoints(firstCollisionIndex,1), ...
    intersectionPoints(firstCollisionIndex,2), ...
    'o', ...
    'MarkerSize', 10, ...
    'LineWidth', 2, ...
    'DisplayName', 'First collision');

xlabel('X [m]');
ylabel('Y [m]');
title('Predicted path and detected intersections');
legend('Location','best');

% Checks

assert(size(intersectionPoints,1) == 2, ...
    'Two intersections should be detected.');

assert(abs(sCollision - 4) < 1e-6, ...
    'The first collision should occur at s = 4 m.');

assert(abs(UxDesired(end)) < 1e-6, ...
    'Terminal speed should be zero.');

fprintf('\nTwo-intersection test:\n');
fprintf('Number of intersections: %d\n', ...
    size(intersectionPoints,1));
fprintf('First collision station: %.3f m\n', ...
    sCollision);
fprintf('Braking start: %.3f m\n', ...
    sBrakingStart);
fprintf('Braking distance: %.3f m\n', ...
    brakingDistance);