%% script_test_fcn_VD_computePredictiveSafetySpeedProfile
%
% This script tests the function
%
%       fcn_VD_computePredictiveSafetySpeedProfile
%
% The function evaluates a predicted vehicle path against one or multiple
% path boundaries. If the predicted path intersects a boundary, the first
% intersection along the path is interpreted as a future safety violation,
% and the existing speed planner is used to compute a braking profile that
% brings the vehicle to rest before that point.
%
% The following cases are evaluated:
%
% The following cases are evaluated:
%
% 1. A predicted path with one boundary intersection.
% 2. A predicted path with no boundary intersections.
% 3. A predicted path with two boundary intersections.
% 4. A curved predicted path, showing that the collision station is based
%    on distance along the path and not on the global x-coordinate.
% 5. A sample path generated using the Path Class Library.
% 6. A friction-coefficient sweep, showing the effect of available friction
%    on braking distance.
% 7. Several independent obstacles represented as multiple boundary paths.
% 8. A mapped pavement-boundary dataset used as a preliminary boundary test.
% 9. A case where the collision is too close for the vehicle to stop from
%    the requested initial speed.
% 10. A Gao-based oval trajectory used as an additional non-trivial
%     predicted-path test case.
% 11. Mapped Test Track boundaries represented as outer and inner limits of
%     the drivable area.
%
% The tests verify the detected collision station, terminal speed,
% braking-start station, braking distance, behavior when no collision is
% present, and behavior when multiple boundaries are provided.
%
% REVISION HISTORY:
%
% 2026_07_09 by Jaime Rodriguez
% - Created test script.
% - Added one-intersection test.
% - Added no-intersection test.
% - Added two-intersection test.
% - Added geometric visualization of the predicted path and boundary.
%
% 2026_07_10 by Jaime Rodriguez
% - Added tests for curved predicted paths.
% - Added tests for different friction coefficients.
% - Added tests for several independent obstacles.
% - Added tests using mapped Test Track boundaries.
% - Added test for infeasible stopping distance.
% - Added test using a Gao-based oval predicted path.
%
% TO-DO:
%

%% Prep the workspace
close all;
clear;
clc;


%% BASIC test 1 - One intersection
% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.6;
plannerSettings.g          = 9.81;

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
% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.6;
plannerSettings.g          = 9.81;

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
% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.6;
plannerSettings.g          = 9.81;

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

%% BASIC test 4 - Curved path with one intersection
% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.6;
plannerSettings.g          = 9.81;

theta = linspace(0, pi/2, 50)';

radius = 10;

predictedPath = [
    radius*sin(theta), ...
    radius*(1-cos(theta))
];

boundaryPath = [
    7 -2
    7 12
];

figNum = 20004;

[~, ~, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

geometryFigNum = 30004;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot(predictedPath(:,1), predictedPath(:,2), ...
    '-o', ...
    'DisplayName', 'Curved predicted path');

plot(boundaryPath(:,1), boundaryPath(:,2), ...
    '-s', ...
    'DisplayName', 'Boundary');

plot(intersectionPoints(:,1), intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersection');

xlabel('X [m]');
ylabel('Y [m]');
title('Curved path with one intersection');
legend('Location','best');

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
    sprintf('Collision: %.2f m', sCollision), ...
    'LineWidth', 2);

fprintf('\nCurved-path test:\n');
fprintf('Collision station: %.3f m\n', sCollision);
fprintf('Collision X-coordinate: %.3f m\n', intersectionPoints(1,1));
fprintf('Braking start: %.3f m\n', sBrakingStart);
fprintf('Braking distance: %.3f m\n', brakingDistance);

%% BASIC test 5 - Sample trajectory from Laps Library

% Shared test settings

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.6;
plannerSettings.g          = 9.81;

% Input

pathsArray = fcn_Path_fillSamplePaths;
predictedPath = pathsArray{1};

boundaryPath = [
    50 10
    50 100
];

figNum = 20005;

[~, ~, intersectionPoints, ...
 sCollision, sBrakingStart, ~, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

geometryFigNum = 30005;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot(predictedPath(:,1), predictedPath(:,2), ...
    '.-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Sample predicted path');

plot(boundaryPath(:,1), boundaryPath(:,2), ...
    '-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Boundary');

plot(intersectionPoints(:,1), intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersections');

xlabel('X [m]');
ylabel('Y [m]');
title('Sample path from Path Class Library');
legend('Location','best');

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
        sprintf('First collision: %.2f m', sCollision), ...
        'LineWidth', 2);
end

%% PARAMETRIC test 6 - Different friction coefficients

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

muValues = [0.2 0.4 0.6 0.8 1.0];

ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

comparisonFigNum = 20006;

figure(comparisonFigNum);
clf;
hold on;
grid on;

brakingStarts = NaN(length(muValues),1);
brakingDistances = NaN(length(muValues),1);

for iMu = 1:length(muValues)

    muValue = muValues(iMu);

    [UxDesired, road, ~, ...
     ~, sBrakingStart, brakingDistance, ...
     ~] = ...
        fcn_VD_computePredictiveSafetySpeedProfile( ...
            predictedPath, ...
            boundaryPath, ...
            muValue, ...
            plannerSettings, ...
            ds, ...
            []);

    brakingStarts(iMu) = sBrakingStart;
    brakingDistances(iMu) = brakingDistance;

    plot( ...
        road.s, ...
        UxDesired, ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('\\mu = %.1f', muValue));

end

xlabel('Station [m]');
ylabel('Speed [m/s]');
title('Effect of friction coefficient on braking profile');
legend('Location','best');

fprintf('\nFriction comparison:\n');
fprintf('mu\tBraking start [m]\tBraking distance [m]\n');

for iMu = 1:length(muValues)

    fprintf( ...
        '%.1f\t%.3f\t\t\t%.3f\n', ...
        muValues(iMu), ...
        brakingStarts(iMu), ...
        brakingDistances(iMu));

end

assert(all(diff(brakingDistances) < 0), ...
    'Braking distance should decrease as friction increases.');

assert(all(diff(brakingStarts) > 0), ...
    'Braking should start later as friction increases.');

%% BASIC test 7 - Several independent obstacles

predictedPath = [
     0  0
     5  0
    10  0
    15  0
];

boundaryPath = {
    [6 -2; 6 2]
    [10 -2; 10 2]
    [13 -2; 13 2]
};

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNum = 20007;

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

assert(size(intersectionPoints,1) == 3, ...
    'Three intersections should be detected.');

assert(abs(sCollision - 6) < 1e-6, ...
    'The first collision should occur at s = 6 m.');

assert(abs(UxDesired(end)) < 1e-6, ...
    'Terminal speed should be zero.');

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

geometryFigNum = 30007;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot(predictedPath(:,1), predictedPath(:,2), ...
    '-o', ...
    'LineWidth', 2, ...
    'DisplayName', 'Predicted path');

for iBoundary = 1:length(boundaryPath)

    currentBoundaryPath = boundaryPath{iBoundary};

    plot( ...
        currentBoundaryPath(:,1), ...
        currentBoundaryPath(:,2), ...
        '-s', ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('Obstacle %d', iBoundary));

end

plot(intersectionPoints(:,1), intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersections');

xlabel('X [m]');
ylabel('Y [m]');
title('Several independent obstacles');
legend('Location','best');

fprintf('\nSeveral-obstacle test:\n');
fprintf('Number of intersections: %d\n', ...
    size(intersectionPoints,1));
fprintf('First collision station: %.3f m\n', ...
    sCollision);
fprintf('Braking start: %.3f m\n', ...
    sBrakingStart);
fprintf('Braking distance: %.3f m\n', ...
    brakingDistance);

%% BASIC test 8 - Mapped test-track boundaries
% 1. Load pavement boundaries
script_loadPavementBoundaries

% 2. Split at NaN rows
nanRows = any(isnan(pavementBoundaries), 2);

segmentStartIndices = find([true; nanRows(1:end-1)]);
segmentEndIndices = find([nanRows(2:end); true]);

boundaryPathsLatLon = cell(length(segmentStartIndices),1);

validBoundaryCount = 0;

for iSegment = 1:length(segmentStartIndices)

    currentSegment = pavementBoundaries( ...
        segmentStartIndices(iSegment):segmentEndIndices(iSegment), :);

    % Remove NaN rows
    currentSegment = currentSegment( ...
        ~any(isnan(currentSegment),2), :);

    % Keep only valid paths with at least two points
    if size(currentSegment,1) >= 2

        validBoundaryCount = validBoundaryCount + 1;

        boundaryPathsLatLon{validBoundaryCount,1} = ...
            currentSegment;

    end

end

% Remove unused empty cells
boundaryPathsLatLon = ...
    boundaryPathsLatLon(1:validBoundaryCount);

fprintf('\nMapped test-track boundaries:\n');
fprintf('Number of independent boundary segments: %d\n', ...
    length(boundaryPathsLatLon));

geometryFigNum = 30008;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

for iBoundary = 1:length(boundaryPathsLatLon)

    currentBoundary = boundaryPathsLatLon{iBoundary};

    plot( ...
        currentBoundary(:,2), ...
        currentBoundary(:,1), ...
        '.-', ...
        'DisplayName', sprintf('Boundary %d', iBoundary));

end

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('Mapped test-track boundary segments');

% 3. Convert boundary segments from Lat/Lon to ENU

referenceLatitude  = 40.86368573;
referenceLongitude = -77.83592832;
referenceAltitude  = 344.189;

referenceEllipsoid = wgs84Ellipsoid('meter');

boundaryPathsXY = cell(size(boundaryPathsLatLon));

for iBoundary = 1:length(boundaryPathsLatLon)

    currentBoundaryLatLon = boundaryPathsLatLon{iBoundary};

    latitude  = currentBoundaryLatLon(:,1);
    longitude = currentBoundaryLatLon(:,2);
    altitude  = referenceAltitude * ones(size(latitude));

    [xEast, yNorth, ~] = geodetic2enu( ...
        latitude, ...
        longitude, ...
        altitude, ...
        referenceLatitude, ...
        referenceLongitude, ...
        referenceAltitude, ...
        referenceEllipsoid);

    boundaryPathsXY{iBoundary} = [
        xEast, ...
        yNorth
    ];

end

% 4. Plot boundary segments in ENU coordinates

geometryFigNum = 30009;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

for iBoundary = 1:length(boundaryPathsXY)

    currentBoundary = boundaryPathsXY{iBoundary};

    plot( ...
        currentBoundary(:,1), ...
        currentBoundary(:,2), ...
        '.-', ...
        'DisplayName', sprintf('Boundary %d', iBoundary));

end

xlabel('East [m]');
ylabel('North [m]');
title('Mapped test-track boundaries in ENU coordinates');

selectedBoundary = boundaryPathsXY{9};

% 5. Create a predicted path crossing the mapped boundary

predictedPath = [
    200  -170
    200  -150
    200  -130
    200  -110
    200   -90
    200   -70
    200   -50
];


% 6. Plot predicted path and selected mapped boundary

geometryFigNum = 30010;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot( ...
    selectedBoundary(:,1), ...
    selectedBoundary(:,2), ...
    '.-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Mapped boundary');

plot( ...
    predictedPath(:,1), ...
    predictedPath(:,2), ...
    '-o', ...
    'LineWidth', 2, ...
    'DisplayName', 'Predicted path');

xlabel('East [m]');
ylabel('North [m]');
title('Predicted path crossing a mapped boundary');
legend('Location','best');

% 7. Compute predictive safety speed profile

boundaryPath = selectedBoundary;

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNum = 20008;

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

% 8. Mark braking start and collision station

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
    sprintf('Collision: %.2f m', sCollision), ...
    'LineWidth', 2);

fprintf('\nMapped-boundary test:\n');
fprintf('Number of intersections: %d\n', ...
    size(intersectionPoints,1));
fprintf('Collision station: %.3f m\n', ...
    sCollision);
fprintf('Braking start: %.3f m\n', ...
    sBrakingStart);
fprintf('Braking distance: %.3f m\n', ...
    brakingDistance);
fprintf('Terminal speed: %.3f m/s\n', ...
    UxDesired(end));

%% BASIC test 9 - Collision too close to stop safely

predictedPath = [
    0 0
    1 0
    2 0
];

boundaryPath = [
    0.5 -2
    0.5  2
];

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNum = 20009;

[UxDesired, road, intersectionPoints, ...
 sCollision, ~, ~, ...
 speedProfileData] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

fprintf('\nBASIC test 9 - Collision too close to stop safely\n');
fprintf('Requested initial speed: %.3f m/s\n', ...
    plannerSettings.U_initial);
fprintf('Maximum feasible initial speed: %.3f m/s\n', ...
    UxDesired(1));
fprintf('Collision station: %.3f m\n', sCollision);
fprintf('Terminal speed: %.3f m/s\n', UxDesired(end));

tolerance = 1e-3;

isFeasible = ...
    plannerSettings.U_initial <= UxDesired(1) + tolerance;

fprintf('Is braking maneuver feasible: %d\n', isFeasible);

assert(~isempty(intersectionPoints), ...
    'An intersection should be detected.');

assert(abs(sCollision - 0.5) < tolerance, ...
    'The collision station should be 0.5 m.');

assert(abs(UxDesired(end)) < tolerance, ...
    'The terminal speed should be zero.');

assert(~isFeasible, ...
    'The braking maneuver should be identified as infeasible.');

%% BASIC test 10 - Gao-based oval trajectory

% Generate the Gao-based oval trajectory
script_VD_reproduceGaoChapter5Figs

% Save the generated centerline
gaoRoad = road;

% Close figures generated internally by the Gao reproduction script
close(21001);
close(21002);
close(21003);
close(21004);

predictedPath = [
    gaoRoad.x, ...
    gaoRoad.y
];

% Boundary placed away from the initial path point
boundaryX = 300;

boundaryPath = [
    boundaryX, min(predictedPath(:,2)) - 50
    boundaryX, max(predictedPath(:,2)) + 50
];
% Planner settings

muValue = 0.6;
ds = 1;

plannerSettings.U_posted   = 15.0;
plannerSettings.U_initial  = 15.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNum = 20010;

% Compute predictive safety speed profile

[UxDesired, plannerRoad, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

% Plot trajectory and detected intersections

geometryFigNum = 30011;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot( ...
    predictedPath(:,1), ...
    predictedPath(:,2), ...
    '-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Gao-based predicted path');

plot( ...
    boundaryPath(:,1), ...
    boundaryPath(:,2), ...
    '--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Artificial boundary');

plot( ...
    intersectionPoints(:,1), ...
    intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Intersections');

xlabel('East [m]');
ylabel('North [m]');
title('Gao-based oval trajectory and artificial boundary');
legend('Location','best');

% Mark braking start and first collision

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

% Print results

fprintf('\nGao-based oval trajectory test:\n');
fprintf('Number of detected intersections: %d\n', ...
    size(intersectionPoints,1));
fprintf('First collision station: %.3f m\n', ...
    sCollision);
fprintf('Braking start: %.3f m\n', ...
    sBrakingStart);
fprintf('Braking distance: %.3f m\n', ...
    brakingDistance);
fprintf('Terminal speed: %.3f m/s\n', ...
    UxDesired(end));

% Checks

assert(~isempty(intersectionPoints), ...
    'At least one intersection should be detected.');

assert(~isnan(sCollision), ...
    'The first collision station should be defined.');

assert(abs(UxDesired(end)) < 1e-6, ...
    'Terminal speed should be zero.');

assert(abs(plannerRoad.s(end) - sCollision) < ds + 1e-6, ...
    'The planner road should end at the first collision station.');

%% BASIC test 11 - Test Track patch boundaries
% The mapped Test Track patch provides the drivable-area boundaries.
% These boundaries are not the vehicle trajectory.
%
% The predictedPath represents a future vehicle trajectory. If this
% predicted path intersects either the outer or inner boundary, the first
% intersection is treated as a potential track-exit event and is passed to
% the speed planner as a zero-terminal-speed constraint.

currentFilePath = which('script_test_fcn_VD_computePredictiveSafetySpeedProfile.m');
[currentFolder, ~, ~] = fileparts(currentFilePath);

repoRoot = fullfile(currentFolder, '..', '..');

dataFolder = fullfile(repoRoot, 'Data');

dataFile = fullfile(dataFolder, 'testTrackBoundariesXY.mat');

fprintf('Loading Test Track boundaries from:\n%s\n', dataFile);

load(dataFile);
% Plot mapped boundaries

geometryFigNum = 30012;

figure(geometryFigNum);
clf;
hold on;
grid on;
axis equal;

plot( ...
    outerBoundaryXY(:,1), ...
    outerBoundaryXY(:,2), ...
    '.-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Outer boundary');

plot( ...
    innerBoundaryXY(:,1), ...
    innerBoundaryXY(:,2), ...
    '.-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Inner boundary');

xlabel('East [m]');
ylabel('North [m]');
title('Mapped Test Track boundaries from patch');
legend('Location','best');

% Create predicted path crossing the Test Track boundaries

predictedPath = [
    -100   0
       0   0
     100   0
     200   0
     300   0
     400   0
     500   0
     600   0
];

plot( ...
    predictedPath(:,1), ...
    predictedPath(:,2), ...
    '-o', ...
    'LineWidth', 2, ...
    'DisplayName', 'Predicted path');

% Compute predictive safety speed profile

muValue = 0.6;
ds = 0.1;

plannerSettings.U_posted   = 5.0;
plannerSettings.U_initial  = 5.0;
plannerSettings.U_terminal = 0.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

figNum = 20012;

[UxDesired, road, intersectionPoints, ...
 sCollision, sBrakingStart, brakingDistance, ...
 ~] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
        predictedPath, ...
        boundaryPath, ...
        muValue, ...
        plannerSettings, ...
        ds, ...
        figNum);

figure(geometryFigNum);
hold on;

plot( ...
    intersectionPoints(:,1), ...
    intersectionPoints(:,2), ...
    'x', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Detected intersections');

legend('Location','best');

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
        sprintf('First collision: %.2f m', sCollision), ...
        'LineWidth', 2);
end

fprintf('\nTest Track patch-boundary test:\n');
fprintf('Number of intersections: %d\n', ...
    size(intersectionPoints,1));
fprintf('First collision station: %.3f m\n', ...
    sCollision);
fprintf('Braking start: %.3f m\n', ...
    sBrakingStart);
fprintf('Braking distance: %.3f m\n', ...
    brakingDistance);
fprintf('Terminal speed: %.3f m/s\n', ...
    UxDesired(end));