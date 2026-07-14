function [UxDesired, road, intersectionPoints, ...
    sCollision, sBrakingStart, brakingDistance, ...
    speedProfileData] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
    predictedPath, boundaryPath, muValue, ...
    plannerSettings, ds, varargin)

% fcn_VD_computePredictiveSafetySpeedProfile
%
% Detects intersections between a predicted vehicle path and one or
% multiple boundaries and computes a longitudinal speed profile.
%
% If an intersection is detected, the function selects the first
% intersection along the predicted path and imposes a zero terminal speed
% at that station.
%
% If no intersection is detected, the speed profile is calculated along
% the complete predicted path without imposing a collision-related stop.
%
% FORMAT:
%
% [UxDesired, road, intersectionPoints, ...
%  sCollision, sBrakingStart, brakingDistance, ...
%  speedProfileData] = ...
%  fcn_VD_computePredictiveSafetySpeedProfile( ...
%  predictedPath, boundaryPath, muValue, ...
%  plannerSettings, ds, (figNum))
%
% INPUTS:
%
% predictedPath:
%     N-by-2 matrix containing the predicted vehicle path as [X Y].
%
% boundaryPath:
%     M-by-2 matrix containing a road boundary or obstacle as [X Y], or a
%     cell array of M-by-2 matrices containing multiple boundaries or
%     obstacles.
%
% muValue:
%     Scalar friction coefficient.
%
% plannerSettings:
%     Structure containing the settings required by
%     fcn_VD_computeSpeedProfile.
%
% ds:
%     Spatial discretization interval in meters.
%
% OPTIONAL INPUT:
%
% figNum:
%     Figure number passed to fcn_VD_computeSpeedProfile.
%
% OUTPUTS:
%
% UxDesired:
%     Desired longitudinal speed profile.
%
% road:
%     Road structure used by the speed planner.
%
% intersectionPoints:
%     Coordinates of all detected intersections.
%
% sCollision:
%     Station of the first intersection along predictedPath. Returns NaN
%     when no intersection is detected.
%
% sBrakingStart:
%     Station where collision-related braking begins. Returns NaN when no
%     collision is detected.
%
% brakingDistance:
%     Distance from the start of braking to the first collision.
%
% speedProfileData:
%     Additional outputs returned by fcn_VD_computeSpeedProfile.
%
% DEPENDENCIES:
%
% fcn_Path_findIntersectionsBetweenPaths
% fcn_VD_buildSampleRoad
% fcn_VD_computeSpeedProfile
%
% REVISION HISTORY:
%
% 2026_07_10 by Jaime Rodriguez
% - Created function.
% - Added first-collision detection.
% - Added collision-related terminal-speed constraint.
% - Added normal speed planning when no collision is detected.
% - Added braking-start and braking-distance calculations.
%
% TO-DO:
%

%% Optional figure input

figNum = [];

if nargin >= 6
    figNum = varargin{1};
end

%% Detect intersections

% Allow a single boundary path or several independent boundary paths
if iscell(boundaryPath)
    boundaryPaths = boundaryPath;
else
    boundaryPaths = {boundaryPath};
end

numberOfBoundaries = numel(boundaryPaths);

intersectionPointsCell = cell(numberOfBoundaries,1);
sCoordinatesCell = cell(numberOfBoundaries,1);

for ithBoundary = 1:numberOfBoundaries

    currentBoundaryPath = boundaryPaths{ithBoundary};

    [currentIntersectionPoints, ...
     currentSCoordinatesInPredictedPath, ...
     ~] = ...
        fcn_Path_findIntersectionsBetweenPaths( ...
            predictedPath, ...
            currentBoundaryPath, ...
            []);

    intersectionPointsCell{ithBoundary} = ...
        currentIntersectionPoints;

    sCoordinatesCell{ithBoundary} = ...
        currentSCoordinatesInPredictedPath;

end

intersectionPoints = ...
    vertcat(intersectionPointsCell{:});

sCoordinatesInPredictedPath = ...
    vertcat(sCoordinatesCell{:});
%% Determine road length and terminal condition

if isempty(intersectionPoints)

    % No collision: use the complete predicted-path length
    pathSegments = diff(predictedPath, 1, 1);

    pathSegmentLengths = sqrt( ...
        pathSegments(:,1).^2 + pathSegments(:,2).^2);

    roadLength = sum(pathSegmentLengths);

    sCollision = NaN;

    % The vehicle does not need to stop at the end of the path
    localPlannerSettings = plannerSettings;
    localPlannerSettings.U_terminal = plannerSettings.U_posted;

else

    % Collision detected: stop at the first collision
    sCollision = min(sCoordinatesInPredictedPath);

    roadLength = sCollision;

    localPlannerSettings = plannerSettings;
    localPlannerSettings.U_terminal = 0;

end

%% Build road ending at collision station

clear roadSpec

roadSpec(1).type   = 'straight';
roadSpec(1).length = roadLength;
roadSpec(1).grade  = 0;
roadSpec(1).angle  = [];
roadSpec(1).radius = [];

clear frictionSpec

frictionSpec(1).mu_default = muValue;
frictionSpec(1).s_start    = [];
frictionSpec(1).s_end      = [];
frictionSpec(1).mu         = [];

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
        localPlannerSettings, ...
        figNum);

%% Find braking start station

if isnan(sCollision)

    % No collision-related braking
    sBrakingStart = NaN;
    brakingDistance = NaN;

else

    tolerance = 1e-3;

    UForward = speedProfileData.U_forward;

    brakingIndex = find( ...
        UxDesired < UForward - tolerance, ...
        1, ...
        'first');

    if isempty(brakingIndex)

        sBrakingStart = NaN;
        brakingDistance = NaN;

    else

        sBrakingStart = road.s(brakingIndex);
        brakingDistance = sCollision - sBrakingStart;

    end

end
end
