function [UxDesired, road, intersectionPoints, ...
    sCollision, sBrakingStart, brakingDistance, ...
    speedProfileData] = ...
    fcn_VD_computePredictiveSafetySpeedProfile( ...
    predictedPath, boundaryPath, muValue, ...
    plannerSettings, ds, varargin)
% fcn_VD_computePredictiveSafetySpeedProfile
%
% Detects intersections between a predicted vehicle path and one or
% multiple boundaries and computes a predictive longitudinal speed profile.
%
% If an intersection is detected, the function selects the first
% intersection along the predicted path and imposes a zero terminal speed
% at that station.
%
% The road geometry used by the speed planner is generated from the
% predicted path, including path coordinates, heading, and curvature.
%
% In addition to the maximum safe-speed profile, the function predicts the
% actual emergency-braking speed profile from the requested initial speed.
% If the available stopping distance is insufficient, the predicted speed
% remains greater than zero at the collision station.
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
%     The path must contain at least two distinct points.
%
% boundaryPath:
%     M-by-2 matrix containing a road boundary or obstacle as [X Y], or a
%     cell array of M-by-2 matrices containing multiple independent
%     boundaries or obstacles.
%
% muValue:
%     Scalar friction coefficient applied along the complete road.
%
% plannerSettings:
%     Structure containing the settings required by
%     fcn_VD_computeSpeedProfile. The structure must include:
%
%     U_posted:
%         Posted or maximum permitted longitudinal speed.
%
%     U_initial:
%         Requested initial vehicle speed.
%
%     U_terminal:
%         Requested terminal speed. When a collision is detected, this
%         value is internally replaced by zero.
%
%     lambda:
%         Friction-utilization margin.
%
%     g:
%         Gravitational acceleration.
%
% ds:
%     Spatial discretization interval in meters.
%
% OPTIONAL INPUT:
%
% figNum:
%     Figure number passed to fcn_VD_computeSpeedProfile. Use [] to disable
%     internal speed-planner plotting.
%
% OUTPUTS:
%
% UxDesired:
%     Maximum safe longitudinal-speed profile. When a collision is
%     detected, the profile reaches zero at the first collision station.
%
% road:
%     Road structure used by the speed planner. The structure includes:
%
%     s:
%         Station along the predicted path.
%
%     x, y:
%         Predicted-path coordinates interpolated onto the planner grid.
%
%     kappa:
%         Signed horizontal path curvature.
%
%     theta:
%         Longitudinal road grade. The current implementation assumes a
%         flat road and sets theta to zero.
%
%     psi:
%         Horizontal path heading.
%
% intersectionPoints:
%     Coordinates of all detected intersections between predictedPath and
%     the supplied boundaries.
%
% sCollision:
%     Station of the first intersection along predictedPath. Returns NaN
%     when no intersection is detected.
%
% sBrakingStart:
%     Station where the collision-related terminal constraint first reduces
%     the speed below the non-obstacle speed limit. Returns NaN when no
%     collision-related braking is detected.
%
% brakingDistance:
%     Distance from sBrakingStart to sCollision. Returns NaN when no
%     collision-related braking is detected.
%
% speedProfileData:
%     Structure containing the outputs returned by
%     fcn_VD_computeSpeedProfile and the following additional fields:
%
%     U_actual_predicted:
%         Predicted vehicle-speed profile obtained by applying the maximum
%         available braking effort from the requested initial speed.
%         Returns NaN values when no collision is detected.
%
%     U_impact_predicted:
%         Predicted residual speed at the collision station. Returns NaN
%         when no collision is detected.
%
%     is_stopping_feasible:
%         Logical value indicating whether maximum braking reduces the
%         predicted vehicle speed to approximately zero at the collision
%         station. Returns NaN when no collision is detected.
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
% 2026_07_18 by Jaime Rodriguez
% - Added predicted-path coordinates, heading, and curvature to the road
%   representation.
% - Added support for short predicted paths with at least two distinct
%   points.
% - Added prediction of the actual emergency-braking speed profile.
% - Added residual impact-speed and stopping-feasibility outputs.
% - Updated braking-start detection to distinguish obstacle-related braking
%   from curvature-related speed limitations.
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

%% Replace straight-road geometry with predicted-path geometry
%
% The predicted path may contain either a small number of points, as in
% the basic tests, or a densely sampled curved trajectory.
%
% Straight and short paths are handled without requiring an arbitrary
% minimum number of samples. Smoothing is applied only when enough path
% points are available.

% Calculate accumulated station along the predicted path

pathSegmentVectors = diff( ...
    predictedPath, ...
    1, ...
    1);

pathSegmentLengths = sqrt(sum( ...
    pathSegmentVectors.^2, ...
    2));

predictedPathStation = [
    0
    cumsum(pathSegmentLengths)
];

% Remove repeated consecutive points

uniquePointMask = [
    true
    diff(predictedPathStation) > 1e-6
];

uniquePathStation = ...
    predictedPathStation(uniquePointMask);

uniquePredictedPath = ...
    predictedPath(uniquePointMask,:);

numberOfUniquePathPoints = ...
    numel(uniquePathStation);

assert(numberOfUniquePathPoints >= 2, ...
    ['The predicted path must contain at least two ', ...
     'distinct points.']);

% Handle a two-point path as a straight road

if numberOfUniquePathPoints == 2

    road.x = interp1( ...
        uniquePathStation, ...
        uniquePredictedPath(:,1), ...
        road.s, ...
        'linear', ...
        'extrap');

    road.y = interp1( ...
        uniquePathStation, ...
        uniquePredictedPath(:,2), ...
        road.s, ...
        'linear', ...
        'extrap');

    constantRoadHeading = atan2( ...
        uniquePredictedPath(2,2) - ...
        uniquePredictedPath(1,2), ...
        uniquePredictedPath(2,1) - ...
        uniquePredictedPath(1,1));

    road.kappa = ...
        zeros(size(road.s));

    road.theta = ...
        zeros(size(road.s));

    road.psi = ...
        constantRoadHeading * ...
        ones(size(road.s));

else

    % Smooth the original path only when sufficient points exist

    pathSmoothingWindow = min( ...
        7, ...
        numberOfUniquePathPoints);

    if mod(pathSmoothingWindow,2) == 0

        pathSmoothingWindow = ...
            pathSmoothingWindow - 1;

    end

    if pathSmoothingWindow >= 3

        smoothedPathX = smoothdata( ...
            uniquePredictedPath(:,1), ...
            'movmean', ...
            pathSmoothingWindow);

        smoothedPathY = smoothdata( ...
            uniquePredictedPath(:,2), ...
            'movmean', ...
            pathSmoothingWindow);

    else

        smoothedPathX = ...
            uniquePredictedPath(:,1);

        smoothedPathY = ...
            uniquePredictedPath(:,2);

    end

    % Calculate first and second derivatives with respect to station

    dx_ds = gradient( ...
        smoothedPathX, ...
        uniquePathStation);

    dy_ds = gradient( ...
        smoothedPathY, ...
        uniquePathStation);

    d2x_ds2 = gradient( ...
        dx_ds, ...
        uniquePathStation);

    d2y_ds2 = gradient( ...
        dy_ds, ...
        uniquePathStation);

    % Calculate signed planar curvature

    curvatureDenominator = ...
        (dx_ds.^2 + dy_ds.^2).^(3/2);

    curvatureDenominator = max( ...
        curvatureDenominator, ...
        1e-9);

    predictedPathCurvature = ...
        (dx_ds .* d2y_ds2 - ...
         dy_ds .* d2x_ds2) ./ ...
        curvatureDenominator;

    % Smooth curvature only when sufficient samples exist

    curvatureSmoothingWindow = min( ...
        9, ...
        numberOfUniquePathPoints);

    if mod(curvatureSmoothingWindow,2) == 0

        curvatureSmoothingWindow = ...
            curvatureSmoothingWindow - 1;

    end

    if curvatureSmoothingWindow >= 3

        predictedPathCurvature = smoothdata( ...
            predictedPathCurvature, ...
            'movmedian', ...
            curvatureSmoothingWindow);

        predictedPathCurvature = smoothdata( ...
            predictedPathCurvature, ...
            'movmean', ...
            curvatureSmoothingWindow);

    end

    % Interpolate coordinates and curvature onto the planner station grid

    road.x = interp1( ...
        uniquePathStation, ...
        smoothedPathX, ...
        road.s, ...
        'pchip', ...
        'extrap');

    road.y = interp1( ...
        uniquePathStation, ...
        smoothedPathY, ...
        road.s, ...
        'pchip', ...
        'extrap');

    road.kappa = interp1( ...
        uniquePathStation, ...
        predictedPathCurvature, ...
        road.s, ...
        'pchip', ...
        'extrap');

    % Calculate horizontal road heading

    roadDxDs = gradient( ...
        road.x, ...
        road.s);

    roadDyDs = gradient( ...
        road.y, ...
        road.s);

    roadHeading = unwrap( ...
        atan2(roadDyDs, roadDxDs));

    % road.theta represents longitudinal road grade.
    % All current predictive-safety tests assume a flat road.

    road.theta = ...
        zeros(size(road.s));

    % road.psi represents the horizontal path heading.

    road.psi = ...
        roadHeading;

end

%% Compute speed profile

[UxDesired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile( ...
        road, ...
        mu, ...
        localPlannerSettings, ...
        figNum);

%% Predict the actual emergency-braking speed profile
%
% UxDesired represents the maximum safe-speed profile that would allow the
% vehicle to stop at the obstacle.
%
% U_actual_predicted starts from the actual requested initial speed and
% predicts the result of applying the maximum available braking effort.
%
% If the available distance is insufficient, the predicted speed remains
% greater than zero at the obstacle.

speedProfileData.U_actual_predicted = ...
    nan(size(road.s));

speedProfileData.U_impact_predicted = ...
    NaN;

speedProfileData.is_stopping_feasible = ...
    NaN;

if ~isnan(sCollision)

    numberOfRoadStations = ...
        numel(road.s);

    actualSpeedSquared = ...
        zeros(numberOfRoadStations,1);

    actualSpeedSquared(1) = ...
        plannerSettings.U_initial^2;

    curvatureVector = ...
        road.kappa(:);

    frictionVector = ...
        mu(:);

    if isfield(road, 'theta') && ...
            ~isempty(road.theta)

        gradeVector = ...
            road.theta(:);

    else

        gradeVector = ...
            zeros(numberOfRoadStations,1);

    end

    lambdaValue = ...
        plannerSettings.lambda;

    if isfield(plannerSettings, 'g') && ...
            ~isempty(plannerSettings.g)

        gravityValue = ...
            plannerSettings.g;

    else

        gravityValue = ...
            9.81;

    end

    % Direct forward integration of maximum braking:
    %
    % d(U^2)/ds =
    % -2*sqrt((lambda*mu*g*cos(theta))^2 - (kappa*U^2)^2)
    % -2*g*sin(theta)

    brakingDerivative = ...
        @(speedSquared, curvatureLocal, ...
          frictionLocal, gradeLocal) ...
        -2 * sqrt(max( ...
            0, ...
            (lambdaValue * frictionLocal * ...
             gravityValue * cos(gradeLocal))^2 - ...
            (curvatureLocal * ...
             max(speedSquared,0))^2)) ...
        -2 * gravityValue * sin(gradeLocal);

    for iStation = 1:(numberOfRoadStations-1)

        currentStepSize = ...
            road.s(iStation+1) - ...
            road.s(iStation);

        curvatureStart = ...
            curvatureVector(iStation);

        curvatureEnd = ...
            curvatureVector(iStation+1);

        curvatureMid = ...
            0.5 * (curvatureStart + curvatureEnd);

        frictionStart = ...
            frictionVector(iStation);

        frictionEnd = ...
            frictionVector(iStation+1);

        frictionMid = ...
            0.5 * (frictionStart + frictionEnd);

        gradeStart = ...
            gradeVector(iStation);

        gradeEnd = ...
            gradeVector(iStation+1);

        gradeMid = ...
            0.5 * (gradeStart + gradeEnd);

        currentSpeedSquared = ...
            actualSpeedSquared(iStation);

        k1 = brakingDerivative( ...
            currentSpeedSquared, ...
            curvatureStart, ...
            frictionStart, ...
            gradeStart);

        k2 = brakingDerivative( ...
            currentSpeedSquared + ...
            0.5 * currentStepSize * k1, ...
            curvatureMid, ...
            frictionMid, ...
            gradeMid);

        k3 = brakingDerivative( ...
            currentSpeedSquared + ...
            0.5 * currentStepSize * k2, ...
            curvatureMid, ...
            frictionMid, ...
            gradeMid);

        k4 = brakingDerivative( ...
            currentSpeedSquared + ...
            currentStepSize * k3, ...
            curvatureEnd, ...
            frictionEnd, ...
            gradeEnd);

        nextSpeedSquared = ...
            currentSpeedSquared + ...
            (currentStepSize / 6) * ...
            (k1 + 2*k2 + 2*k3 + k4);

        actualSpeedSquared(iStation+1) = ...
            max(0, nextSpeedSquared);

    end

    actualSpeedPredicted = ...
        sqrt(actualSpeedSquared);

    predictedImpactSpeed = ...
        actualSpeedPredicted(end);

    feasibilityTolerance = ...
        1e-3;

    speedProfileData.U_actual_predicted = ...
        actualSpeedPredicted;

    speedProfileData.U_impact_predicted = ...
        predictedImpactSpeed;

    speedProfileData.is_stopping_feasible = ...
        predictedImpactSpeed <= ...
        feasibilityTolerance;

end


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
