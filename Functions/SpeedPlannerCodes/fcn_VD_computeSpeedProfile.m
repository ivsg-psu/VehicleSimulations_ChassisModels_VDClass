function [Ux_desired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile(road, mu, speedPlannerSpecs, varargin)
%% fcn_VD_computeSpeedProfile
%
% Computes the desired longitudinal speed profile Ux,d(s) along a road,
% given the road curvature kappa(s) and friction supply mu(s). Uses Gao's
% algorithm: (1) a pointwise cornering speed limit, (2) a backward pass for
% braking feasibility, (3) a forward pass for acceleration feasibility. The
% forward and backward passes are integrated independently with fixed step
% RK4. Then the final velocity is obtained by taking the minimum velocity
% of both forward and backward passes.
% 
% Vehicle model: point mass with a single lumped friction circle of radius
% lambda * mu * g (Gao Eq. 5.7 with quasi-steady assumptions Eq. 5.5: Uy ~
% 0, dot Uy = 0, dot r = 0, r ~ kappa Ux). No per-axle resolution, no load
% transfer, no tire model.
% 
% Governing ODE (Gao Eq. 5.11, full form with grade theta):
% 
%   d(Ux^2)/ds = +/- 2 * sqrt( (lambda*mu*g*cos(theta))^2 - (kappa*Ux^2)^2)
%                - 2*g*sin(theta)
% 
% Grade sign convention: positive theta = uphill; the -2 g sin(theta) term
% always subtracts along the direction of motion (uphill costs longitudinal
% budget, downhill adds to it) regardless of pass sign.
% 
% The "+/-" is for the forward (acceleration) pass and the backward
% (braking) pass respectively when integrating in s. The backward pass is
% implemented by flipping the kappa/mu/theta arrays and integrating in the
% reversed coordinate sigma = L - s; under this coordinate change the grade
% term acquires a "+" sign because dV/d(sigma) = -dV/ds.
% 
% Grade is read from road.theta if present; if the field is missing, the
% road is treated as flat (theta = 0) for backward compatibility.
%
% FORMAT:
% 
% [Ux_desired, speedProfileData] = ...
%   fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (figNum))
% 
% INPUTS:
% 
% road: struct with the fields
%       s [Nx1] arc length [m] (uniform spacing assumed)
%       kappa [Nx1] curvature [1/m]
% 
% mu: [Nx1] friction supply on the same s grid
% 
% speedPlannerSpecs: struct with fields
%       U_posted: posted speed limit [m/s]
%       U_initial: initial speed at s(1) [m/s]
%       U_terminal: terminal speed at s(end) [m/s]
%       lambda: friction margin parameter, in (0, 1]
%       g: gravity [m/s^2] (default 9.81 if absent)
% 
% (OPTIONAL INPUTS)
% 
% figNum: a figure number to plot results. If set to -1, skips any input
% checking or debugging, no figures will be generated, and sets up code to
% maximize speed.
%
% OUTPUTS:
%
% Ux_desired: [Nx1] desired speed profile [m/s] on the same s grid as road.s
% 
% speedProfileData: struct with fields
%       U_corner   [Nx1] pointwise cornering limit [m/s]
%       U_forward  [Nx1] forward-pass profile [m/s]
%       U_backward [Nx1] backward-pass profile [m/s]
%
% DEPENDENCIES:
%
%   fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%   See the script: script_test_fcn_VD_computeSpeedProfile for a full test
%   suite.
%
% This function was written on 2026_06_07 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu

% REVISION HISTORY:
% 2026_06_07 by Aneesh Batchu, abb6486@psu.edu
% - In fcn_VD_computeSpeedProfile
%   % * Wrote the code originally

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the figNum variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
flag_max_speed = 0;
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; %   %   Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; %   %   Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS");
    MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG     = getenv("MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    % debug_figNum = 999978;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (0 == flag_max_speed)
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(3, MAX_NARGIN);

        % Check the road AND speedPlannerSpecs STRUCT input
        road_template = struct('s', {}, 'kappa', {});
        speedPlannerSpecs_template = struct( ...
            'U_posted', {}, 'U_initial', {}, 'U_terminal', {}, 'lambda', {});

        fcn_DebugTools_checkInputsToFunctions(road, 'likestructure', road_template);
        fcn_DebugTools_checkInputsToFunctions(speedPlannerSpecs, ...
            'likestructure', speedPlannerSpecs_template);
        fcn_DebugTools_checkInputsToFunctions(mu, '1column_of_numbers');


        % mu must match road.s in length
        if(numel(mu) ~= numel(road.s))
            error('mu must be the same length as road.s');
        end
        if(numel(road.kappa) ~= numel(road.s))
            error('road.kappa must be the same length as road.s');
        end

    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to show plots
figNum = []; % Empty by default
if (0==flag_max_speed) && (MAX_NARGIN == nargin)
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        figNum = temp;
        flag_do_plots = 1;
    end
end


%% Start of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the inputs into local variables 
s = road.s(:);
kappa = road.kappa(:); % Curvature
mu = mu(:);
Nrows = numel(s);

% Grade: read road.theta if present, else default to zeros.
if isfield(road, 'theta') && ~isempty(road.theta)
    theta = road.theta(:);
    if(numel(theta) ~= Nrows)
        error('road.theta must have the same length as road.s');
    end
else
    theta = zeros(Nrows, 1);
end

U_posted   = speedPlannerSpecs.U_posted;
U_initial  = speedPlannerSpecs.U_initial;
U_terminal = speedPlannerSpecs.U_terminal;
lambda     = speedPlannerSpecs.lambda;
if isfield(speedPlannerSpecs, 'g') && ~isempty(speedPlannerSpecs.g)
    g = speedPlannerSpecs.g;
else
    g = 9.81;
end

% Step size is assumed uniform. So, using the median to be robust to small
% rounding.
ds = median(diff(s));

%% Pass 1: pointwise cornering speed limit
%
% In steady cornering with dot(Ux) = 0, the tire must still supply Fx =
% m*g*sin(theta) to hold the vehicle against gravity along the road
% tangent. The friction circle then binds on the remaining lateral budget,
% giving (Gao Eq. 5.14, generalized to include grade):
%
%   Ux_corner(s)^2 = (1/|kappa|) * sqrt( (lambda*mu*g*cos(theta))^2
%                                       - (g*sin(theta))^2 )
%
% For flat road (theta = 0), this reduces to lambda*mu*g/|kappa|. The grade
% correction to the cornering limit is O(theta^2) and negligible for
% highway grades; the first-order effect of grade lives in the ODE.
%
% For straight road (kappa = 0), the lateral term is +inf and U_posted
% binds. For now, we treat |kappa| < kappa_floor as straight to avoid
% division by zero.
% 
% Note: Modify kappa_floor based on real road data. We might pick a larger
% value whern using real data (may be 1e-6)
% 
% Mask is to find the segments with road curvature more than 0 (roads that
% are not straight). It's impportant to do this as U_corner does not
% converge when road curvature is zero (denominator cannot be zero).

% Treats near zero curvature as straight to avoid dividing by zero
kappa_floor = 1e-12; 

U_corner_lateral = inf(Nrows, 1);

% Vector of stations where curvature is nonzero.
mask = abs(kappa) > kappa_floor;  

budget_sq_static = (lambda * mu(mask) * g .* cos(theta(mask))).^2 ...
                  - (g * sin(theta(mask))).^2;

U_corner_lateral(mask) = sqrt( sqrt(max(0, budget_sq_static)) ./ abs(kappa(mask)) );

% Selects the minimum of the U_posted and U_corner_lateral
U_corner = min(U_posted, U_corner_lateral);


%% Pass 2: forward pass (acceleration)
%
% Integrate dV/ds = +2*sqrt((lambda*mu*g)^2 - (kappa*V)^2), V = Ux^2,
% from s(1) with V(1) = U_initial^2. 
% 
% At each step, clip V to U_corner(s)^2 so the integrated speed never
% exceeds the cornering limit.

U_forward = fcn_INTERNAL_RK4Integrator( ...
    kappa, mu, theta, U_corner, U_initial, ds, lambda, g, +1);

%% Pass 3: backward pass (braking)
% 
% Flip the road and friction arrays, integrate the same ODE forward in the
% reversed coordinate sigma = L - s, then flip the result.
%
% This is algebraically the same as integrating dV/ds =
% -2*sqrt((lambda*mu*g*cos(theta))^2 - (kappa*V)^2) - 2*g*sin(theta)
% backward in s from s(end) with V(end) = U_terminal^2. Under the change of
% variables sigma = L - s, dV/d(sigma) = -dV/ds, so the sqrt term keeps its
% + sign but the grade term flips sign. Passing pass_sign = -1 to the
% derivative flips the grade term correctly.

U_backward_flipped = fcn_INTERNAL_RK4Integrator( ...
    flipud(kappa), flipud(mu), flipud(theta), flipud(U_corner), U_terminal, ds, lambda, g, -1);
U_backward = flipud(U_backward_flipped);

%% Combine: elementwise minimum of forward and backward profiles
%
% The vehicle must respect both feasibility constraints, so the executable
% speed profile is the smaller of the two at every station distance (s).

Ux_desired = min(U_forward, U_backward);


%% Finding accelerations from the ODE (avoids finite-diff spikes)
%
% At each station, determine which pass binds (forward, backward, or
% cornering/posted), then evaluate the ODE-based longitudinal acceleration:
%
%   forward binding:   a_x = +sqrt((lambda*mu*g*cos(theta))^2 - (kappa*Ux^2)^2)
%                             - g*sin(theta)
%   backward binding:  a_x = -sqrt(...) - g*sin(theta)
%   cornering/posted:  a_x = 0        (Ux constant along s; steady state)
%
% Lateral: a_y = kappa * Ux^2 always (quasi-steady cornering).

tol = 1e-6;
budget_sq  = (lambda * mu .* g .* cos(theta)).^2;
lateral_sq = (kappa .* Ux_desired.^2).^2;
sqrt_term  = sqrt(max(0, budget_sq - lateral_sq));

% fwd_binds is true when Ux_desired equals U_forward (within tolerance),
% and U_forward is strictly less than both the cornering limit and U_backward
fwd_binds  = (abs(Ux_desired - U_forward)  < tol) & (Ux_desired < U_corner - tol) ...
             & (U_forward < U_backward - tol);
bwd_binds  = (abs(Ux_desired - U_backward) < tol) & (Ux_desired < U_corner - tol) ...
             & (U_backward < U_forward - tol);

a_x = zeros(Nrows, 1);
a_x(fwd_binds) = +sqrt_term(fwd_binds) - g * sin(theta(fwd_binds));
a_x(bwd_binds) = -sqrt_term(bwd_binds) - g * sin(theta(bwd_binds));

% All other points (cornering/posted binds or tied): a_x = 0
a_y = kappa .* Ux_desired.^2;

%% Create speedProfileData

speedProfileData.U_corner = U_corner;
speedProfileData.U_forward = U_forward;
speedProfileData.U_backward = U_backward;
speedProfileData.a_x = a_x;
speedProfileData.a_y = a_y;

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

    % check whether the figure already has data
    temp_h = figure(figNum);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    hold on;
    grid on;

    plot(s, U_corner,   'k--', 'LineWidth', 4, ...
        'DisplayName', 'Cornering limit');
    plot(s, U_forward,  'b-',  'LineWidth', 3, ...
        'DisplayName', 'Forward pass');
    plot(s, U_backward, 'r-',  'LineWidth', 2, ...
        'DisplayName', 'Backward pass');
    plot(s, Ux_desired, 'g-',  'LineWidth', 3, ...
        'DisplayName', 'Final U_{x,d}(s)');

    yline(U_posted, ':', sprintf('U_{posted} = %.1f m/s', U_posted), ...
        'Color', [0.4 0.4 0.4], 'DisplayName', 'U_{posted}');

    xlabel('Station (m)');
    ylabel('Speed (m/s)');
    % title('Speed planner: pointwise, forward, backward and final profile');
    legend('Location', 'best');

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U_pass = fcn_INTERNAL_RK4Integrator(...
    kappa, mu, theta, U_corner, U_initial, ds, lambda, g, pass_sign)
% 
% Integrate dV/ds = +2*sqrt((lambda*mu*g*cos(theta))^2 - (kappa*V)^2)
%                  - 2*g*sin(theta)*pass_sign
% with fixed step RK4 on the given uniform grid. (Here, V = Ux^2)
% 
% pass_sign is +1 for forward-pass integration in s and -1 for
% backward-pass integration in the flipped coordinate sigma = L - s.
% Only the grade term is sign-sensitive; the sqrt term retains its "+"
% because it represents the magnitude of the friction-limited
% longitudinal acceleration budget, which is direction-symmetric.
%
% After each step, clip V to U_corner^2 at the new station. 
% Returns U_pass = sqrt(V) on the grid.
%
% Inputs are assumed to be column vectors of the same length, with
% V_initial = U_init^2 applied at index 1.

Nrows = numel(kappa);
V = zeros(Nrows, 1);
V(1) = min(U_initial, U_corner(1))^2; % V_initial

% Runge-Kutta Integrator (RK4)
for ith_step = 1:(Nrows-1)
    % Local values at the start (a) and end (b) of the step
    kappa_a = kappa(ith_step);
    mu_a = mu(ith_step);
    theta_a = theta(ith_step);
    kappa_b = kappa(ith_step + 1);
    mu_b = mu(ith_step + 1);
    theta_b = theta(ith_step + 1);

    % Midpoint values via linear interpolation
    kappa_m = 0.5 * (kappa_a + kappa_b);
    mu_m    = 0.5 * (mu_a + mu_b);
    theta_m = 0.5 * (theta_a + theta_b);

    V_curr = V(ith_step);

    % RK4 stages on dV/ds = f(s, V)
    k1 = fcn_INTERNAL_speedODE(V_curr, kappa_a, mu_a, theta_a, lambda, g, pass_sign);
    k2 = fcn_INTERNAL_speedODE(V_curr + 0.5*ds*k1, kappa_m, mu_m, theta_m, lambda, g, pass_sign);
    k3 = fcn_INTERNAL_speedODE(V_curr + 0.5*ds*k2, kappa_m, mu_m, theta_m, lambda, g, pass_sign);
    k4 = fcn_INTERNAL_speedODE(V_curr + ds*k3, kappa_b, mu_b, theta_b, lambda, g, pass_sign);

    V_next = V_curr + (ds/6) * (k1 + 2*k2 + 2*k3 + k4);

    % Clip to the cornering limit at the new station
    V(ith_step + 1) = min(V_next, U_corner(ith_step + 1)^2);
end

% V cannot be negative. Numerical noise near saturation can produce
% small negatives, so floor at zero before sqrt.
U_pass = sqrt(max(V, 0));

end % Ends fcn_INTERNAL_RK4Integrator



function dVds = fcn_INTERNAL_speedODE(V, kappa_local, mu_local, theta_local, lambda, g, pass_sign)
% 
% Right-hand side of the speed-planner ODE (Gao Eq. 5.11, full form):
%
%   dV/ds = 2*sqrt((lambda*mu*g*cos(theta))^2 - (kappa*V)^2)
%           - 2*g*sin(theta)*pass_sign
%
% pass_sign convention: +1 for the forward pass in s (acceleration),
% -1 for the backward pass integrated in the flipped coordinate.
% Under the flipped coordinate change dV/d(sigma) = -dV/ds, the grade
% term flips sign; the sqrt (budget magnitude) does not.
% 
% The square root argument is clipped at zero so that saturation 
% (kappa*V == lambda*mu*g*cos(theta)) yields zero longitudinal 
% acceleration budget rather than complex values.

budget_sq  = (lambda * mu_local * g * cos(theta_local))^2;
lateral_sq = (kappa_local * V)^2;
dVds = 2 * sqrt(max(0, budget_sq - lateral_sq)) - 2 * g * sin(theta_local) * pass_sign;

end % Ends fcn_INTERNAL_speedODE


