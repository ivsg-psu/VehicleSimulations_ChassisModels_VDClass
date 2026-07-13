function [road, mu] = fcn_VD_buildSampleRoad(road_specifications, friction_specifications, ds, varargin)
%% fcn_VD_buildSampleRoad
% 
% Builds a road geometry and a friction-supply profile on a uniform
% arc-length grid. Designed to be a single source of (s, kappa, mu)
% for the downstream speed planner and MPC pipeline.
% 
% FORMAT:
%
%   [road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds)
%
% INPUTS:
%
% road_specification: a struct array describing road segments in order. 
% 
%       type: 'straight' or 'arc'
%       length: segment arc length [m]   (for 'straight')
%       radius: signed radius [m]        (for 'arc'; +ve = left turn)
%       angle: sweep angle [rad]        (for 'arc')
%       grade: road grade theta [rad] over the segment.
%              Positive theta = uphill. Defaults to 0 if not given.
%              Applied as a constant grade over the segment
%     For an 'arc', the segment length is radius*abs(angle) and is
%     computed internally. Only two of {radius, angle, length} need
%     to be supplied for an arc, but for clarity user can give both.
% 
% friction_specifications: a struct array describing the friction supply. 
% 
%   s_start: start arc length of segment [m]
%   s_end: end arc length of segment [m]
%   mu: friction coefficient in segment
%   L_transition: (optional) length [m] of a linear ramp at
%                 each end of the patch, smoothing the
%                 transition between the patch mu and the
%                 surrounding background mu.
% Segments are applied in order. Any (s) not covered by any spec element
% receives the default mu (last field of frictionSpec(1) labeled
% 'mu_default' if present, else
% 0.9 hard-coded fallback).
% 
% ds: arc-length discretization step [m]
%
% OUTPUTS:
%
% road: struct with fields
%       s [Nx1] arc length [m]
%       kappa [Nx1] curvature [1/m]
%       theta [Nx1] grade [rad]
%       psi [Nx1] heading [rad]
%       x [Nx1] global x [m]
%       y [Nx1] global y [m]
%       segments: a struct array of segment boundaries
% 
% mu: [Nx1] friction supply on the same s grid
%
% DEPENDENCIES:
%
%  fcn_DebugTools_checkInputsToFunctions    
% 
% EXAMPLES:
%
% See the script: script_test_fcn_VD_buildSampleRoad for a full test suite.
%
% This function was written on 2026_06_02 by Aneesh Batchu
% Questions or comments? abb6486@psu.edu 

% REVISION HISTORY:
% 2026_06_02 by Aneesh Batchu, abb6486@psu.edu
% - In fcn_VD_buildSampleRoad
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
    MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG = getenv("MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VEHICLEDYNAMICS_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VEHICLEDYNAMICS_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_figNum = 999978; %#ok<NASGU>
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


if (0==flag_max_speed)
    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        narginchk(3,MAX_NARGIN);

        road_specifications_template = struct('type',{},'length',{},'radius',{},'angle',{},'grade',{});
        friction_specifications_template = struct('s_start',{},'s_end',{},'mu',{});

    
        % Check the road_specifications input 
        fcn_DebugTools_checkInputsToFunctions(...
            road_specifications, 'likestructure',road_specifications_template);
        
        % Check the friction_specifications input 
        fcn_DebugTools_checkInputsToFunctions(...
            friction_specifications, 'likestructure',friction_specifications_template);

        %    % Check the ds input 
        % fcn_DebugTools_checkInputsToFunctions(...
        %     ds, '1column_of_numbers',[1 1]);



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
%% Use the road specification to build (s, kappa) 
% 
% For each segment, compute its arc length and curvature, append a uniform
% sub set to the global set. Within each segment we use exactly round(L/ds)
% sub steps so the local step is L/round(L/ds), which equals ds to within
% rounding

if isempty(ds)
    ds = 1; % segement size
end

s_all = 0;
kappa_all = 0;
theta_all = 0;
segBounds = zeros(numel(road_specifications), 2);   % [s_start s_end] per segment

s_current = 0;

for ith_road_segment = 1:numel(road_specifications)
    seg = road_specifications(ith_road_segment);

    switch lower(seg.type)
        case 'straight'
            L = seg.length;
            % kap = 0;
        case 'arc'
            if isfield(seg, 'length') && ~isempty(seg.length)
                L = seg.length;
            else
                L = abs(seg.radius) * abs(seg.angle);
            end
            % kap = sign(seg.angle) / abs(seg.radius);
        case 'spiral'
            L = seg.length;
        otherwise
            error('Unknown segment type: %s', seg.type);
    end

    % Build a subset for this segment
    nSteps = max(1, round(L / ds));
    s_local = s_current + linspace(0, L, nSteps + 1).';

    switch lower(seg.type)
        case 'straight'
            kappa_seg = zeros(nSteps, 1);
        case 'arc'
            k_arc = sign(seg.angle) / abs(seg.radius);
            kappa_seg = k_arc * ones(nSteps, 1);
        case 'spiral'
            % Curvature varies linearly with arc length from k_start to
            % k_end. This is the OpenDRIVE 'spiral' (clothoid) geometry,
            % required for any realistic transition between segments of
            % different curvature.
            
            s_param = linspace(0, L, nSteps + 1).';
            kappa_full = seg.k_start + (seg.k_end - seg.k_start) * s_param / L;
            kappa_seg = kappa_full(2:end);
    end

    % Grade: constant over the segment. Default to 0 if the field is
    % missing or empty.
    if isfield(seg, 'grade') && ~isempty(seg.grade)
        theta_seg = seg.grade * ones(nSteps, 1);
    else
        theta_seg = zeros(nSteps, 1);
    end

    % First point of this segment overlaps with last point of previous
    % segment (already in s_all), so drop the first point to avoid
    % duplication.
    s_all = [s_all; s_local(2:end)];          %#ok<AGROW>
    kappa_all = [kappa_all; kappa_seg];   %#ok<AGROW>
    theta_all = [theta_all; theta_seg];   %#ok<AGROW>

    segBounds(ith_road_segment,:) = [s_current, s_current + L];

    s_current = s_current + L;
end

s = s_all;
kappa = kappa_all;
theta = theta_all;

%% Integrate heading and global position
% 
psi = cumtrapz(s, kappa);
x = cumtrapz(s, cos(psi));
y = cumtrapz(s, sin(psi));

%% Apply friction spec 
% 
% Default mu read from element 1's mu_default if present, else 0.9. Each
% subsequent element overwrites mu in its [s_start, s_end] window. Later
% elements override earlier ones in any overlap.
if isfield(friction_specifications(1), 'mu_default') && ~isempty(friction_specifications(1).mu_default)
    mu_default = friction_specifications(1).mu_default;
else
    mu_default = 0.9;
end

mu = mu_default * ones(size(s));

for kth_friction_specification = 1:numel(friction_specifications)
    seg = friction_specifications(kth_friction_specification);
    if ~isfield(seg, 's_start') || isempty(seg.s_start)
        continue;
    end

    % % Determine the transition length for this patch
    % if isfield(seg, 'L_transition') && ~isempty(seg.L_transition)
    %     L_trans = seg.L_transition;
    % else
    %     L_trans = 0;
    % end

    % % Background mu on each side of this patch is whatever mu has been
    % % set so far at the patch boundaries. Read it from the current mu
    % % array rather than assuming mu_default, so that overlapping or
    % % adjacent patches behave sensibly.
    % [~, idx_start] = min(abs(s - seg.s_start));
    % [~, idx_end]   = min(abs(s - seg.s_end));
    % mu_bg_left  = mu(max(1, idx_start - 1));
    % mu_bg_right = mu(min(numel(s), idx_end + 1));

    
    in_seg = (s >= seg.s_start) & (s <= seg.s_end);
    mu(in_seg) = seg.mu;


    % % Apply linear ramps on either side if given. The ramps overwrite
    % % whatever was there in the windows [s_start - L_trans, s_start] and
    % % [s_end, s_end + L_trans].
    % if L_trans > 0
    %     % Left ramp: from mu_bg_left at (s_start - L_trans) down to seg.mu at s_start
    %     in_left = (s >= seg.s_start - L_trans) & (s < seg.s_start);
    %     if any(in_left)
    %         alpha = (s(in_left) - (seg.s_start - L_trans)) / L_trans;
    %         mu(in_left) = mu_bg_left + alpha .* (seg.mu - mu_bg_left);
    %     end
    % 
    %     % Right ramp: from seg.mu at s_end up to mu_bg_right at (s_end + L_trans)
    %     in_right = (s > seg.s_end) & (s <= seg.s_end + L_trans);
    %     if any(in_right)
    %         alpha = (s(in_right) - seg.s_end) / L_trans;
    %         mu(in_right) = seg.mu + alpha .* (mu_bg_right - seg.mu);
    %     end
    % end


end

%% Road output structure

road.s = s;
road.kappa = kappa;
road.theta = theta;
road.psi = psi;
road.x = x;
road.y = y;
road.segments = segBounds;

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
    % % check whether the figure already has data
    % temp_h = figure(figNum);
    % flag_rescale_axis = 0;
    % if isempty(get(temp_h,'Children'))
    %     flag_rescale_axis = 1;
    % end

    % Threshold below which a friction value is considered "low" for
    % figures. 
    LOW_MU_THRESHOLD = 0.5;

    % Detect low friction regions for the overlay. This works whether or
    % not a patch was specified, and works for multiple patches with
    % arbitrary mu values.
    in_low_mu = (mu <= LOW_MU_THRESHOLD);

    % Curvature ylim: scale to the largest kappa we actually have.
    max_kappa = max(abs(kappa));
    if max_kappa == 0
        kappa_ylim_top = 0.001;
    else
        kappa_ylim_top = 1.2 * max_kappa;
    end

    temp_h = figure(figNum);

    set(temp_h, 'Position', [100 100 900 800]);

    % Curvature vs arc length
    subplot(4,1,1);
    hold on;
    grid on;
    
    yesSpiral = 0; 
    for ith_road_segment = 1:numel(road_specifications)

        seg = road_specifications(ith_road_segment);

        s_start = segBounds(ith_road_segment,1);
        s_end = segBounds(ith_road_segment,2);
         
        % if ith_road_segment < numel(road_specifications)
        %     in_this_segment = (s >= s_start) & (s < s_end);
        % else
        %     in_this_segment = (s >= s_start) & (s <= s_end);
        % end
    
        in_this_segment = (s >= s_start) & (s <= s_end);

        switch lower(seg.type)
            case 'straight'
                this_color = 'g';   % green
            case 'spiral'
                this_color = 'b';   % blue
                yesSpiral = 1; 
            case 'arc'
                this_color = 'm';   % magenta
            otherwise
                this_color = 'k';
        end

        plot(s(in_this_segment), kappa(in_this_segment), ...
            'Color', this_color, ...
            'LineWidth', 2.5);

    end

    xlabel('Arc length s (m)');
    ylabel('\kappa(s) (1/m)');
    % title('Road curvature');
    ylim([-kappa_ylim_top, kappa_ylim_top]);

    h_straight = plot(nan, nan, 'g', 'LineWidth', 2.5);
    
    h_arc = plot(nan, nan, 'm', 'LineWidth', 2.5);
    
    if yesSpiral
        h_spiral = plot(nan, nan, 'b', 'LineWidth', 2.5);
        legend([h_straight, h_spiral, h_arc], {'Straight', 'Spiral', 'Arc'}, 'Location', 'best');

    else
        legend([h_straight, h_arc], {'Straight', 'Arc'}, 'Location', 'best');
    end

   
    % Friction supply vs arc length with low mu overlay
    subplot(4,1,2);
    hold on;
    grid on;

    % Shade low mu regions first (so the line draws on top).
    if any(in_low_mu)
        yl_patch = [0, 1.0];

        % Find low mu regions and shade each one.
        lowMU_regions = fcn_INTERNAL_findLowMURegions(in_low_mu);
        for rth_lowMuReg = 1:size(lowMU_regions, 1)
            sx = [s(lowMU_regions(rth_lowMuReg,1)), s(lowMU_regions(rth_lowMuReg,2)), s(lowMU_regions(rth_lowMuReg,2)), s(lowMU_regions(rth_lowMuReg,1))];
            sy = [yl_patch(1), yl_patch(1), yl_patch(2), yl_patch(2)];
            patch(sx, sy, [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        end
    end

    % plot(s, mu, 'LineWidth', 1.5, 'Color','g');
    
    yesSpiral = 0;
    for ith_road_segment = 1:numel(road_specifications)

        seg = road_specifications(ith_road_segment);

        s_start = segBounds(ith_road_segment,1);
        s_end = segBounds(ith_road_segment,2);

        if ith_road_segment < numel(road_specifications)
            in_this_segment = (s >= s_start) & (s < s_end);
        else
            in_this_segment = (s >= s_start) & (s <= s_end);
        end

        switch lower(seg.type)
            case 'straight'
                this_color = 'g';   % green
            case 'spiral'
                this_color = 'b';   % blue
                yesSpiral = 1;
            case 'arc'
                this_color = 'm';   % magenta
            otherwise
                this_color = 'k';
        end

        plot(s(in_this_segment), mu(in_this_segment), ...
            'Color', this_color, ...
            'LineWidth', 2.5);

    end

    h_straight = plot(nan, nan, 'g', 'LineWidth', 2.5);

    h_arc = plot(nan, nan, 'm', 'LineWidth', 2.5);

    if yesSpiral
        h_spiral = plot(nan, nan, 'b', 'LineWidth', 2.5);
        legend([h_straight, h_spiral, h_arc], {'Straight', 'Spiral', 'Arc'}, 'Location', 'best');

    else
        legend([h_straight, h_arc], {'Straight', 'Arc'}, 'Location', 'best');
    end

    xlabel('Arc length s (m)');
    ylabel('\mu(s)');
    % title('Friction supply (shaded = low-friction region)');
    ylim([0, 1.0]);

    % Road view in (x, y) (Topview)
    subplot(4,1,3);
    hold on;
    grid on;
    
    yesSpiral = 0;
    for ith_road_segment = 1:numel(road_specifications)

        seg = road_specifications(ith_road_segment);

        s_start = segBounds(ith_road_segment,1);
        s_end = segBounds(ith_road_segment,2);

        if ith_road_segment < numel(road_specifications)
            in_this_segment = (s >= s_start) & (s < s_end);
        else
            in_this_segment = (s >= s_start) & (s <= s_end);
        end

        switch lower(seg.type)
            case 'straight'
                this_color = 'g';   % green
            case 'spiral'
                this_color = 'b';   % blue
                yesSpiral = 1;
            case 'arc'
                this_color = 'm';   % magenta
            otherwise
                this_color = 'k';
        end

        plot(s(in_this_segment), theta(in_this_segment), ...
            'Color', this_color, ...
            'LineWidth', 2.5);

    end

    h_straight = plot(nan, nan, 'g', 'LineWidth', 2.5);

    h_arc = plot(nan, nan, 'm', 'LineWidth', 2.5);

    if yesSpiral
        h_spiral = plot(nan, nan, 'b', 'LineWidth', 2.5);
        legend([h_straight, h_spiral, h_arc], {'Straight', 'Spiral', 'Arc'}, 'Location', 'best');

    else
        legend([h_straight, h_arc], {'Straight', 'Arc'}, 'Location', 'best');
    end

    xlabel('Arc length s (m)');
    ylabel('Grade');

    % Road view in (x, y) (Topview)
    subplot(4,1,4);
    hold on; 
    grid on;

    plot(x, y, 'LineWidth', 2);
    plot(x(1), y(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot(x(end), y(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    % Overlay each low-mu region on the road geometry
    if any(in_low_mu)
        plot(x(in_low_mu), y(in_low_mu), 'r.', 'MarkerSize', 8);
        % legend_entries = {'Road centerline', 'Start', 'End', ...
        %                   sprintf('Low-\\mu region (\\mu \\leq %.2f)', LOW_MU_THRESHOLD)};

        legend_entries = {'Road centerline', 'Start', 'End'};
    else
        legend_entries = {'Road centerline', 'Start', 'End'};
    end

    axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    title('Topview of road geometry');
    legend(legend_entries, 'Location', 'best');

    % % Make axis slightly larger?
    % if flag_rescale_axis
    %     temp = axis;
    %     axis_range_x = temp(2)-temp(1);
    %     axis_range_y = temp(4)-temp(3);
    %     percent_larger = 0.3;
    %     axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
    % end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
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

function runs = fcn_INTERNAL_findLowMURegions(logicalVec)
% Returns an [Mx2] array of [startIdx endIdx] for each run of true values
% in the input logical column vector.

% Make sure the segement indices are in a column matrix
logicalVec = logicalVec(:);
d = diff([false; logicalVec; false]);

% Find the start and end indices of the segment
startIdx = find(d == 1);
endIdx = find(d == -1) - 1;
runs = [startIdx, endIdx];
end
