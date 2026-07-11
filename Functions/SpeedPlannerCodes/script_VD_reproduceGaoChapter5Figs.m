%% script_VD_reproduceGaoChapter5Figs
% Reproduction of Gao dissertation Figures 5-4, 5-5, and 5-6.
%
% Geometry parameters chosen to match Gao Fig 5-4:
%   - Turn radius R = 190 m 
%   - Spiral (clothoid) length L_s = 150 m each
%   - Bottom straight 300 m split by the start point at s = 150 m
%     (150 m before start, 150 m after end to close the loop)
%   - Top straight 300 m
%   - Total oval length approximately 2100 m (2093 m - tbp) (matches Fig 5-5 s-axis)
%
% Friction profile (Fig 5-5 middle plot): background mu = 0.85 with
% four friction patches at mu = 0.3, 0.4, 0.15, 0.3 in sequence,
% distributed across the second half of turn 1, the middle straight,
% and the entry to turn 2.
%
% Grade profile (Fig 5-5 bottom panel): +6% (approx +3.43 deg) from
% s = 100 to s = 350 (uphill entering turn 1); -6% (approx -3.43 deg)
% from s = 700 to s = 1050 (downhill exiting turn 1). Zero every where else.
%
% Speed planner parameters (Gao Section 5.1.4):
%   U_posted = 50 m/s, lambda = 0.95
%   Initial and terminal speeds at 50 m/s.
%

% REVISION HISTORY:
% 2026_06_28 by Aneesh Batchu, abb6486@psu.edu
% - In fcn_VD_buildSampleRoad
%   % * Wrote the code originally

%% Set up the workspace
close all

%% Oval geometry (Approximated based on Fig 5-5(a))

R_turningRadius = 190;   % turn radius of arc in m
kappa_turningRadius = 1 / R_turningRadius; % arc curvature approx 5e-3
L_spiral = 150; % spiral length in m
dHeadingAngle_spiral = L_spiral * (0 + kappa_turningRadius) / 2; % heading change per spiral
arcHeadingAngle = pi - 2 * dHeadingAngle_spiral;  % arc angle for a full pi turn

L_botStraight_pre = 150; % from start heading east
L_topStraight = 300; % top straight (west)
L_botStraight_post = 150; % closing the loop back to start

%% Road specification: 9 segments building a closed oval, starting in
% the middle of the bottom straight and ending where it began.
% Grade is set to 0 here. We override road.theta below after building.

clear roadSpec

roadSpec(1).type   = 'straight';   
roadSpec(1).length = L_botStraight_pre;
roadSpec(1).grade  = [];

roadSpec(2).type    = 'spiral';    
roadSpec(2).length  = L_spiral;
roadSpec(2).k_start = 0;           
roadSpec(2).k_end   = kappa_turningRadius;
roadSpec(2).grade   = [];

roadSpec(3).type    = 'arc';        
roadSpec(3).radius  = R_turningRadius;
roadSpec(3).angle   = arcHeadingAngle;
roadSpec(3).grade   = [];

roadSpec(4).type    = 'spiral';    
roadSpec(4).length  = L_spiral;
roadSpec(4).k_start = kappa_turningRadius;      
roadSpec(4).k_end   = 0;
roadSpec(4).grade   = [];

roadSpec(5).type   = 'straight';   
roadSpec(5).length = L_topStraight;
roadSpec(5).grade  = [];

roadSpec(6).type    = 'spiral';    
roadSpec(6).length  = L_spiral;
roadSpec(6).k_start = 0;           
roadSpec(6).k_end   = kappa_turningRadius;
roadSpec(6).grade   = [];

roadSpec(7).type   = 'arc';        
roadSpec(7).radius = R_turningRadius;
roadSpec(7).angle  = arcHeadingAngle;
roadSpec(7).grade  = [];

roadSpec(8).type    = 'spiral';    
roadSpec(8).length  = L_spiral;
roadSpec(8).k_start = kappa_turningRadius;      
roadSpec(8).k_end   = 0;
roadSpec(8).grade   = [];

roadSpec(9).type   = 'straight';   
roadSpec(9).length = L_botStraight_post;
roadSpec(9).grade  = [];


%% Friction specification: match Gao Fig 5-5 middle panel
% Background dry mu = 0.85, with four patches with different mu values.

clear frictionSpec
frictionSpec(1).mu_default = 0.85;
frictionSpec(1).s_start    = 420;
frictionSpec(1).s_end      = 790;
frictionSpec(1).mu         = 0.30;

frictionSpec(2).s_start    = 790;
frictionSpec(2).s_end      = 1180;
frictionSpec(2).mu         = 0.40;

frictionSpec(3).s_start    = 1180;
frictionSpec(3).s_end      = 1480;
frictionSpec(3).mu         = 0.15;

frictionSpec(4).s_start    = 1480;
frictionSpec(4).s_end      = 1790;
frictionSpec(4).mu         = 0.35;


%% Build the road (with zero grade), then overwrite theta (grade) to match Gao

ds = 1;
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, -1);

L_total = road.s(end);

fprintf(1,'Total oval length: %.1f m (Gao Fig 5-5 s-axis reaches ~2100)\n', L_total);

% Grade override: regions matching Gao Fig 5-5(c).
% +6% grade approximately = +3.43 deg = +0.0599 rad.
road.theta = zeros(size(road.s));
road.theta(road.s >= 100 & road.s <= 350) = +0.06;   % uphill
road.theta(road.s >= 850 & road.s <= 1050) = -0.06;   % downhill


%% Segment type coloring (for reproducing Figures 5-4 and 5-5 plots)

segIdx = struct('line', false(size(road.s)), ...
                'spiral', false(size(road.s)), ...
                'arc', false(size(road.s)));

for ith_seg = 1:numel(roadSpec)

    s_start = road.segments(ith_seg, 1);
    s_end = road.segments(ith_seg, 2);

    % Select the segement indices
    in_this = (road.s >= s_start) & (road.s <= s_end); 

    % Color based on the type (same as Gao)
    switch lower(roadSpec(ith_seg).type)
        case 'straight'
            segIdx.line(in_this) = true;
        case 'spiral'
            segIdx.spiral(in_this) = true;
        case 'arc'
            segIdx.arc(in_this) = true;
    end
end

% Modify the colors accordingly
col_line = [0.20 0.90 0.20];
col_spiral = [0.00 0.20 0.90];
col_arc = [0.95 0.15 0.85];

%% Figure 5-4 reproduction: oval in (x, y) with segment coloring

figNum = 21001;
titleString = sprintf('Figure 5-4 reproduction: oval in (x, y) with segment coloring');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

set(gcf, 'Position', [80 80 850 600]);

hold on
grid on 
axis equal;

% Shift to Gao's coordinate frame: start at (150, -200)
x_shift = 150;
y_shift = -200;

x_plot = road.x + x_shift;
y_plot = road.y + y_shift;

fcn_INTERNAL_localPlotBySeg(x_plot, y_plot, segIdx.line, col_line, '-',  3.5);
fcn_INTERNAL_localPlotBySeg(x_plot, y_plot, segIdx.spiral, col_spiral, '-.', 3.5);
fcn_INTERNAL_localPlotBySeg(x_plot, y_plot, segIdx.arc, col_arc, ':', 3.5);

plot(x_plot(1), y_plot(1), 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
text(x_plot(1) + 20, y_plot(1) + 20, 'start', 'FontSize', 12);

xlim([-300 600])

xlabel('xEast (m)'); 
ylabel('yNorth (m)');

% title('Reproduction of Gao Figure 5-4: circular oval sample-path');
h1 = plot(nan, nan, '-', 'Color', col_line, 'LineWidth', 3.5);
h2 = plot(nan, nan, '-.', 'Color', col_spiral, 'LineWidth', 3.5);
h3 = plot(nan, nan, ':', 'Color', col_arc, 'LineWidth', 3.5);
legend([h1 h2 h3], {'Line', 'Spiral', 'Arc'}, 'Location', 'best');

sgtitle(titleString, 'Interpreter','none');

%% Figure 5-5 reproduction: kappa(s), mu(s), theta(s) with segment coloring

figNum = 21002;
titleString = sprintf('Figure 5-5 reproduction: kappa(s), mu(s), theta(s) with segment coloring');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

set(gcf, 'Position', [80 80 900 800]);

subplot(3, 1, 1)
hold on; 
grid on;

fcn_INTERNAL_localPlotBySeg(road.s, road.kappa, segIdx.line, col_line, '-',  2.5);
fcn_INTERNAL_localPlotBySeg(road.s, road.kappa, segIdx.spiral, col_spiral, '-',  2.5);
fcn_INTERNAL_localPlotBySeg(road.s, road.kappa, segIdx.arc, col_arc, '-',  2.5);

xlabel('Station (m)'); 
ylabel('Curvature');

h1 = plot(nan, nan, '-', 'Color', col_line, 'LineWidth', 2.5);
h2 = plot(nan, nan, '-', 'Color', col_spiral, 'LineWidth', 2.5);
h3 = plot(nan, nan, '-', 'Color', col_arc, 'LineWidth', 2.5);

legend([h1 h2 h3], {'Line', 'Spiral', 'Arc'}, 'Location', 'north', 'Orientation', 'horizontal');
ylim([0, 6e-3]);
xlim([0, L_total]);

subplot(3, 1, 2)
hold on; 
grid on;

fcn_INTERNAL_localPlotBySeg(road.s, mu, segIdx.line, col_line, '-', 2.5);
fcn_INTERNAL_localPlotBySeg(road.s, mu, segIdx.spiral, col_spiral, '-', 2.5);
fcn_INTERNAL_localPlotBySeg(road.s, mu, segIdx.arc, col_arc, '-', 2.5);

xlabel('Station (m)'); 
ylabel('Friction Coefficient');
ylim([0, 1.0]);
xlim([0, L_total]);


subplot(3, 1, 3)
hold on; 
grid on;

grade_pct = tan(road.theta) * 100;

fcn_INTERNAL_localPlotBySeg(road.s, grade_pct, segIdx.line, col_line, '-', 2.5);
fcn_INTERNAL_localPlotBySeg(road.s, grade_pct, segIdx.spiral, col_spiral, '-', 2.5);
fcn_INTERNAL_localPlotBySeg(road.s, grade_pct, segIdx.arc, col_arc, '-', 2.5);

xlabel('Station (m)'); 
ylabel('Grade (%)');

ylim([-10, 10]);
xlim([0, L_total]);

sgtitle(titleString, 'Interpreter','none');


%% Speed planner with Gao's parameters
speedPlannerSpecs.U_posted   = 50.0;
speedPlannerSpecs.U_initial  = 50.0;
speedPlannerSpecs.U_terminal = 50.0;
speedPlannerSpecs.lambda     = 0.95;
speedPlannerSpecs.g          = 9.81;

[Ux_desired, speedData] = fcn_VD_computeSpeedProfile(road, mu, speedPlannerSpecs, -1);

%% Figure 5-6 reproduction: speed profile and accelerations (bottom)

figNum = 21003;
titleString = sprintf('Figure 5-6 reproduction: speed profile and accelerations');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

set(gcf, 'Position', [80 80 800 800]);

subplot(2, 1, 1)
hold on;
grid on;

plot(road.s, speedData.U_corner,   '-',  'LineWidth', 5.0, ...
    'Color', [0.85 0.70 0.10], 'DisplayName', 'curve limit speed');

plot(road.s, speedData.U_forward,  '--', 'LineWidth', 2.0, ...
    'Color', [0.20 0.85 0.20], 'DisplayName', 'forward limit');

plot(road.s, speedData.U_backward, '--', 'LineWidth', 2.0, ...
    'Color', [0.95 0.15 0.85], 'DisplayName', 'backward limit');

plot(road.s, Ux_desired, 'b-', 'LineWidth', 1.5, 'DisplayName', 'combined results');

xlabel('Station (m)'); ylabel('U_x (m/s)');
title('(a)');
legend('Location', 'south', 'Orientation', 'horizontal', 'NumColumns', 2);
xlim([0, L_total]);
ylim([0, 55]);


ax_profile = speedData.a_x;
ay_profile = speedData.a_y;
a_total = sqrt(ax_profile.^2 + ay_profile.^2); % Gao Eq. 5.15
a_max  = mu * speedPlannerSpecs.g .* cos(road.theta) + speedPlannerSpecs.g * abs(sin(road.theta)); % Gao Eq. 5.16

% % This will give some noise
% dUx_ds  = gradient(Ux_desired, road.s);
% ax_prof = Ux_desired .* dUx_ds;
% ay_prof = road.kappa .* (Ux_desired .^ 2);
% a_total = sqrt(ax_prof.^2 + ay_prof.^2);

% % with lambda
% a_max = speedPlannerSpecs.lambda * mu * speedPlannerSpecs.g .* cos(road.theta) + speedPlannerSpecs.g * abs(sin(road.theta));

subplot(2, 1, 2)

hold on; 
grid on;

plot(road.s, ax_profile, '-.', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.9], 'DisplayName', 'a_x');

plot(road.s, ay_profile, '--', 'LineWidth', 1.5, 'Color', [0.2 0.7 0.2], 'DisplayName', 'a_y');

plot(road.s, a_total, ':', 'LineWidth', 1.8, 'Color', [0.9 0.1 0.1], 'DisplayName', 'a_{total}');

plot(road.s, a_max, '-', 'LineWidth', 1.5, 'Color', 'k', 'DisplayName', 'a_{max}');

xlabel('Station (m)'); 
ylabel('Acceleration (m/s^2)');
title('(b)');
legend('Location', 'south', 'Orientation', 'horizontal', 'NumColumns', 4);
xlim([0, L_total]);
ylim([-10, 12]);

sgtitle(titleString, 'Interpreter','none');

%% Check the actual friction ellipse constraint
% Physical constraint: sqrt((a_x + g sin theta)^2 + a_y^2) <= lambda*mu*g*cos(theta)
% This is what the speed planner ODE enforces. a_total <= a_max in Fig
% 5-6(b) is a weaker one-directional comparison that can be violated on
% downhill without any physical infeasibility.

a_x_eff = ax_profile + speedPlannerSpecs.g * sin(road.theta);
a_total_eff = sqrt(a_x_eff.^2 + ay_profile.^2);
a_bound_ellipse = speedPlannerSpecs.lambda * mu * speedPlannerSpecs.g .* cos(road.theta);

figNum = 21004;
titleString = sprintf('Feasibility check: gravity augmented total demand vs friction ellipse bound');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

hold on
grid on

plot(road.s, a_total_eff, 'r-', 'LineWidth', 1.5, 'DisplayName', '\surd((a_x + g sin\theta)^2 + a_y^2)');

plot(road.s, a_bound_ellipse, 'k--', 'LineWidth', 1.5, 'DisplayName', '\lambda\mug cos\theta');

xlabel('Station (m)') 
ylabel('Acceleration (m/s^2)')

legend('Location', 'best');
xlim([0, L_total]);

assert(all(a_total_eff <= a_bound_ellipse + 5e-2));


sgtitle(titleString, 'Interpreter','none');

%% Assertions: Check the outputs

fprintf(1, 'Curve limit peak (should be U_posted=50): %.2f m/s\n', max(speedData.U_corner));

mask_dryArc = (road.kappa >= 0.9 * kappa_turningRadius) & (mu >= 0.8);
if any(mask_dryArc)
    fprintf(1, 'Curve limit on dry arc (mu=0.85, R=200): %.2f m/s\n', min(speedData.U_corner(mask_dryArc)));
end

mask_lowFricPatch = (mu < 0.2);
if any(mask_lowFricPatch)
    fprintf(1, 'Curve limit on low friction patch (mu=0.15): %.2f m/s\n', min(speedData.U_corner(mask_lowFricPatch)));
end

fprintf(1, 'Max a_total / a_max ratio: %.3f (must be <= 1)\n', max(a_total ./ max(a_max, eps)));

assert(all(a_total <= a_max + 1e-2));
assert(all(Ux_desired <= speedData.U_forward + 1e-6));
assert(all(Ux_desired <= speedData.U_backward + 1e-6));


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

function fcn_INTERNAL_localPlotBySeg(x, y, segIndices, color, style, line_width)

% Make sure the segement indices are in a column matrix
segIndices = segIndices(:);

d = diff([false; segIndices; false]);
% Find the start and end indices of the segment
start_idx = find(d == 1);
end_idx  = find(d == -1) - 1;

% Plot the segment and color it based on the type
for k = 1:numel(start_idx)
    seg_idx = start_idx(k):end_idx(k);
    plot(x(seg_idx), y(seg_idx), style, 'Color', color, 'LineWidth', line_width, 'HandleVisibility', 'off');
end

end % Ends fcn_INTERNAL_localPlotBySeg

