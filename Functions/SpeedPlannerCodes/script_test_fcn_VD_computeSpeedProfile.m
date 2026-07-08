%% script_test_fcn_VD_computeSpeedProfile
% Tests: fcn_VD_computeSpeedProfile
%
% REVISION HISTORY:
% 2026_06_07 by Aneesh Batchu, abb6486@psu.edu
% - In script_test_fcn_VD_computeSpeedProfile
%   % * Wrote the code originally

%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
fprintf(1, 'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: Sample road with one low-friction patch in the curve
figNum = 10001;
titleString = sprintf('DEMO case: Speed profile on a sample road with one low friction patch');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

% Build the sample road: 200 m straight + 90-deg left turn at R = 200 m
% + 200 m straight. One low-friction patch inside the curve.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];


roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum + 1000));

% Planner settings: posted speed 32 m/s, friction margin lambda = 0.9,
% initial and terminal speeds at the posted speed (cruising entry/exit).

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

% Call the function with plotting on
[Ux_desired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (figNum));

sgtitle(titleString, 'Interpreter', 'none');

% Variable types
assert(isnumeric(Ux_desired));
assert(isstruct(speedProfileData));

% Variable sizes
assert(isequal(size(Ux_desired), size(road.s)));
assert(isequal(size(speedProfileData.U_corner), size(road.s)));
assert(isequal(size(speedProfileData.U_forward), size(road.s)));
assert(isequal(size(speedProfileData.U_backward), size(road.s)));

% Variable values: U_xd must not exceed posted speed anywhere
assert(all(Ux_desired <= plannerSettings.U_posted + 1e-9));

% On the low-friction patch the speed must be slower than the posted speed
in_patch = (road.s >= 350) & (road.s <= 450);
assert(max(Ux_desired(in_patch)) < plannerSettings.U_posted);

% Make sure plot opened up
assert(isequal(get(gcf, 'Number'), figNum));

%% DEMO case: straight road only - speed should hold at posted speed

figNum = 10002;
titleString = sprintf('DEMO case: straight road only - speed should hold at posted speed');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 500;
roadSpec(1).radius = [];
roadSpec(1).angle  = [];
roadSpec(1).grade  = [];

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = [];
frictionSpec(1).s_end      = [];
frictionSpec(1).mu         = [];

% Build a road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, 1, (figNum+1000));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

% Call the function
Ux_desired = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (figNum));

% On a straight road with cruising entry/exit, the speed should be
% exactly U_posted everywhere.
assert(all(abs(Ux_desired - plannerSettings.U_posted) < 1e-9));

%% DEMO case: Speed profile on a sample road with one very low friction (0.05) patch

figNum = 10003;
titleString = sprintf('DEMO case: Speed profile on a sample road with one very low friction (0.05) patch');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;


clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.05;

ds = 1;
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum+1000));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

% Call the function with plotting on
[Ux_desired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (figNum));

sgtitle(titleString, 'Interpreter', 'none');

% Variable types
assert(isnumeric(Ux_desired));
assert(isstruct(speedProfileData));

% Variable sizes
assert(isequal(size(Ux_desired), size(road.s)));
assert(isequal(size(speedProfileData.U_corner),   size(road.s)));
assert(isequal(size(speedProfileData.U_forward),  size(road.s)));
assert(isequal(size(speedProfileData.U_backward), size(road.s)));

% Variable values: U_xd must not exceed posted speed anywhere
assert(all(Ux_desired <= plannerSettings.U_posted + 1e-9));

% On the low-friction patch the speed must be slower than the posted speed
in_patch = (road.s >= 350) & (road.s <= 450);
assert(max(Ux_desired(in_patch)) < plannerSettings.U_posted);

% Make sure plot opened up
assert(isequal(get(gcf, 'Number'), figNum));


%% DEMO case: Speed profile on a sample road with two low friction patches and 0.6 default friction
 
figNum = 10004;
titleString = sprintf('DEMO case: Speed profile on a sample road with two low friction patches and 0.6 default friction');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.6;
frictionSpec(1).s_start    = 150;
frictionSpec(1).s_end      = 250;
frictionSpec(1).mu         = 0.4;

frictionSpec(2).mu_default = 0.6;
frictionSpec(2).s_start    = 450;
frictionSpec(2).s_end      = 550;
frictionSpec(2).mu         = 0.1;

ds = 1;
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum+1000));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

% Call the function 
[Ux_desired, speedProfileData] = ...
    fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (figNum));

sgtitle(titleString, 'Interpreter', 'none');

% Variable types
assert(isnumeric(Ux_desired));
assert(isstruct(speedProfileData));

% Variable sizes
assert(isequal(size(Ux_desired), size(road.s)));
assert(isequal(size(speedProfileData.U_corner),   size(road.s)));
assert(isequal(size(speedProfileData.U_forward),  size(road.s)));
assert(isequal(size(speedProfileData.U_backward), size(road.s)));

% Variable values: U_xd must not exceed posted speed anywhere
assert(all(Ux_desired <= plannerSettings.U_posted + 1e-9));

% On the low-friction patch 1 the speed must be slower than the posted speed
in_patch1 = (road.s >= 150) & (road.s <= 250);
assert(max(Ux_desired(in_patch1)) <= plannerSettings.U_posted);

% On the low-friction patch 2 the speed must be slower than the posted speed
in_patch2 = (road.s >= 450) & (road.s <= 550);
assert(max(Ux_desired(in_patch2)) <= plannerSettings.U_posted);

% Make sure plot opened up
assert(isequal(get(gcf, 'Number'), figNum));


%% Test cases start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
fprintf(1, 'Figure: 2XXXXXX: TEST mode cases\n');

%% TEST case: Effect of friction margin lambda

figNum = 20001;
titleString = sprintf('TEST case: Speed profile sweeping lambda (friction margin)');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, -1);

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.g = 9.81;

lambdas = [1.0, 0.9, 0.8, 0.7];
colors  = lines(numel(lambdas)); % Easiest way to assign different colors

hold on; 
grid on;

for ith_lambda = 1:numel(lambdas)
    plannerSettings.lambda = lambdas(ith_lambda);
    Uxd_i = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, -1);
    plot(road.s, Uxd_i, 'LineWidth', 1.5, ...
         'Color', colors(ith_lambda,:), ...
         'DisplayName', sprintf('\\lambda = %.1f', lambdas(ith_lambda)));
end
yline(plannerSettings.U_posted, ':', 'U_{posted}', 'Color', [0.4 0.4 0.4], 'DisplayName','U_{posted}');

xlabel('Arc length s (m)');
ylabel('U_{x,d}(s) (m/s)');
% title('Speed profile dependence on lambda');
legend('Location', 'best');

sgtitle(titleString, 'Interpreter', 'none');


%% TEST case: Effect of posted speed

figNum = 20002;
titleString = sprintf('TEST case: Speed profile sweeping U_posted');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);
figure(figNum); clf;

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, -1);

clear plannerSettings
plannerSettings.g = 9.81;
plannerSettings.lambda = 0.9;

posted_speeds = [25, 32, 40];
colors = lines(numel(posted_speeds));

hold on; 
grid on;
for jth_speed = 1:numel(posted_speeds)
    plannerSettings.U_posted = posted_speeds(jth_speed);
    plannerSettings.U_initial = posted_speeds(jth_speed);
    plannerSettings.U_terminal = posted_speeds(jth_speed);

    Uxd_j = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, []);
    plot(road.s, Uxd_j, 'LineWidth', 1.5, ...
        'Color', colors(jth_speed,:), ...
        'DisplayName', sprintf('U_{posted} = %.0f m/s', posted_speeds(jth_speed)));
end

temp = axis;
axis_range_x = temp(2)-temp(1);
axis_range_y = temp(4)-temp(3);
percent_larger = 0.1;
axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);

xlabel('Arc length s (m)');
ylabel('U_{x,d}(s) (m/s)');
% title('Speed profile dependence on posted speed limit');
legend('Location', 'best');

sgtitle(titleString, 'Interpreter', 'none');

%% TEST case: Convergence of RK4 under grid refinement

figNum = 20003;
titleString = sprintf('DEMO case: Convergence of RK4 under grid refinement');
fprintf(1, 'Figure %.0f: %s\n', figNum, titleString);

clear roadSpec
roadSpec(1).type   = 'arc';
roadSpec(1).radius = 100;           % tighter curve so cornering limit binds
roadSpec(1).length = [];
roadSpec(1).angle  = +pi/2;
roadSpec(1).grade  = [];

clear frictionSpec
frictionSpec(1).mu_default = 0.4;   % low friction so cornering limit binds
frictionSpec(1).s_start = [];
frictionSpec(1).s_end   = [];
frictionSpec(1).mu      = [];
% No patches; the friction is uniform at 0.7.

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 5.0;   % start slow so forward pass works hard
plannerSettings.U_terminal = 100.0; % effectively unbounded
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

ds_grid = [4, 2, 1, 0.5, 0.25];
end_values = zeros(size(ds_grid));

for ith_grid = 1:length(ds_grid)
    [roadArc, muArc] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds_grid(ith_grid), -1);
    Uxd_arc = fcn_VD_computeSpeedProfile(roadArc, muArc, plannerSettings, -1);
    end_values(ith_grid) = Uxd_arc(end);
end

fprintf(1, '\nConvergence values:\n');
for ith_grid = 1:length(ds_grid)
    fprintf(1, '  ds = %5.2f m  -->  U_x(s_end) = %.6f m/s\n', ds_grid(ith_grid), end_values(ith_grid));
end

figure(figNum); clf;
% plot(ds_grid, end_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogx(ds_grid, end_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Step size ds (m)');
ylabel('U_{x,d}(s_{end}) (m/s)');
title('Convergence of RK4 end-of-arc speed under grid refinement');
ylim([min(end_values) - 1e-5, max(end_values) + 1e-5]);
% sgtitle(titleString, 'Interpreter', 'none');

% Variable values: results at ds = 1, 0.5, 0.25 should agree to within
% 1e-3 m/s (RK4 has O(ds^4) error).
assert(abs(end_values(3) - end_values(5)) < 1e-3);

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
fprintf(1, 'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE

figNum = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

[UxdFast, ~] = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, ([]));
assert(isnumeric(UxdFast));
assert(isequal(size(UxdFast), size(road.s)));

figHandles = get(groot, 'Children');
assert(~any(figHandles == figNum));


%% Basic fast mode - NO FIGURE, FAST MODE

figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

[UxdFast2, ~] = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (-1));
assert(isnumeric(UxdFast2));
assert(isequal(size(UxdFast2), size(road.s)));

figHandles = get(groot, 'Children');
assert(~any(figHandles == figNum));


%% Compare speeds of normal mode vs fast mode
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum); close(figNum);

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = [];

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

clear plannerSettings
plannerSettings.U_posted   = 32.0;
plannerSettings.U_initial  = 32.0;
plannerSettings.U_terminal = 32.0;
plannerSettings.lambda     = 0.9;
plannerSettings.g          = 9.81;

Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations

    % Call the function
     [~, ~] = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;

for ith_test = 1:Niterations

    % Call the function
    [~, ~] = fcn_VD_computeSpeedProfile(road, mu, plannerSettings, (-1));

end
fast_method = toc;

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')

% Assertions
assert(isnumeric(UxdFast));
assert(isequal(size(UxdFast), size(road.s)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BUG

%% Fail conditions
if 1 == 0
    %% Fails because mu length doesn't match road.s

    figNum = 90001;
    fprintf(1,'Figure: %.0f:Bug case\n',figNum);
    figure(figNum); close(figNum);

    clear roadSpec
    roadSpec(1).type   = 'straight';
    roadSpec(1).length = 200;
    roadSpec(1).grade  = [];

    roadSpec(2).type   = 'arc';
    roadSpec(2).radius = 200;
    roadSpec(2).angle  = +pi/2;

    roadSpec(3).type   = 'straight';
    roadSpec(3).length = 200;

    clear frictionSpec
    frictionSpec(1).mu_default = 0.9;
    frictionSpec(1).s_start    = 350;
    frictionSpec(1).s_end      = 450;
    frictionSpec(1).mu         = 0.4;

    ds = 1;

    % Build the sample road
    [road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

    clear plannerSettings
    plannerSettings.U_posted = 32.0;
    plannerSettings.U_initial = 32.0;
    plannerSettings.U_terminal = 32.0;
    plannerSettings.lambda = 0.9;
    plannerSettings.g = 9.81;

    [~, ~] = fcn_VD_computeSpeedProfile(road, mu(1:end-5), plannerSettings, (figNum));

    % Make sure plot did NOT open up
    figHandles = get(groot, 'Children');
    assert(~any(figHandles==figNum));

    %% Fails because plannerSettings is missing required fields

    figNum = 90002;
    fprintf(1,'Figure: %.0f:Bug case\n',figNum);
    figure(figNum); close(figNum);

    clear roadSpec
    roadSpec(1).type   = 'straight';
    roadSpec(1).length = 200;
    roadSpec(1).grade  = [];

    roadSpec(2).type   = 'arc';
    roadSpec(2).radius = 200;
    roadSpec(2).angle  = +pi/2;

    roadSpec(3).type   = 'straight';
    roadSpec(3).length = 200;

    clear frictionSpec
    frictionSpec(1).mu_default = 0.9;
    frictionSpec(1).s_start    = 350;
    frictionSpec(1).s_end      = 450;
    frictionSpec(1).mu         = 0.4;

    ds = 1;

    % Build the sample road
    [road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

    badSettings = struct('U_posted', 32);
    [~, ~] = fcn_VD_computeSpeedProfile(road, mu, badSettings, []);

    % Make sure plot did NOT open up
    figHandles = get(groot, 'Children');
    assert(~any(figHandles==figNum));
end


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

% (no local helper functions in this script)