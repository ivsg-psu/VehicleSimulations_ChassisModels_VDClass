
%% script_test_fcn_VD_buildSampleRoad
% tests fcn_VD_buildSampleRoad.m
% 
% REVISION HISTORY:
% 2026_06_04 by Aneesh Batchu, abb6486@psu.edu
% - In script_test_fcn_VD_buildSampleRoad
%   % * Wrote the code originally

% TO-DO:
%
% 2026_06_04 by Aneesh Batchu, abb6486@psu.edu
% - (fill in items here)

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
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: A sample road with a positive curvature and a low friction patch
figNum = 10001;
titleString = sprintf('DEMO case: A sample road with a positive curvature and low friction patches');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Road specifications
% Three segments: 200 m straight, 90-degree left turn at R = 200 m,
% then 200 m straight.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;

% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.3;

% Step size
ds = 1; 

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));

%% DEMO case: A sample road with a positive curvature and low friction patche with a curvature

figNum = 10002;
titleString = sprintf('DEMO case: A sample road with a positive curvature and low friction patche with a curvature');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Road specifications
% Three segments: 200 m straight, 90-degree left turn at R = 200 m,
% then 200 m straight.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;


roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;


% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default   = 0.9;
frictionSpec(1).s_start      = 380;
frictionSpec(1).s_end        = 480;
frictionSpec(1).mu           = 0.1;
frictionSpec(1).L_transition = 15;   % 15 m linear ramp on each side

% Step size
ds = 1; 

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% DEMO case: A sample road with a positive curvature and low friction patches
figNum = 10003;
titleString = sprintf('DEMO case: A sample road with a positive curvature and low friction patches');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

% Road specifications
% Three segments: 200 m straight, 90-degree left turn at R = 200 m,
% then 200 m straight.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;


% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 100;
frictionSpec(1).s_end      = 200;
frictionSpec(1).mu         = 0.4;

frictionSpec(2).s_start    = 350;
frictionSpec(2).s_end      = 450;
frictionSpec(2).mu         = 0.3;

% Step size
ds = []; 

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


%% Test cases start here. These are very simple, usually trivial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _______ ______  _____ _______ _____
% |__   __|  ____|/ ____|__   __/ ____|
%    | |  | |__  | (___    | | | (___
%    | |  |  __|  \___ \   | |  \___ \
%    | |  | |____ ____) |  | |  ____) |
%    |_|  |______|_____/   |_| |_____/
%
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 2

close all;
fprintf(1,'Figure: 2XXXXXX: TEST mode cases\n');

%% TEST case: A sample road with a positive curvature connected by a spiral and a low friction patch
figNum = 20001;
titleString = sprintf('TEST case: A sample road with a positive curvature connected by a spiral and a low friction patch');
fprintf(1,'Figure %.0f: %s\n',figNum, titleString);
figure(figNum); clf;

R_curve  = 300;
k_arc    = 1/R_curve;
L_spiral = 60;
dHead_spiral = L_spiral * (0 + k_arc) / 2;   % heading change per spiral
arc_angle = pi/2 - 2 * dHead_spiral;       % so total turn = pi/2

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type    = 'spiral';
roadSpec(2).length  = L_spiral;
roadSpec(2).k_start = 0;
roadSpec(2).k_end   = k_arc;
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'arc';
roadSpec(3).radius = R_curve;
roadSpec(3).angle  = arc_angle;
roadSpec(3).grade  = 0;

roadSpec(4).type    = 'spiral';
roadSpec(4).length  = L_spiral;
roadSpec(4).k_start = k_arc;
roadSpec(4).k_end   = 0;
roadSpec(4).grade  = 0;

roadSpec(5).type   = 'straight';
roadSpec(5).length = 200;
roadSpec(5).grade  = 0;

% Patch sits inside the arc segment (between s ~= 260 and s ~= 514).
clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 380;
frictionSpec(1).s_end      = 480;
frictionSpec(1).mu         = 0.4;

ds = 1;

% Build a sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (figNum));


sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + L_spiral + ceil(roadSpec(3).radius * roadSpec(3).angle) + L_spiral + roadSpec(5).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),figNum));


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
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
figNum = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty figNum\n',figNum);
figure(figNum); close(figNum);

% Road specifications
% Three segments: 200 m straight, 90-degree left turn at R = 200 m,
% then 200 m straight.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;


% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.3;

% Step size
ds = 1; 

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, ([]));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% % Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Basic fast mode - NO FIGURE, FAST MODE
figNum = 80002;
fprintf(1,'Figure: %.0f: FAST mode, figNum=-1\n',figNum);
figure(figNum); close(figNum);

% Road specifications
% Three segments: 200 m straight, 90-degree left turn at R = 200 m,
% then 200 m straight.

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;

% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.3;

% Step size
ds = 1; 

% Build the sample road
[road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% % Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==figNum));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
figNum = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',figNum);
figure(figNum);
close(figNum);

clear roadSpec
roadSpec(1).type   = 'straight';
roadSpec(1).length = 200;
roadSpec(1).grade  = 0;

roadSpec(2).type   = 'arc';
roadSpec(2).radius = 200;
roadSpec(2).angle  = +pi/2;     % positive = left turn
roadSpec(2).grade  = 0;

roadSpec(3).type   = 'straight';
roadSpec(3).length = 200;
roadSpec(3).grade  = 0;

% Friction specifications

clear frictionSpec
frictionSpec(1).mu_default = 0.9;
frictionSpec(1).s_start    = 350;
frictionSpec(1).s_end      = 450;
frictionSpec(1).mu         = 0.3;

% Step size
ds = 1; 

Niterations = 50;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    [road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    [road, mu] = fcn_VD_buildSampleRoad(roadSpec, frictionSpec, ds, (-1));
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

% Check variable types
assert(isstruct(road));
assert(ismatrix(mu));

% Check variable sizes
assert(isequal(size(road.s,1), roadSpec(1).length + ceil(roadSpec(2).radius * roadSpec(2).angle) + roadSpec(3).length)); 
assert(isequal(size(mu,1), size(road.s,1))); 
assert(size(mu,2)==1); 

% % Check variable values
% assert(isequal(orientation,1));
% assert(isconvex);
% assert(isequal(Nencirclements,1));


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
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG 

%% Fail conditions
if 1==0
    %
        %% Fails because start_definition is not correct type
        


        
        %% Warning because start_definition is 3D not 2D
        % Start_zone definition is a 3D point [radius num_points X Y Z]
        


        
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
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


