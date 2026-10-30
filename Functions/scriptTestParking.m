p1 = [1 1; 4 2];
p2 = [2 5; 5 4];

% Create figure and attach drag callbacks
fig = figure('Name','Draggable Line Segments', ...
    'WindowButtonDownFcn',@startDrag, ...
    'WindowButtonMotionFcn',@drag, ...
    'WindowButtonUpFcn',@stopDrag);

% Set up axes for interactive display
ax = axes(fig);
hold(ax,'on');
axis(ax,'equal');
grid(ax,'on');

% Plot the two movable segments and their connector
h1 = plot(ax,p1(:,1),p1(:,2),'b-','LineWidth',6);
h2 = plot(ax,p2(:,1),p2(:,2),'r-','LineWidth',6);
hc = plot(ax,[p1(2,1) p2(2,1)],[p1(2,2) p2(2,2)], ...
    'k-','LineWidth',1);

% Create hidden distance label at the connector midpoint
ht = text(ax,mean([p1(2,1) p2(2,1)]),mean([p1(2,2) p2(2,2)]),'', ...
    'HorizontalAlignment','center','BackgroundColor','w', ...
    'Visible','off');

% Store handles and drag state in figure app data
d.p = {p1,p2};
d.lines = [h1 h2];
d.connector = hc;
d.label = ht;
d.active = [];
d.mode = [];
d.offset = [];
guidata(fig,d);

function startDrag(src,~)
d = guidata(src);
q = src.CurrentAxes.CurrentPoint(1,1:2);

% Find the nearest segment endpoint under the cursor
ends = [d.p{1}; d.p{2}];
[dist,idx] = min(vecnorm(ends-q,2,2));

% Ignore clicks that are too far from any endpoint
if dist > 0.5
    return
end

% Select the target segment and whether its start or end is moving
d.active = ceil(idx/2);
d.mode = mod(idx-1,2)+1;

% Preserve the grabbed-point offset when dragging whole segments
if d.mode == 1
    d.offset = q - d.p{d.active}(1,:);
else
    d.offset = [];
end

guidata(src,d);
end

function drag(src,~)
d = guidata(src);

% Exit early unless a drag is currently active
if isempty(d.active)
    return
end

q = src.CurrentAxes.CurrentPoint(1,1:2);
p = d.p{d.active};

% Update either the whole segment or just its free endpoint
if d.mode == 1
    delta = q - d.offset - p(1,:);
    p = p + delta;
else
    p(2,:) = q;
end

% Refresh segment, connector, and distance label
d.p{d.active} = p;
set(d.lines(d.active),'XData',p(:,1),'YData',p(:,2));

a = d.p{1}(2,:);
b = d.p{2}(2,:);
set(d.connector,'XData',[a(1) b(1)],'YData',[a(2) b(2)]);

m = (a+b)/2;
dist = norm(b-a);
set(d.label,'Position',[m 0],'String',sprintf('%.2f',dist),'Visible','on');

guidata(src,d);
drawnow limitrate
end

function stopDrag(src,~)
d = guidata(src);
d.active = [];
d.mode = [];
d.offset = [];
guidata(src,d);
end