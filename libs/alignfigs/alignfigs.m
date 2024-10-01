function alignfigs(nx)
%alignfigs   Arrange all open figure windows nicely.
%   ALIGNFIGS() arranges all the open current figure windows to fit nicely on
%   the screen in a grid as large as possible.
%
%   ALIGNFIGS(NX) specifies that the grid should have NX columns.
%
%   Note that if not all windows have the same size, they will be adjusted to be
%   so (to the smallest open window).
%
%   Examples:
%       close all, for k = 1:8, figure, end, alignfigs()
%       close all, for k = 1:15, figure, end, alignfigs(5)

% Nick Hale - v1.0 Dec 2014.
%           - v2.0 Jan 2019.

padding = true;
toolbars = true;
menubars = true;

h = findobj('Type', 'figure');
nh = numel(h);

% Strip toolbars and menubars if requested:
if ( ~toolbars )
    set(h, 'toolbar', 'none')
end
if ( ~menubars )
    set(h, 'Menubar', 'none')
end

% Get the screen size:
ss = get(0, 'ScreenSize');
ss = ss(3:4);
sd = get(0, 'ScreenDepth');
ss(1) = ss(1) - sd;
ss(2) = ss(2) - 2*sd;

if ( nh == 0 )
    % Nothing to do!
    return
elseif ( nh == 1 )
    p = zeros(1,4);
    p(1,3) = ss(1);
    p(1,4) = ss(2)-sd;
    p(1:2) = [sd/2 3*sd/2];
    set(h, 'OuterPosition', p);
    figure(1)
    return
end

% Sort into the right order:
[~, idx] = sort(cell2mat(get(h, 'Number')));
h = h(idx);

% Get the sizes:
r = {};
for k = 1:nh
    r{k} = get(h(k), 'Resize');
    set(h(k), 'Resize', 'off')
    p(k,:) = get(h(k), 'OuterPosition');
end

% Determine the grid size:
if ( nargin == 0 )
    nx = ceil(sqrt(nh));
end
ny = ceil(numel(h)/nx);
if ( padding )
    p(:,3) = .95*ss(1)/nx;
    p(:,4) = .95*ss(2)/ny;
else
    p(:,3) = ss(1)/nx;
    p(:,4) = ss(2)/ny;
end
p = round(p);

% Horizontal padding:
if ( padding )
    padx = round(ss(1)/(nx+1) - nx*p(1,3)/(nx+1));
    pady = round(ss(2)/(ny+1) - ny*p(1,4)/(ny+1));
else
    padx = 0;
    pady = 0;
end

% Place the first figure in the top left:
px = sd/2 + padx;
py = ss(2) - p(1,4) + sd;
p(1,1:2) = [px, py];
set(h(1), 'OuterPosition', p(1,:));

% Initialise px and py:
if ( nx > 1 )
    px = p(1,1) + p(1,3) + padx;
    py = p(1,2);
else
    % Need to treat single column as a special case.
    px = sd/2 + padx;
    py = p(1,2) - p(1,4) - pady;
end

% Loop over the remaining windows:
for k = 2:numel(h)
    % Set to the new px and py values
    p(k,1:2) = [px, py];
    set(h(k), 'OuterPosition', p(k,:));
   
    if ( mod(k, nx) )
        % Move to the next column:
        px = px + p(k-1,3) + padx;
    else
        % Move to the next row. 
        px = sd/2 +  padx;
        % Update py:
        py = py - p(k-1,4) - pady;
    end        
end

% Bring all the windows to the front:
for k = 1:nh
    figure(h(k))
    set(h(k), 'Resize', r{k});
end

end
