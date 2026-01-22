function space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints)
% Create Cartesian 2D space grid
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function creates a 2d Cartesian space grid to model transverse
% field.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% xrange = [-15e-6, 15e-6];  
% yrange = [-15e-6, 15e-6]; 
% nxpoints = 2001;
% nypoints = 2001;
% space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
% [space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% xrange            range of x values, in m [2-element real vector]
%
% yrange            range of x values, in m [2-element real vector]
%
% nxpoints          number of points along the x dimension [integer scalar]
%
% nypoints          number of points along the y dimension [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% space_grid        grid elements [structure]
%
%                       space_grid.x
%                           x-values [real vector]
%
%                       space_grid.y
%                           y-values [real vector]
%
%                       space_grid.dx
%                           step along x-dimension [real scalar]
%
%                       space_grid.dy    
%                           step along y-dimension [real scalar]
%
%                       space_grid.X
%                           x-values of all points in the grid 
%                           [real matrix]
%
%                       space_grid.Y     
%                           y-values of all points in the grid
%                           [real matrix]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

space_grid.x = linspace(xrange(1),xrange(2),nxpoints); % x values
space_grid.y = linspace(yrange(1),yrange(2),nypoints); % y values
space_grid.dx = space_grid.x(2) - space_grid.x(1);     % dx step
space_grid.dy = space_grid.y(2) - space_grid.y(1);     % dy step

[space_grid.X,space_grid.Y] = meshgrid(space_grid.x,space_grid.y); 
% Cartesian grid

end