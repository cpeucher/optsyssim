function I = num_int2d_simpson(F)
% Double integral of a function defined on a Cartesian coordinate system by Simpson method
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function integrates the function F defined by its values on the 
% rectangular grid with equal spacings space_grid using Simpson's method
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% F                 values of the function to integrate, defined as a 
%                       matrix on the Cartesian grid space_grid 
%                       [real matrix]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% I                value of the integral [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% space_grid       2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% Exploits functions from:
% Introduction to Numerical Methods and Matlab Programming for Engineers
% Todd Young & Martin J. Mohlenkamp
% Copyright © 2008, 2009, 2014, 2017, 2020, 2021 Todd Young and Martin J. Mohlenkamp.
% http://www.ohiouniversityfaculty.com/youngt/IntNumMeth [retrieved
% 2021-11-24]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global space_grid

[m,n] = size(space_grid.X);
% Extract space grid size.

if ~rem(m,2)
    error('num_int2d_Simpson: the number of points on the grid should be odd.')
end

m = m - 1;        % m is the number of intervals along the x axis.
n = n - 1;        % n is the number of intervals along the y axis.

xmin = space_grid.x(1);
xmax = space_grid.x(end);
% Extract min and max values of the grid along the x axis.

ymin = space_grid.y(1);
ymax = space_grid.y(end);
% Extract min and max values of the grid along the y axis.

W = mydblsimpweights(m,n);
% Weights of the Riemann sum.

I=(ymax-ymin)*(xmax-xmin)*sum(sum(W.*F))/(9*m*n);
% Integral according to Simpson rule.

end