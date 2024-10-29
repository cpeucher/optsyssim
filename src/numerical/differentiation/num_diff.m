function [xn,dy] = num_diff(order,x,y)
% Approximation of first or second-order derivatives by central finite differences
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Constant-step x vector is assumed. 
% Quick 'n dirty calculation of the derivatives at the midpoints xn
% The derivatives are estimated according to:
% dy/dx(xn)     = (y(xn + h) - y(xn - h))/(2h)
% d^2y/dx^2(xn) = (y(xn + h) - 2y(xn) + y(xn - h))/h^2
% where h is the step size.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [xn,dy] = num_diff(order,x,y);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% order             order of the derivative [integer]
%                       order = 1 or order = 2
%
% x                 x values [real vector]
%                       The x values should be evenly spaced.
%   
% y                 y = f(x) [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% xn                values at which the derivatives will be estimated 
%                       [real vector]
%
%                       The xn are the midpoints of the input x vector.
%
% dy                estimate of the derivative y' = df/dx or y''=d^2f/dx^2 
%                       at the points xn [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

dy = diff(y);
dx = x(2) - x(1);
% x step.

dys = circshift(dy,[0 -1]);

switch order
    % Switch over derivative order
    
    case 1
        % First-order central derivative
        
        dy = (dy + dys)/dx/2;
        
    case 2
        % Second-order central derivative     
        
        dy = (dys - dy)/dx^2;
        
    otherwise
        
        error('num_diff: function only defined for orders 1 and 2.');
end
% End of switch over derivative order.

dy = dy(1:end-1);
xn = x(2:end-1);

end