function integral = num_int1d_simpson(func,dx)
% Numerical integration using Simpson rule
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates one dimensional integrals using Simpson rule.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% integral = num_int1d_simpson(func,dx)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% func              values of the function to integrate evaluated on a grid
%                       with constant step dx [real vector]
%                       The vector should contain an even number of 
%                       samples.
%
% dx                constant step size of the grid over which the function 
%                       to integrate is defined [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% integral         value of the integral [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sumeven = 0;
sumodd = 0;
M = length(func);

if rem(M,2) == 1
    error('num_int1d_simpson: the input vector should have an even number of elements');
end
% Check that the input vector has an even number of values

N = M - 1;
for i = 1:(N-1)/2
    sumeven = sumeven + func(2*i);
end

for i = 1:(N-3)/2
    sumodd = sumodd + func(2*i+1);    
end

integral = (dx/3)*(func(1) + func(N) + 4*sumeven + 2*sumodd);

end