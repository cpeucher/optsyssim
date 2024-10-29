function y = func_sinc(x)
% sinc function that can be applied to vectors containing zeros
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a sinc function that is compatible with vectors
% containing zeros. The sinc function is defined here according to:
% sinc(x) = sin(x)/x.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% y = func_sinc(x); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input values [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% x                 output values [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if ~isempty(find(x == Inf))    
    error('func_sinc: the input array should not contain Inf.');    
end
% First we check that we do not attempt to calculate sinc(Inf).

% In this case the only option to have NaN in the output array is when the
% input is equal to 0.

x(x == 0) = eps(0);
% Replace 0 by something very small in the input array

y = sin(x)./x;
% Calculate sin(x)/x. This will return 1 for x == eps(0).

end