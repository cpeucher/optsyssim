function dfunc = num_diff_1d_pb(func,h)
% Differentiate evenly spaced 1D data with periodic boundary conditions
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the values of the derivative of a function whose
% values are specified as an input vector. The derivative is calculated as
% a finite difference
% f'(x)=(f(x+h)-f(x-h))/2h,
% where h is a constant step size.
% The calculation is performed via convolution of an extension of the input
% to 3 periods by a filter of the type [1 0 -1]/2/h.
% This calculation is only valid under the following assumptions:
% 1. Periodic boundary conditions are applied.
% 2. The samples are taken on the fixed step grid with spacing h.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% dfunc = num_diff_1d_pb(func,h);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% func              vector containing the values of the function to
%                       integrate [real vector]
%
%                       The values are specified on a constant-step grid.
%
% h                 step of the grid on which the input vector is
%                       specified [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% dfunc             values of the derivative of the function whose values
%                       are provided as input [real vector]
%
%                        The values of the derivative are specified on the
%                        same grid as the input data.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

func = func(:).';
% Force the input data to a line vector.

nsamples = length(func);
% Number of samples in the input data.

func = repmat(func,1,3);
% Extend the input array to three periods.

dfunc = conv(func,[1 0 -1]/(2*h),'full');
% Convolution between the extended array and the [1 0 -1] filter.

dfunc = dfunc(nsamples + 2:2*nsamples + 1);
% Extract the derivatives from the convolution, so that the derivative data
% is specified on the same grid as the input data.

end
% End of function