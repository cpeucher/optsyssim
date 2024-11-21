function [y,z] = dsp_fir_linear(x,b,z)
% Finite impulse response filter (linear buffer implementation)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% Linear buffer implementation of a finite impulse response (FIR) filter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [y,z] = dsp_fir_linear(x,b,z);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input samples [real vector]
%
% b                 tap coefficients [real vector of length L]
%
% z                 initial state of the buffer [real vector of length L]
%                       The first element z[1] corresponds to the input.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% y                 output samples [real vector]
%
%                       y[n] = sum_{k=1}^L b[k] x[n - k + 1]
%
% z                 final state of the buffer [real vector of length L]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if ~isvector(x)
    error('dsp_fir_linear: input signal should be a vector');
end

if ~isvector(b)
    error('dsp_fir_linear: tap coefficients should be expressed as a vector');
end

if ~isvector(z)
    error('dsp_fir_linear: initial buffer state should be provided as a vector')
end

if length(b) ~= length(z)
    error('dsp_fir_linear: the number of tap coefficients should be equal to the size of the buffer')
end


z = z(:);
% Ensure z is a column vector

b = b(:).';
% Ensure b is a line vector


y = zeros(1,length(x));
% Pre-initialise output vector

for ii = 1:length(x)
    z = [x(ii);z(1:end-1)];
    % New state of the buffer
    y(ii) = b*z;
    % Filter output
end

end