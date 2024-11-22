function y = dsp_upsample(x,L)
% Upsample by an integer factor by adding L - 1 zeros between samples
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function upsamples (expands) the input sequence by a factor L by
% adding L -  1 zeros between the input samples, where L is an integer.
% It returns a vector of the same type as the input vector (line or
% column).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input sequence [complex vector]
%
% L                 upsample factor [integer]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% y                 uspsampled sequence [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if mod(L,1) ~= 0
    error('dsp_upsample: the upsampling factor should be an integer.');
end
% The upsampling factor should be an integer.

input_size = size(x);
% Determine size of input sequence

if input_size(1) > 1 && input_size(2) > 1
    error('dsp_upsample: the input sequence should be provided as a vector.');
end

x = x(:).';
% Force the input sequence to a line vector

y = [x;zeros(L-1,length(x))];
% Add L -1 lines of zeros

y = reshape(y,[1,length(x)*L]);
% Interleave the samples from the lines

if input_size(1) > 1
    % The input sequence is a column vector.
    y = y.';
    % Then we ensure the output sequence is also a column vector.
end

end