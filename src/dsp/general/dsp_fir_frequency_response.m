function [wh,H] = dsp_fir_frequency_response(h,varargin)
% Calculate the frequency response of a FIR filter from its impulse response / tap coefficients
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the frequency response of a FIR filter from its
% impulse response / tap coefficients
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [wh,H] = dsp_fir_frequency_response(h);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% h                 impulse response / tap coefficients of the FIR filter
%                       [vector]
%
%                       h[n], n = 0, ..., N-1, where N is the filter
%                       length.
%
% varargin          optional input argument [integer scaler]
%
%                       Number of points over which the frequency response
%                       will be computed.
%                       If not provided, a default value of 1000 is used.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% wh               normalized angular frequency range over which the
%                       frequency response has been calculated, in
%                       radians-per-sample.
%
% H                frequency response over the frequency range 
%                       normalized frequency: [0, 1/2] cycles-per-sample
%                       normalized angular frequency: [0, pi]
%                           radians-per-sample
%                       Thus corresponding to the real frequency [0,fs/2]
%                       (first Nyquist zone).
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

switch nargin
    case 1
        npoints = 1000;
    case 2
        npoints = varargin{1};
    otherwise
        error('dsp_fir_frequency_response: wrong number of input arguments.') 
end

h = h(:).';
% Ensures the impulse response is input as a line vector

wh = 2*pi*linspace(0,0.5,npoints);
% Normalized angular frequency over [0,pi] interval
n = [0:1:length(h) - 1].';
H = h*exp(-1i*n*wh);
% Matrix product:
% Impulse response h (line vector of N elements, where N is the filter length).
% Terms exp(-1i*wh*n) (matrix with N lines and npoints columns).
% We therefore have
% HH[k] = sum_{n= 0}^{N - 1} h[n] exp(-1i*wh[k]*n)



end