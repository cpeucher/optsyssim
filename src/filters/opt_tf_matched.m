function tf = opt_tf_matched(pulse,freq)
% Transfer function of optical matched filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer function of an optical matched
% filter based on the pulse shape.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tf = opt_tf_matched(freq,pulse);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% pulse             pulse shape [real vector]
%   
%                       The pulse shape is defined on the standard
%                       time_array axis.
%
% freq              frequency values at which the transfer function is 
%                       calculated, in Hz
% 
%                       These are relative frequencies, i.e. the transfer
%                       function is centered at 0, corresponding to
%                       reference_frequency
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                values of the transfer function that applies to the 
%                       field [complex vector]
%
%                       Square for power transfer function.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array

tf = conj(fftshift(fft(pulse)));
% Calculate values of TF for the points in the frequency_array grid.
% In increasing frequency order.
    
tf = interp1(frequency_array,tf,freq);
% Interpolate at frequencies in freq.

end