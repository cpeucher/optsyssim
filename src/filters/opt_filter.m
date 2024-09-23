function sig = opt_filter(sig,tf)
% Optical filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a calculated optical transfer function to an 
% optical signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = opt_filter(sig,tf);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% tf                filter (complex) transfer function defined on the
%                       frequency_array frequency grid [complex vector]
%
%                       tf is specified in increasing frequency order.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output filtered optical signal 
%                       [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.x = ifft(fft(sig.x).* fftshift(tf));
sig.y = ifft(fft(sig.y).* fftshift(tf));
% Applies the filter to the input signal.

end