function sig = elec_pulse_rc(time,t0,ts,alpha)
% Raised-cosine pulse shape
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates a raised-cosine pulse with adjustable duration
% and roll-off factor.
% With the raised-cosine pulse definition used in this function:
% Pulse peak amplitude: 1
% Pulse energy: (4 - alpha)*ts/4
% Spectrum (Fourier transform of pulse) max amplitude: ts
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_rc(time_array,pulse_position,1/symbol_rate,alpha);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values over which the pulse will be calculated
%                       [real vector]
%
% t0                position of the pulse, in s [real scalar]
%
% ts                sampling interval, in s [real scalar]
%                       
%                       ISI-free pulse when ts = 1/Rs, where Rs is the
%                       symbol rate.
%
% alpha             roll-of factor [real scalar]
%
%                       The roll-off factor takes values in [0,1].
%                       When alpha = 0; the pulse corresponds to a 
%                       sinc pulse.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% pulse             generated pulse [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

a = pi*(time - t0)/ts;
% Argument of the function.

denom = (1 - (2*alpha*a/pi).^2);
sig = func_sinc(a).*cos(alpha*a)./denom;
% Raised-cosine pulse

find_zero = find(abs(denom) < eps);
sig(find_zero)= pi/4*func_sinc(a(find_zero));
% We need to take special care when the denominator is zero.

end