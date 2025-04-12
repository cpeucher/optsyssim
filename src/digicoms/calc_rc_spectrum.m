function hh = calc_rc_spectrum(freq,ts,roll_off)
% Calculation of raised-cosine pulse spectrum
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the spectrum of a raised cosine pulse.
% The spectrum of the root raised-cosine is the square root of the spectrum
% returned by this function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% roll_off = 0.1;
% freq = linspace(0,fs/2,1000);
% H = calc_rc_spectrum(freq,symbol_rate,roll_off); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              frequency values at which the spectrum will be 
%                       calculated, in Hz [real vector]
%
% ts                sampling interval, in s [real scalar]
%
% roll_off          roll-off factor [real]
%
%                       0 <= roll_off <= 1
%                       For roll_off = 0, the pulse shape corresponds to a
%                       sinc pulse.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% hh                pulse spectrum [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

hh = ts*(1 + cos(pi*ts/roll_off*(abs(freq) - (1 - roll_off)/2/ts)))/2;
hh(abs(freq) < (1 - roll_off)/2/ts) = ts;
hh(abs(freq) > (1 + roll_off)/2/ts) = 0; 

end
