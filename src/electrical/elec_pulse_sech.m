function sig = elec_pulse_sech(time,t0,tfwhm)
% Electrical hyperbolic secant pulse 
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a time domain sech pulse.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_sech(time_array,pulse_position,pulse_duration); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values at which points the pulse will be
%                       calculated, in s [real vector]
%
% tfwhm             pulse centre position, in s [real scalar]
%
% duration          pulse FWHM, in s [real scalar]
%
%                       This is the FWHM of the sech pulse itself, not of 
%                       its square. 
%                       Observe the difference with the opt_pulse_sech 
%                       function
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               electrical pulse [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

ts = tfwhm/(2*log(2 + sqrt(3)));
% Conversion from FWHM

sig = sech((time - t0)/ts);
% Pulse shape

end