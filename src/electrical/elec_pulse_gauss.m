function sig = elec_pulse_gauss(time,t0,tfwhm,order)
% Electrical Gaussian pulse shape
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a time domain Gaussian pulse.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_gauss(time_array,pulse_position,pulse_duration,pulse_order); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values at which points the pulse will be
%                   calculated, in s [real vector]
%
% t0                pulse centre position, in s [real scalar]
%
% tfwhm             pulse FWHM, in s [real scalar]
%
% order             Gaussian order [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               electrical pulse shape [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig = exp(-log(2)*(2*(time - t0)/tfwhm).^(2*order));
% Gaussian pulse

end