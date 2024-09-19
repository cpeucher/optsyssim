function sig = elec_pulse_gaussian(time,position,duration,order)
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
% sig = elec_pulse_gaussian(time_array,pulse_position,pulse_duration,pulse_order); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values at which points the pulse will be
%                   calculated, in s [real vector]
%
% position          pulse centre position, in s [real scalar]
%
% duration          pulse FWHM, in s [real scalar]
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

sig = exp(-log(2)*(2*(time - position)/duration).^(2*order));
% Gaussian pulse.

end