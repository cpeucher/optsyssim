function sig = elec_pulse_rectangular(time,position,duration)
% Rectangular pulse shape
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a single rectangular pulse of specified duration.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_rectangular(time_array,pulse_position,pulse_duration);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values at which points the pulse will be 
%                       calculated, in s [real vector]
%
% position          pulse centre position, in s [real scalar]
%
% duration          pulse duration, in s [real scalar]
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

sig = (time >= position - duration/2) & (time < position + duration/2);
% Rectangular pulse

end
