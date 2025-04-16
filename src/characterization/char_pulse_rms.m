function rms = char_pulse_rms(power)
% rms width of an optical pulse
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the root-mean-square (rms) width of an (optical)
% pulse.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% rms = char_pulse_rms(abs(sig.x).^2); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% power             signal to characterise [real vector]
%
%                       If optical, it should be the power, abs(sig.x).^2
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% rms               rms width, in s [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt 
global time_array


D = num_int1d_simpson(power,dt);
% Integration of the pulse power <P (t)>, i.e. pulse energy

A = num_int1d_simpson(time_array.*time_array.*power,dt)/D;
% Normalised 2nd order moment
% <t^2 P(t)> / <P(t)>

B = num_int1d_simpson(time_array.*power,dt)/D;
% Normalised 1st order moment
% <t P(t)> / <P(t)>


rms = sqrt(A-B^2);
% rms pulse width

end