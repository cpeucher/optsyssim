function pulse = opt_pulse_sech(time,peak_power,position,duration,chirp)
% Hyperbolic secant pulse shape with linear chirp
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates an optical hyberbolic secant pulse.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% pulse = opt_pulse_sech(time_array,peak_power,position,duration,chirp);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time vector at which points the pulse will be
%                       evaluated, in s [real vector]
%
% peak_power        pulse peak power, in W [real scalar]
%
% position          pulse centre position, in s [real scalar]
%
% duration          pulse FWHM duration (defined on the intensity), in s
%                       [real scalar]
%
% chirp             linear chirp parameter of the pulse, no unit 
%                       [real scalar]
%                       In case the chirp is of the form e^(0.5*Cp*t^2),
%                       as is often the case in the literature, 
%                       the correspondance with our chirp parameter
%                       definition is:
%                       Cp = C*(4*log(2)/duration^2)
% 
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% pulse             complex envelope of the pulse [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

t0 = duration/(2*log(1 + sqrt(2)));
% Calculate t0 from the expected FWHM.
% see e.g. Agrawal, "Nonlinear Fiber Optics," 3rd edition, Academic Press, 
% pp. 71-72, 2001.

pulse = sqrt(peak_power)*sech((time - position)/t0).*exp(0.5*1i*chirp*((time - position)/t0).^2);

end