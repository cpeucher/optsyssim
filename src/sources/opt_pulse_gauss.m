function pulse = opt_pulse_gauss(time,order,peak_power,position,duration,chirp)
% Gaussian pulse shape with linear chirp
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a linearly-chirped optical Gaussian pulse.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% pulse = opt_pulse_gauss(time_array,order,peak_power,position,duration,chirp);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time vector at which points the pulse will be
%                       evaluated, in s [real vector]
%
% order             Gaussian order, no unit [integer scalar]
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
%                       the correspondance with our chirp parameter C
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

pulse = sqrt(peak_power)*exp(-0.5*log(2)*(2*(time - position)/duration).^(2*order)).*exp(1i*0.5*chirp*log(2)*(2*(time - position)/duration).^2); 

end