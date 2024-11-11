function sig = elec_pulse_rrc(time,t0,ts,alpha)
% Root-raised-cosine pulse shape
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates a root raised-cosine pulse with adjustable 
% duration and roll-off factor.
% With the root-raised-cosine pulse definition used in this function:
% Pulse peak amplitude: (1 - alpha + 4*alpha/pi)/sqrt(ts)
% Pulse energy: 1
% Spectrum (Fourier transform of pulse) max amplitude: sqrt(ts)
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = elec_pulse_rrc(time_array,time_array(nsamples/2),1/symbol_rate,alpha); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% time              time values over which the pulse will be calculated
%                       [real vector]
%
% t0                position of the pulse, in s [real scalar]
%
% ts                symbol duration, in s [real scalar]
%
%                       This corresponds to 1/symbol_rate when the pulse is
%                       used to encode data at the symbol_rate.
%
% alpha             roll-of factor [real scalar]
%
%                       The roll-off factor takes values in [0,1].
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               generated pulse [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

a = pi*(time - t0)/ts;
% Argument of the function. Normalised time.

numerator = sin((1 - alpha)*a) + 4*alpha*a/pi.*cos((1 + alpha)*a);
% Numerator

denominator = 1 - (4*alpha*a/pi).^2;
% Element of the denominator

sig = numerator./denominator./a;
% Root raised-cosine pulse (not normalised)

% This may be calculated at some singular points: 
% (t - t0)/Ts = +1/alpha/4 or, equivalently, a = +pi/alpha/4;
% (t - t0)/Ts = -1/alpha/4 or, equivalently, a = -pi/alpha/4
% t = t0

% We manage these singularities:
tolerance = 1e-7;
sig(abs(a - pi/alpha/4) < tolerance) = alpha/sqrt(2)*(( 1 + 2/pi)*sin(pi/4/alpha)+(1 - 2/pi)*cos(pi/4/alpha));
sig(abs(a + pi/alpha/4) < tolerance) = alpha/sqrt(2)*(( 1 + 2/pi)*sin(pi/4/alpha)+(1 - 2/pi)*cos(pi/4/alpha));
sig(abs(a) <  tolerance) = (1 - alpha) + 4*alpha/pi;

% The pulse defined so far has:
% - an energy equal to ts
% - a peak value equal to (1 - alpha + 4*alpha/pi) 
% - a spectrum with maximum magnitude ts

% However, the raised-cosine pulse as defined in elec_pulse_rc and
% calc_rc_spectrum has a spectrum with magnitude ts. If we want to have a
% compatible root-raised-cosine pulse, the  maximum magnitude of the
% spectrum should be sqrt(ts).
% We therefore proceed o the normalisation:

sig = sig/sqrt(ts);
% Now the RRC pulse generated by elec_pulse_rrc is the inverse Fourier
% transform of the square root of the Fourier transform of the pulse 
% generated by elec_pulse_rc.
% In practice this does not really matter, but at least we have consistent
% defintions of the raised-cosine (elec_pulse_rc) and root-raised-cosine
% (elec_pulse_rrc) pulses.
% The pulse returned by this function has 
% - an energy equal to 1
% - a peak value equal to (1 - alpha + 4*alpha/pi)/sqrt(ts) 
% - a spectrum with maximum magnitude sqrt(ts)

end