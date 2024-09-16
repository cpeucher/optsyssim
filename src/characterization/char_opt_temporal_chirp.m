function chirp = char_opt_temporal_chirp(sig)
% Calculate temporal chirp of an optical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the instantaneous frequency (with respect to the
% reference_frequency) of an optical signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% chirp = char_opt_temporal_chirp(sig);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% chiro              instantaneous frequency variations with respect to
%                       reference_frequency, in Hz [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% When the periodic boundary condition is not satisfied, one needs to
% truncate the chirp vector and the corresponding time vector
% chirp = chirp(2:end -1);
% time_chirp = time_array(2:end -1);
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt

phi = -unwrap(angle(sig.x));
% Signal phase.
% We employ the phase convention exp(-1j*phi(t)) while angle in matlab
% returns the argument of a complex number.
% This is as good as unwrap is...

chirp = -num_diff_1d_pb(phi,dt)/2/pi;
% Instantaneous frequency

end