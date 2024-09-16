function sig = opt_phase_shift(sig,phase_shift)
% Optical phase shifter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a phase shift to a signal. The same phase shift is
% applied to both -x and -y polarisation components, i.e. the state of
% polarisation of the signal is not modified when calling this function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = opt_phase_shift(sig,phase_shift); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
% 
% phase_shift       phase shift to apply, in rad [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               phase shifted optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.x = sig.x.*exp(-1i*phase_shift);
sig.y = sig.y.*exp(-1i*phase_shift);

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------