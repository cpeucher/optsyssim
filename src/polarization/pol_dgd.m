function [sig,dgd_actual] = pol_dgd(sig,include_phase_shift,dgd_target)
% Differential group delay
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a delay between the two polarisations components
% of the envelope of an optical signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% dgd_target = 50e-12;
% include_phase_shift = 0;
% [sig,dgd_actual] = pol_dgd(sig,include_phase_shift,dgd_target);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                    
% dgd_target        target differential group delay, in s [real scalar]
%                       If dgd_target is positive, then the -x 
%                       polarisation is delayed relatively to the -y 
%                       polarisation.
%                       If dgd_target is negative, then the -x 
%                       polarisation is advanced relatively to the -y 
%                       polarisation.
%
% include_phase_shift   specifies if a phase shift corresponding to the 
%                           delay is applied to the carrier of the signal
%                           [0/1]
%
%                           include_phase_shift = 1
%                               a phase shift of the carrier of
%                               -pi*reference_frequency*dgd_actual 
%                               is also applied to -x polarisation of the 
%                               signal, while a phase shift of 
%                               +pi*reference_frequency*dgd_actual
%                               is applied to the -y polarisation 
%                               component.
%
%                           include_phase_shift = 0
%                               no phase shifts are applied to the 
%                               polarisation components of the signal
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
% 
% dgd_actual        actual differential group delay, in s [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                    time samples separation, in s [real scalar]
%
% frequency_array       relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Christophe Peucheret (christophe.peucheret@univ-rennes1.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt 
global reference_frequency


delay_samples = round(dgd_target/dt/2);
% Calculate the integer number of samples resulting in a delay which is
% the closest to the target delay divided by 2 (since we will apply 
% \pm DGD/2 to the two polarisation components).
dgd_actual = 2*delay_samples*dt;
% Convert this number of samples into real delay.
shift_size = [0 delay_samples];
% Create array indicating that the column of the line vector where the
% signal is stored should be shifted by the proper number of samples.

if include_phase_shift == 1
    sig.x = circshift(sig.x,shift_size)*exp(-1i*2*pi*reference_frequency*dgd_actual/2);
    sig.y = circshift(sig.y,-shift_size)*exp(+1i*2*pi*reference_frequency*dgd_actual/2);
    % Shift the signal array by the right amount and apply the corresponding
    % phase shift.
elseif include_phase_shift == 0
    sig.x = circshift(sig.x,shift_size);
    sig.y = circshift(sig.y,-shift_size);
    % Shift the signal array by the right amount.
else
    error('pol_dgd: differential group delay option not implemented.')
end


end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------