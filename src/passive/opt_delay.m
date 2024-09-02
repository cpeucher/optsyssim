function [sig,delay_actual] = opt_delay(sig,include_phase_shift,delay_target)
% Optical delay
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function delays the complex envelope of a signal by the closest
% allowed value (in terms of signal samples) from a certain target delay.
% The corresponding phase shift may or may not be applied. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% include_phase_shift = 0;%1;
% delay_target = 100e-12;
% [sig,delay_actual] = opt_delay(sig,include_phase_shift,delay_target); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig                   input optical signal [optical signal structure]
%
% include_phase_shift   specifies if a phase shift is applied to the 
%                           signal [0/1]
%                               
%                           include_phase_shift = 1
%                               a phase shift of
%                               -2*pi*reference_frequency*delay is also 
%                               applied to the signal. 
%                               This behaviour corresponds to that of a 
%                               physical optical delay achieved by 
%                               propagation over a certain length of a 
%                               medium.
%
%                           include_phase_shift = 0
%                               no phase shift is applied to the signal
%
% delay_target          target delay, in s [real scalar]
%                           If delay_target is positive, the signal is
%                           delayed.
%                           If delay_target is negative, the signal is 
%                           advanced.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% delay_actual      actual delay that has been applied, corresponding to an 
%                       integer number of samples [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                    time samples separation, in s [real scalar]
%
% reference_frequency   reference frequency, in Hz [real scalar]
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


delay_samples = round(delay_target/dt);
% Calculate the integer number of samples resulting in a delay which is
% the closest to the target delay.
delay_actual = delay_samples*dt;
% Convert this number of samples into real delay.
shift_size = [0 delay_samples];
% Create array indicating that the column of the line vector where the
% signal is stored should be shifted by the proper number of samples.

if include_phase_shift == 1
    sig.x = circshift(sig.x,shift_size)*exp(-1i*2*pi*reference_frequency*delay_actual);
    sig.y = circshift(sig.y,shift_size)*exp(-1i*2*pi*reference_frequency*delay_actual);
    % Shift the signal array by the right amount and apply the 
    % corresponding phase shift.
elseif include_phase_shift == 0
    sig.x = circshift(sig.x,shift_size);
    sig.y = circshift(sig.y,shift_size);
    % Shifts the signal array by the right amount
else
    disp('opt_delay: optical delay option not implemented.')
end


end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------