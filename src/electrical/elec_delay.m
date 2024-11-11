function [sig,actual_delay] = elec_delay(sig,target_delay)
% Electrical delay line
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function delays an electrical signal by the closest allowed value
% (in terms of signal samples) from a certain target delay.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% target_delay = 1/symbol_rate/2;
% [sig,actual_delay] = elec_delay(sig,target_delay);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               original electrical signal [real vector]
%
% target_delay      target delay, in s [real scalar]
%                       If target_delay is positive the signal is delayed.
%                       If target_delay is negative the signal is advanced.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               delayed electrical signal [real vector]
%
% actual_delay      actual delay that has been applied, corresponding to an
%                       integer number of samples [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt


sig = sig(:).';
% Ensure the signal is a line vector.

delay_samples = round(target_delay/dt);
% Calculates the integer number of samples resulting in a delay which is
% the closest to the target delay.

actual_delay = delay_samples*dt;
% Convert this number of samples into real delay.

sig = circshift(sig,[0 delay_samples]);
% Shifts the signal array by the right amounts.

end