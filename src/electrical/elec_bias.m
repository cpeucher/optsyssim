function sig = elec_bias(sig,dc,pp)
% Applies bias and set peak-to-peak value to electrical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function adjusts an electrical signal to desired peak-to-peak 
% and dc levels.
% Warning: this function uses the elec_dc_block function, which works for
% "well balanced" signal. The number of marks and spaces should be ~ the
% same.
% Does not work with noisy or strongly distorted signals.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% dc_level = 0.0;
% peak_peak = 1.0;
% sig = elec_bias(sig,dc_level,peak_peak); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% dc                dc level of the output signal [real scalar]
%
% pp                peak-to-peak level of the output signal [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output electrical signal [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[sig, ~] = elec_dc_block(sig);
% Remove the dc content from the input signal

ppin = max(sig) - min(sig);
% Extract its peak-to-peak value

sig = dc + pp/ppin*sig;
% Adjust the dc and peak-to-peak levels

end