function sig = mod_pm(sig,drive,vpi,loss)
% Phase modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models a phase modulator. No bandwidth limitation is
% introduced. Therefore, if needed, such limitation needs to be modelled by
% low pass filtering the driving signals.
% Only the -x polarisation of the incoming signal is modulated. 
% The -y polarisation is blocked (equivalent to having a polariser // x at
% the input of the modulator).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% vpi = 1.0;
% loss = 0;
% sig = mod_pm(sig,nrz_data_sig,vpi,loss); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% drive             electrical driving signal applied to the modulator
%                       [real vector]
%
% vpi               half-wave voltage of the modulator [real scalar]
%
% loss              modulator loss, in dB (positive number) [real scalar]              
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.x = 10^(-loss/20)*sig.x.*exp(-1i*pi*(drive/vpi));
sig.y = zeros(1,length(sig.y));

end