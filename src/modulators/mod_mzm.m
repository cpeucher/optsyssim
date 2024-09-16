function sig = mod_mzm(sig,drive1,drive2,bias1,bias2,vpi,split1,split2,loss)
% Mach-Zehnder modulator 
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models a Mach-Zehnder modulator. No bandwidth limitation is
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
% driving_signal_1 = vpi/2*(nrz_data_sig-0.5);
% driving_signal_2 = -vpi/2*(nrz_data_sig-0.5);
% bias_1 = 1.5*vpi;
% bias_2 = 0;
% split_in = 0.5;
% split_out = 0.5;
% loss = 0;
% sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,split_in,split_out,loss);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% drive1            electrical driving signal applied to the upper arm
%                       [real vector]
%
% drive2            electrical driving signal applied to the lower arm
%                       [real vector]
%
% bias1             dc bias applied to the the upper arm [real scalar]
%
% bias2             dc bias applied to the lower arm [real scalar]
%
% vpi               half wave voltage of the modulator [real scalar]
%
% split1            power splitting ratio of the input Y-junction
%                       [real scalar]
%
% split2            power splitting ratio of the output Y-junction
%                       [real scalar]
%
% loss              modulator insertion loss in dB (positive number)
%                       [real scalar]
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

a = sqrt(split1*split2);
b = sqrt((1 - split1)*(1 - split2));

sig.x = 10^(-loss/20)*sig.x.*(a*exp(-1i*pi*(drive1 + bias1)/vpi) + b*exp(-1i*pi*(drive2 + bias2)/vpi));
sig.y = zeros(1,length(sig.y));

end