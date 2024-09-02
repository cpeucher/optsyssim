function [sigout1,sigout2,sigout3,sigout4] = opt_hybrid90_2x4(sigin1,sigin2)
% 90 degree 2x4 optical hybrid
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function simulates a 90 degree 2x4 optical hybrid. The structure of
% the hybrid consists of 2 input 3 dB couplers, 2 output 3 dB couplers and
% a 90 degree phase shifter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [sigout1,sigout2,sigout3,sigout4] = opt_hybrid90_2x4(sig,sig_lo); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sigin1            input optical signal 1 [optical signal structure]
%
% sigin2            input optical signal 2 [optical signal structure]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sigout1           output optical signal 1 [optical signal structure]
%
% sigout2           output optical signal 2 [optical signal structure]
%
% sigout3           output optical signal 3 [optical signal structure]
%
% sigout4           output optical signal 4 [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
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

[sig11,sig12] = opt_coupler_2x2(sigin1,opt_no_sig,'lin',0.5);
% Input 3 dB coupler (top).
[sig21,sig22] = opt_coupler_2x2(opt_no_sig,sigin2,'lin',0.5);
% Input 3 dB coupler (bottom).
[sigout1,sigout2] = opt_coupler_2x2(sig11,sig21,'lin',0.5);
% Output 3 dB coupler(top).
[sigout3,sigout4] = opt_coupler_2x2(sig12,opt_phase_shift(sig22,pi/2),'lin',0.5);
% Output 3 dB coupler (bottom).

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------