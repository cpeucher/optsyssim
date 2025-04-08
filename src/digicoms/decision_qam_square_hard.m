function symbs_cx = decision_qam_square_hard(symbs_rx,m)
% Hard decision for square QAM
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs hard-decision for square QAM constellations.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs_cx = decision_qam_square_hard(symbs_rx,m);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs_rx          received signal samples [complex vector]
%
%                       The received signal has already been normalised 
%                       to the standard constellation so that the
%                       decision thresholds take values within
%                       {...,-6,-4,-2,0,2,4,6,...}
%
% m                 number of elements in the constellation
%                       [integer scalar]
%
%                       Only square constellations are considered
%                       therefore m = 2^k with k even
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs_cx          symbols after hard decision [complex vector]  
%
%                       Their real and imaginary parts take values 
%                       within {...,-7,-5,-3,-1,1,3,5,7,...}
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

MM = sqrt(m)-1;

re = (real(symbs_rx) >= 0).*min(2*fix(real(symbs_rx)/2) + 1,MM) + (real(symbs_rx) < 0).*max(2*fix(real(symbs_rx)/2) - 1,-MM);
im = (imag(symbs_rx) >= 0).*min(2*fix(imag(symbs_rx)/2) + 1,MM) + (imag(symbs_rx) < 0).*max(2*fix(imag(symbs_rx)/2) - 1,-MM);

symbs_cx = re + 1i*im;


end