function [v1,v2] = iq_gsop(u1,u2)
% Gram-Schmidt orthogonalization procedure for IQ imbalance compensation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the Gram-Schmidt orthogonalization procedure 
%  (GSOP) for IQ imbalance compensation.
% We follow:
% I. Fatadin, S. J. Savory, and D. Ives, "Compensation of quadrature 
% imbalance in an optical QPSK coherent receiver," IEEE Photon. Technol. 
% Lett. 20, 1733 (2008) doi: 10.1109/LPT.2008.2004630
% 
% The output vectors have unit power, i.e. mean(v1.^2) = mean(v2.^2) = 1.
% v1 is proportional to u1, but v1 and v2 are orthogonal, i.e.
% mean(v1.*v2) = 0;
%
% Statistics of the signal are calculated over the total length of the
% input vectors. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [samps_cmp_r,samps_cmp_i] = iq_gsop(real(samps),imag(samps));
% samps_cmp = samps_cmp_r + 1i*samps_cmp_i; 
% samps_cmp = normalise_constellation(samps_cmp,norm_es);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% u1                first input vector [real vector]
%
%                       In the IQ imbalance compensation usage, this
%                       corresponds to real samples at the output of the 
%                       in-phase balanced detector.
% 
% u2                second input vector [real vector]
%
%                       In the IQ imbalance compensation usage, this
%                       corresponds to real samples at the output of the 
%                       quadrature balanced detector.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% v1                first normalised vector [real vector]
%               
%                       v1 is proportional to u1, but it has been 
%                       normalised to a unit vector, i.e.
%
%                       mean(v1.^2) = 1
%
% v2                second orthonormalised vector [real vector]
%
%                       It is orthogonal to v1 and is also an unit vector, 
%                       i.e.
%
%                       mean(v2.^2) = 1
%                       mean(v1.*v2) = 0
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

v1 = u1/sqrt(mean(u1.^2));
a = u2 - u1*mean(u1.*u2)/mean(u1.^2);
v2 = a/sqrt(mean(a.^2));

end