function [v1,v2] = iq_lop(u1,u2)
% Löwdin orthogonalization procedure for IQ imbalance compensation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the Löwdin orthogonalization procedure 
%  (LOP) for IQ imbalance compensation.
% We follow:
% Md.S. Faruk and S.J. Savory, "Digital signal processing for coherent 
% transceivers employing multilevel formats," 
% J. Lightwave Technol. 35, 1125 (2017) doi: 10.1109/JLT.2017.2662319.
% 
% The output vectors have unit power, i.e. mean(v1.^2) = mean(v2.^2) = 1.
% v1 and v2 are orthogonal, i.e. mean(v1.*v2) = 0;
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
% v1                first orthonormalised vector [real vector]
%               
%                       v1 is a unit vector, i.e.
%
%                       mean(v1.^2) = 1
%
% v2                second orthonormalised vector [real vector]
%
%                       v2 is orthogonal to v1 and is also a unit vector, 
%                       i.e.
%
%                       mean(v2.^2) = 1
%                       mean(v1.*v2) = 0
%
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

u1 = u1/sqrt(mean(u1.^2));
u2 = u2/sqrt(mean(u2.^2));

aa = mean(u1.*u2);

l11 = 0.5*(1/sqrt(1 + aa) + 1/sqrt(1 - aa));
l12 = 0.5*(1/sqrt(1 + aa) - 1/sqrt(1 - aa));

v1 = l11*u1 + l12*u2;
v2 = l12*u1 + l11*u2;

end