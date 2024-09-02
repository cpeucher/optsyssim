function S = char_opt_stokes(sig)
% Characterization of instantaneous Stokes parameters
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the instantaneous Stokes parameters of a light
% beam using the conventional method:
% S0 is calculated from the sum of the intensities mesured after a
% horizontal and a vertical linear polarizer
% S1 is calculated from the difference between the intensities measured
% after a horizontal and a vertical linear polarizer
% S2 is calculated from the difference between the intensities measured
% after a linear polarizer at +45 degrees and a linear polarizer at +135
% degrees
% S3 is calculated as the difference between the intensities measured after
% a right hand circular and left hand circular polarizer. These
% polarization components can be obtained after a quarter wave plate
% followed by a linear polarizer at 45 degrees (R) or ar 135 degrees (L).
% The Stokes parameters are "instantaneous" in the sense that no averaging
% is performed in this function. If a time-varying signal is input, then
% the Stokes parameters are calculated as a function of time.
% This also implies that the detector has infinite bandwidth.
% Some time averaging will be necessary in order to calculate the Stokes
% vector in case the input signal is not fully polarized. This time
% averaging needs to be performed after calling this function, for instance
% using the cal_stokes_integrate function.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% S = char_opt_stokes(sig);
% S = S./S(1,:);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% S                 instantaneous Stokes parameters [real matrix]
%                       
%                       S is an 4xN matrix, where N is the number of
%                       samples:
%                       S(1,:) contains samples of Stokes parameter S0
%                       S(2,:) contains samples of Stokes parameter S1
%                       S(3,:) contains samples of Stokes parameter S2
%                       S(4,:) contains samples of Stokes parameter S3
%
%                       The Stokes parameters can be normalized according
%                       to: 
%                       S = S./S(1,:);
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
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig_000 = pol_polarizer(sig,0,0,Inf);
% Linear polarizer // x

sig_090 = pol_polarizer(sig,pi/2,0,Inf);
% Linear polarizer // y

sig_045 = pol_polarizer(sig,pi/4,0,Inf);
% Linear polarizer at 45 degrees with respect to x

sig_135 = pol_polarizer(sig,3*pi/4,0,Inf);
% Linear polarizer at 135 degrees with respect to x

sig_R = pol_polarizer(pol_retarder(sig,0,pi/2),pi/4,0,Inf);
% Quarter wave plate followed by linear polarizer at 45 degrees with
% respect to -x
% Equivalent to right-handed circular polarizer

sig_L = pol_polarizer(pol_retarder(sig,0,pi/2),3*pi/4,0,Inf);
% Quarter wave plate followed by linear polarizer at 135 degrees with
% respect to -x
% Equivalent to Left-handed circular polarizer

I_000 = abs(sig_000.x).^2 + abs(sig_000.y).^2;
I_090 = abs(sig_090.x).^2 + abs(sig_090.y).^2;
I_045 = abs(sig_045.x).^2 + abs(sig_045.y).^2;
I_135 = abs(sig_135.x).^2 + abs(sig_135.y).^2;
I_R = abs(sig_R.x).^2 + abs(sig_R.y).^2;
I_L = abs(sig_L.x).^2 + abs(sig_L.y).^2;

S(1,:) = I_000 + I_090;
S(2,:) = I_000 - I_090;
S(3,:) = I_045 - I_135;
S(4,:) = I_R - I_L;
% Instantaneous Stokes parameters

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------