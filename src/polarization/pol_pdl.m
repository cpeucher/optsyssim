function sig = pol_pdl(sig,loss,pdl)
% polarization dependent loss
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements polarization dependent loss (PDL) to a signal. 
% The -x and -y polarizations of the input signal are attenuated by 
% different amounts.
% Loss and PDL can vary with time.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% il = 0;       % Insertion loss, in dB (positive number)
% pdl = 3;      % Polarization dependent loss, in dB (positive number)
% sig = pol_pdl(sig,insertion_loss,pdl); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% loss              attenuation experienced by the -x polarization, in dB
%                       [real vector]
%
%                       For attenuation a positive number should be 
%                       provided. Otherwise the function will actually 
%                       amplify the signal.
%                       
%                       For static loss, loss can be a real scalar
%
% pdl               extra attenuation experienced by the -y polarization,
%                       in dB [real vector]
%
%                       If this number is positive, the -y polarization 
%                       experiences more attenuation than the -x 
%                       polarization
%
%                       For static PDL, pdl can be a real scalar
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

attx = loss;
atty = loss + pdl;

sig.x = sig.x.*10.^(-attx/20);
sig.y = sig.y.*10.^(-atty/20);

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------