function sig = mod_polm(sig,drive,vpi,loss)
% Ideal polarisation modulator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function models an ideal polarisation modulator. The modulator is 
% ideal in the sense that it does not create any phase modulation on top of
% the wished polarisation modulation. Such a phase modulation would arise
% in electro-optic polarisation modulators unless the sum of the
% electro-optic induced refractive index changes along the -x and -y
% directions cancel out (can be achieved in GaAs, for instance). Here the
% polarisation is switched between -x (for V=0) and -y (for V=vpi). 
% No bandwidth limitation is introduced. Therefore, if needed, such 
% limitation needs to be modelled by low pass filtering the driving
% signals.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% vpi = 1.0;
% loss = 0;
% sig = mod_polm(sig,nrz_data_sig,vpi,loss); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       The signal is assumed to be polarised along -x for
%                       standard operation. This is currently not checked.
%
% drive             electrical driving signal applied to the modulator
%                       [real vector]
%
% vpi               half wave voltage of the modulator [real scalar]
%
% loss              modulator loss in dB (positive number) [real scalar]
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
% Christophe Peucheret (christophe.peucheret@univ-rennes1.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.y = 10^(-loss/20)*sig.x.*sin(0.5*pi*(drive/vpi));
sig.x = 10^(-loss/20)*sig.x.*cos(0.5*pi*(drive/vpi));
% Keep the lines above in this order to avoid applying the modulation
% twice to sig.x!


end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------