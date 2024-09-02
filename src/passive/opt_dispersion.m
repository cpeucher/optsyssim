function sig = opt_dispersion(sig,params)
% Linear and lossless dispersive element
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a dispersion value (possibly including dispersion
% slope) to the input signal. It can be used for e.g. investigation of 
% dispersion tolerance where the real behaviour of a fibre (including loss 
% and nonlinearities) is of no interest, but the focus is exlusively on
% dispersion.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_dispersion.dispersion = 17;           % Dispersion in ps/nm
% params_dispersion.dispersion_slope = 0.058;  % Dispersion slope in ps/nm^2.
% params_dispersion.dispersion_spec_frequency = reference_frequency;
% sig = opt_dispersion(sig,params_dispersion); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       This is it.'
%
% params            dispersion values [structure]
%
%                       params.dispersion
%                           value of dispersion to apply to the signal, 
%                           in ps/nm [real scalar]
%
%                       params.dispersion_slope
%                           dispersion slope, in ps/nm^2 [real scalar]
%
%                       params.dispersion_spec_frequency
%                           frequency at which the params.dispersion value 
%                           is specified, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%                       This is it.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
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

global frequency_array 
global reference_frequency


length = 1000;
% The fibre length is arbitrarily set to 1 km.
% Then the dispersion value given in ps/nm corresponds to the same
% dispersion per unit length value in ps/nm/km.

dispersion_si = params.dispersion*1.0e-6;
% Convert dispersion from ps/nm/km to s/m^2.
dispersion_slope_si = params.dispersion_slope*1.0e3;
% Convert dispersion slope from ps/nm^2/km to s/m^3.

beta = dispersion_conv_d_beta([dispersion_si dispersion_slope_si],'to_beta','SI','SI',params.dispersion_spec_frequency);
% Convert the dispersion and dispersion slope to to beta_2 and beta_3 
% (elements beta(2) and beta(3), respectively).

ref_freq = params.dispersion_spec_frequency - reference_frequency;
% Dispersion specification frequency relative to the reference frequency.

tf = exp(-1i*(0.5*beta(2).*(2*pi*(frequency_array - ref_freq)).^2+beta(3)/6*(2*pi*(frequency_array - ref_freq)).^3)*length);
% All-pass filter due to dispersion.

sig.x = ifft(fft(sig.x).* fftshift(tf));
sig.y = ifft(fft(sig.y).* fftshift(tf));
% Apply the filter to the input signal.



end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------