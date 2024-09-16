function sig = opt_dispersion(sig,params,z)
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
% params_dispersion.dispersion = 17;           % Dispersion in ps/nm/km
% params_dispersion.dispersion_slope = 0.058;  % Dispersion slope in ps/nm^2/km
% params_dispersion.dispersion_curvature = 0;  % Dispersion curvature, in ps/nm^3/km
% params_dispersion.dispersion_spec_frequency = reference_frequency;
% z = 1e3; 
% sig = opt_dispersion(sig,params_dispersion,z); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% params            dispersion parameters [structure]
%
%                       params.dispersion
%                           dispersion, in ps/nm/km [real scalar]
%
%                       params.dispersion_slope
%                           dispersion slope, in ps/nm^2/km [real scalar]
% 
%                       params.dispersion_curvature
%                           dispersion curvature, in ps/nm^3/km 
%                           [real scalar]
%
%                       params.dispersion_spec_frequency
%                           frequency at which the params.dispersion value 
%                           is specified, in Hz [real scalar]
%
% z                 distance at which the dispersed signal is calculated,  
%                       in m [real scalar]
% 
%                       When z = 1000, then params.dispersion corresponds 
%                       to the accumulated dispersion in ps/nm.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array 
global reference_frequency

beta = dispersion_conv_d_beta([params.dispersion, params.dispersion_slope, params.dispersion_curvature],'to_beta','eng','SI',params.dispersion_spec_frequency);
% Convert the dispersion and dispersion slope to beta_2 and beta_3 

ref_freq = params.dispersion_spec_frequency - reference_frequency;
% Dispersion specification frequency relative to the reference frequency.

tf = exp(-1i*(0.5*beta(2)*(2*pi*(frequency_array - ref_freq)).^2 +...
    beta(3)/6*(2*pi*(frequency_array - ref_freq)).^3)*z);
% All-pass filter due to dispersion.

sig.x = ifft(fft(sig.x).* fftshift(tf));
sig.y = ifft(fft(sig.y).* fftshift(tf));
% Apply the filter to the input signal.

end