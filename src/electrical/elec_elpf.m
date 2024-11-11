function sig = elec_elpf(sig,params)
% Electrical low-pass filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a standard electrical low pass filter to an input
% electrical signal.
% For backward compatibility.
% Equivalent to:
% tf = elec_tf_elpf(params_elpf,frequency_array);
% sig = elec_filter(sig,tf);
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
% params_elpf.order = 4;
% params_elpf.f3dB = 0.75*symbol_rate;
% params_elpf.roll_off = 0.6;            % for 'raised_cosine' and 'root_raised_cosine'
% params_elpf.symbol_rate = symbol_rate; % for 'raised_cosine' and 'root_raised_cosine'
% sig = elec_elpf(sig,params_elpf); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% params            parameters of the low-pass filter [structure]
%
%                       See the description of the elec_tf_elpf function.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               low pass filtered electrical signal [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array   relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array

tf = elec_tf_elpf(params,frequency_array);
% Calculate transfer function of the low pass filter

sig = real(ifft(fft(sig).* fftshift(tf)));
% Applies the filter to the input signal

end
