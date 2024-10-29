function sig = elec_filter(sig,tf)
% Electrical filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a calculated transfer function to an electrical
% signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
% params_elpf.order = 4;
% params_elpf.f3dB = 0.75*symbol_rate;
% params_elpf.roll_off = 0.15;
% params_elpf.symbol_rate = symbol_rate;
% params_elpf.hold_time = 1/symbol_rate/2;
% tf = elec_tf_elpf(params_elpf,frequency_array);
% sig = elec_filter(sig,tf);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [complex vector]
% 
% tf                filter transfer function, defined on the same 
%                       frequency grid as the input signal [complex vector]
%
%                       The transfer function is specified in increasing
%                       frequency order.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               filtered electrical signal [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig = ifft(fft(sig).* fftshift(tf));
% Applies the filter to the input signal

end