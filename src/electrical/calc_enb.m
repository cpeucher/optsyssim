function enb = calc_enb(type,tf,tf0,df)
% Equivalent noise bandwidth of optical or electrical filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the noise equivalent bandwidth of electrical and
% optical filters.
% The transfer function is calculated outside the function and is provided
% as an input.
% The transfer function needs to be defined over a frequency grid with
% fixed step df.
% One effectively distinguishes between optical and electrical low-pass
% filters through the 'bandpass' vs 'lowpass' parameter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tf = elec_tf_elpf(params_elpf,frequency_array);
% enb = calc_enb_tf('lowpass',tf,max(abs(tf)),df);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% type              type of filter [string]
%
%                       type = 'bandpass'
%                           bandpass filter (typically optical filter)
%
%                       type = 'lowpass';
%                           lowpass filter (typically electrical
%                               post-detection filter)
%
% tf                filter transfer function specified on a fixed-step
%                       frequency grid [complex vector]
%
%                       The step of the frequency grid is df.
%                       The vector tf is specified in increasing frequency
%                       order.
%
% tf0               reference level
%                       Typically, for well behaved electrical low-pass
%                       filter or optical bandpass filters whose transfer
%                       function is centered on the reference frequency of
%                       the baseband transformation:
%                       tf = tf(frequency = 0),
%
%                       If tf is defined over the global frequency_array
%                           grid:
%
%                           tf_fft = fftshift(tf);   % sort tf to FFT order
%                           tf0 = tf_fft(1);         % value at dc (elec) /
%                                                    % centre freq. (opt)
%
%                       In case the max of |tf| is not at 0 frequency,
%                           one could also choose:
%
%                           tf0 = max(abs(tf));
%
%                       For many filters, whose transfer function peaks at
%                       dc (elec) / reference_frequency (opt), the
%                       distinction is irrelevant.
%
% df                step of the frequency grid over which the transfer
%                       function is specified, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% enb               equivalent noise bandwidth, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(type,'lowpass')
    
    enb = num_int1d_simpson(abs(tf).^2,df)/(2*abs(tf0)^2);
    % Noise equivalent bandwidth of the electrical low-pass filter.    
    
elseif strcmp(type,'bandpass')
    
    enb = num_int1d_simpson(abs(tf).^2,df)/abs(tf0)^2;
    % Noise equivalent bandwidth of the optical bandpass filter.
    
end

end