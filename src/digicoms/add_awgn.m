function symbs = add_awgn(symbs,esn0)
% Addition of complex AWGN to signal samples in the digital domain
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function adds samples of (complex, circularly-symmetric) white
% Gaussian noise to a discrete time signal according to a specified
% signal-to-noise ratio, defined in terms of energy per symbol divided by 
% white noise power spectral density, Es/N0.
% The energy per symbol is calculated from the input signal samples. It is
% therefore understood that the length of the input signal should be
% sufficient for this calculation to make sense (i.e. the signal should
% describe the entire constellation and each symbol of the constellation 
% should be present with approximately the same number of occurences).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% symbs = add_awgn(symbs,esn0_db); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             samples of the signal [real or complex vector]
%                       The signal is assumed to be noiseless.
%
% esn0              target signal-to-noise ratio, in dB [real scalar]
%                       The signal-to-noise ratio is expressed in terms of
%                       energy per symbol divided by the additive white 
%                       Gaussian noise spectral density, Es/N0.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs             samples of the noisy signal [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsymbols = length(symbs);
% Number of symbols

es = mean(abs(symbs).^2);
% Mean symbol energy

n0 = es / 10^(esn0/10);
% Noise spectral density

symbs = symbs + sqrt(n0/2)*(randn(1,nsymbols)-1i*randn(1,nsymbols));
% Noisy signal

end