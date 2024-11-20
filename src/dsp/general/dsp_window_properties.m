function [nenb,window_spectrum,ang_freq_norm] = dsp_window_properties(window,fft_length)
% Calculate some properties of DSP windows
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates some properties of DSP windows.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [nenb,spectrum,norm_ang_freq] = dsp_window_properties(w,fft_length);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% window            window as generate by the dsp_window function 
%                       [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% nenb              window normalised equivalent noise bandwidth (NENB),
%                       in frequency bins [real scalar]
%                       
%                       Should take values:
%                           Hann window:            1.5
%
%                       The equivalent noise bandwidth ENB is equal to:
%                            ENB = NENB * fres
%                            where fres is the frequency resolution 
%                           (fs / number of bins)
%
% window_spectrum    magnitude of the window spectrum, in dB [real vector]
% 
% ang_freq_norm      normalised angular frequency, in radian/sample 
%                       [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

window_spectrum = abs(fftshift(fft(window,fft_length))).^2;
window_spectrum = 10*log10(window_spectrum/max(window_spectrum));
% Magnitude of the window spectrum, in dB

ang_freq_norm = (-1:2/fft_length:1 - 2/fft_length)*pi;
% Normalised angular frequency in radian/sample

nenb = length(window)*sum(abs(window).^2)/abs(sum(window))^2;
% Normalised noise equivalent bandwidth, in frequency bins


end