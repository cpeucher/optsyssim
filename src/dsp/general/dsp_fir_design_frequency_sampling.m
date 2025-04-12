function h = dsp_fir_design_frequency_sampling(freq,tf,fs,ntaps,display)
% FIR filter coefficients determination by frequency sampling
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function synthesize a digital FIR filter by the frequency sampling
% method. The analogue transfer function is provided. It is them sampled by
% linear interpolation with a specified number of samples, which
% corresponds to the number of taps of the filter.
% No check on whether the definition of the provided analogue response
% matches the sampling frequency is performed. We are targeting an educated
% audience.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% roll_off = 0.1; 
% % Raised-cosine roll-off factor.
% ts = 1/symbol_rate; 
% % Symbol rate, in baud.
% nos = 4;
% % Oversampling factor. Even number.
% fsa = nos*symbol_rate;
% % Sampling frequency.
% nlengthsymbs = 16;
% % Filter depth, in terms of number of symbols.
% % 16 means 8 symbols on each side.
% % Even number.
% ntaps = nlengthsymbs*nos + 1;
% % FIR filter length. Ensure it is an odd number.
% freq_sampling = linspace(0,fsa/2,1000);
% spectrum_sampling = calc_rc_spectrum(freq_sampling,ts,roll_off);
% % Analog spectrum of raise-cosine filter,
% h = dsp_fir_design_frequency_sampling(freq_sampling,spectrum_sampling,fsa,ntaps,1);
% % Tap coefficients. 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              analogue frequency vector [real vector]
%
%                       Contains the positive frequencies over which the
%                       analogue transfer function is defined.
%
% tf                analogue transfer function [real vector]
%
%                       Contains the values of the analogue filter one
%                       wishes to implement digitally.
%
% fs                sampling frequency, in Hz [real scalar]
%
% ntaps             number of taps of the FIR filter [integer]
%
%                       This corresponds to the number of samples in the
%                       frequency domain or the length of the filter.
%
% display           plot the tap coefficients and filter response [boolean]
%
%                       display = 0,1.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% h                 tap coefficient [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


if nargin == 4
    display = 0;
    % Do not display the taps and transfer functions if no display argument
    % is passed to the function.
end





fnorm = freq/fs;
% Normalized frequency axis

fsamp_norm = [0:1:(ntaps-1)/2]/ntaps;
% Normalized sampling frequencies over [0,1/2] interval

fsamp = fsamp_norm*fs; 
% Corresponding frequencies

Hsamp = interp1(freq,tf,fsamp);
% Sample the analogue transfer function by interpolation

Hsamp = [Hsamp,fliplr(conj(Hsamp(2:end)))];
% Reconstitute samples of full frequency response over [0,1[ normalized 
% frequency interval

h = ifft(Hsamp);
h = fftshift(h);
% Filter tap coefficients

if display    
    figure('Name','tap coefficients')
    stem([0:ntaps - 1],h)
    xlabel('n')
    ylabel('tap coefficient')    
end

wh = [2*pi*fnorm];
% Normalized angular frequency over [0,pi] interval
kk = [0:1:ntaps-1].';
H = h*exp(-1i*kk*wh);
% Calculate the frequency response over [0,1/2] normalized frequency
% interval

if display
    figure('Name','analog and digital filter responses')
    plot(2*fnorm,abs(tf).^2,'b-')
    hold on
    plot(2*fsamp_norm,abs(Hsamp(1:1+(ntaps-1)/2)).^2,'ro')
    plot(2*fnorm,abs(H).^2,'r-')
    xlabel('normalized angular frequency (x\pi radians-per-sample)')
    ylabel('|H|^2')
    legend('desired analog filter response','samples','digital filter response')
end



end