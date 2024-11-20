% -------------------------------------------------------------------------
% Test of Welch periodogram function psd_welch.m
%
% We test our own implementation of Welch periodogram for power spectral
% density (PSD) estimation.
% We reproduce the examples:
% Use Matlab Function pwelch to Find Power Spectral Density – or Do It 
% Yourself by Neil Robertson, January 13, 2019
% https://www.dsprelated.com/showarticle/1221.php
%
% 2021-04-16
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all


% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Reset random number generator
% -------------------------------------------------------------------------
rng default

% -------------------------------------------------------------------------
% Signal generation
% -------------------------------------------------------------------------
fs = 4000;
% Sample rate, in Hz.
f0 = 500;
% Sinusoidal signal frequency, in Hz.
A = sqrt(2);
% Peak amplitude of the sinusoid in units of V
% The rms value is A/sqrt(2) V.
% The average power is A^2/2 = 1 W into 1 ohm
noise_variance = 0.1^2;
% Gaussian noise variance.
% The average noise power is equal to noiseVar.
nsamples = 1024*8;
% Number of samples in the signal.
time_array = (0:1:nsamples-1)/fs;
% Time axis.
df = fs/nsamples;
% Frequency interval between consecutive samples.
frequency_array = (-nsamples/2:nsamples/2-1)*df;
% Frequency axis.
x = A*cos(2*pi*f0*time_array) + sqrt(noise_variance)* randn(size(time_array));
% Noisy sinusoidal signal.

% In the double sideband spectrum, the amplitude of each frequency tone is
% A/2.
% Its power is A^2/4, i.e. 0.5 W, or -3 dBW.


% -------------------------------------------------------------------------
% Waveform
% -------------------------------------------------------------------------
figure('Name','Noisy signal waveform')
plot(time_array,x,'b');
xlabel('time (s)');
ylabel('voltage (V)')



% -------------------------------------------------------------------------
% Fourier transform of the waveform
% -------------------------------------------------------------------------
X = fftshift(fft(x))/nsamples;

fprintf('\n\n%s\t%3.2f\t%s\n','Signal average power= ',A^2/2,'W');
fprintf('%s\t%3.2f\t%s\n','Signal average power= ',10*log10(A^2/2),'dBW');

fprintf('\n\n%s\t%3.2f\t%s\n','Expected power of each tone in double sideband spectrum= ',A^2/4,'W');
fprintf('%s\t%3.2f\t%s\n','Expected power of each tone in double sideband spectrum= ',10*log10(A^2/4),'dBW');

fprintf('\n\n%s\t%3.2f\t%s\n','Retrieved peak power in double sideband fft spectrum= ',max(abs(X).^2),'W');
fprintf('%s\t%3.2f\t%s\n','Retrieved peak power in double sideband fft spectrum= ',10*log10(max(abs(X).^2)),'dBW');


figure('Name','spectrum of noisy sinusoidal signal using fft')
plot(frequency_array,10*log10(abs(X).^2),'b-')
xlabel('frequency (Hz)');
ylabel('power (dBW)');
ylim([-90 0]);

% Summary
% -------
% We are able to retrieve the power of the sinusoidal signal in the fft
% spectrum...


% -------------------------------------------------------------------------
% Power spectrum using psd_welch
% Single block
% Rectangular Window
% No overlap
%--------------------------------------------------------------------------
% We expect something pretty close to what we obtained by simply
% calculating the fft in the previous stage...

block_size = nsamples;
overlap_samples = 0;
fft_length = nsamples;
window_type = 'rectangular';

[psd_welch_pow_1,freq_1] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','power',fs);


fprintf('\n\n%s\t%3.2f\t%s','Peak power in single-sided psd_welch power spectrum= ',max(psd_welch_pow_1),'W');
fprintf('\n%s\t%3.2f\t%s\n','Peak power in single-sided psd_welch power spectrum= ',10*log10(max(psd_welch_pow_1)),'dBW');


figure('Name','psd_welch - power - single rectangular window');
plot(freq_1,10*log10(psd_welch_pow_1),'b-');
xlabel('frequency (Hz)')
ylabel('power (dBW)')
ylim([-65 3]);
grid on

% Summary
% -------
% This works. We indeed retrieve the power of the signal in the
% single-sided periodogram using a single window...


% -------------------------------------------------------------------------
% Power spectrum using psd_welch
% Block size = nsamples/8
% Hann window
% Overlap = 1/2 block size
% -------------------------------------------------------------------------
block_size = 1024;
overlap_samples = block_size/2;
fft_length = block_size;
window_type = 'hann';


[psd_welch_pow_2,freq_2] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','power',fs);

fprintf('\n\n%s\t%3.2f\t%s','Peak power in single-sided psd_welch power spectrum= ',max(psd_welch_pow_2),'W');
fprintf('\n%s\t%3.2f\t%s\n','Peak power in single-sided psd_welch power spectrum= ',10*log10(max(psd_welch_pow_2)),'dBW');

figure('Name','psd_welch - power - overlapping Hann windows');
plot(freq_2,10*log10(psd_welch_pow_2),'b-');
xlabel('frequency (Hz)')
ylabel('power (dBW)')
%ylim([-65 3]);
grid on

% Summary
% -------
% This still works. We still recover the power of the signal


% -------------------------------------------------------------------------
% Power spectral density using psd_welch
% Block size = nsamples/8
% Hann window
% Overlap = 1/2 block size
% -------------------------------------------------------------------------

[psd_welch_psd_3,freq_3] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd',fs);

fprintf('\n\n%s\t%3.2f\t%s','Peak power in single-sided psd_welch psd spectrum= ',max(psd_welch_psd_3),'W');
fprintf('\n%s\t%3.2f\t%s\n','Peak power in single-sided psd_welch psd spectrum= ',10*log10(max(psd_welch_psd_3)),'dBW');

% The single-sided noise spectral density is 
N0 = 2*noise_variance/fs;
fprintf('\n\n%s\t%3.2f\t%s\n','Calculated single-sided noise spectral density = ',10*log10(N0),'dB/Hz');

figure('Name','psd_welch - psd - overlapping Hann windows');
plot(freq_3,10*log10(psd_welch_psd_3),'b-');
xlabel('frequency (Hz)')
ylabel('psd (dB/Hz)')
hline(10*log10(N0),'r--')
%ylim([-65 3]);
grid on

% Summary
% -------
% We retrieve the correct power spectral density of the additive (filtered)
% white Gaussian noise.