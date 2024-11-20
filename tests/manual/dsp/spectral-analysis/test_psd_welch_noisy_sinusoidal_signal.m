% -------------------------------------------------------------------------
% Test of white noise addition to a real signal and retrieval of the noise
% PSD with the Welch periodogram method
%
% 2021-04-20
%
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all


% -------------------------------------------------------------------------
% Reset random number generator.
% -------------------------------------------------------------------------
rng default


% -------------------------------------------------------------------------
% Generate sinusoidal signal with white noise
% -------------------------------------------------------------------------
fs = 10000;
% Sample rate, in Hz.
signal_frequency = 100;
% Sinusoidal signal frequency, in Hz.
signal_peak_amplitude = sqrt(2);
% Peak amplitude of the sinusoid in units of U
% The rms value is signal_peak_amplitude/sqrt(2)
% The average power is signal_peak_amplitude^2/2
noise_variance = 0.1;
% Gaussian noise variance.
% The average noise power is equal to noise_variance.
time_array = 0:1/fs:5-1/fs;
% Time array.
sig = signal_peak_amplitude*cos(2*pi*signal_frequency*time_array) + sqrt(noise_variance)* randn(size(time_array));
% Noisy sinusoidal signal.

figure('Name','noisy sinusoidal signal')
plot(time_array,sig,'b-');
xlabel('time (s)');
ylabel('amplitude (U)');

% We recover the total power in the time domain:
estimated_total_power_mean_time = mean(sig.^2);

estimated_noise_power_mean_time = estimated_total_power_mean_time - signal_peak_amplitude^2/2;
% That should be the noise power.

fprintf('\n\n\n%s\n','Analysis of noisy sinusoidal signal')
fprintf('%s\n','===================================')

fprintf('\n%s\n','Signal and additive white Gaussian noise generation')
fprintf('%s','---------------------------------------------------')

fprintf('\n\n%s\t%3.2f\t%s\n','Set signal average power= ',signal_peak_amplitude^2/2,'W');
fprintf('%s\t%3.2f\t%s\n','Set signal average power= ',10*log10(signal_peak_amplitude^2/2),'dBW');

fprintf('\n%s%3.5f%s\n','Total power estimated in the time domain= ',estimated_total_power_mean_time,' W');
fprintf('%s%3.5f%s\n','Noise power estimated in the time domain= ',estimated_noise_power_mean_time,' W');
fprintf('%s%3.5f%s\n','Set noise variance= ',noise_variance,' W');
fprintf('%s%3.5f%s\n','Ratio between set and estimated noise power= ',noise_variance/estimated_noise_power_mean_time,' (-)');
fprintf('%s%3.5f%s\n','Ratio between set and estimated noise power= ',10*log10(noise_variance/estimated_noise_power_mean_time),' dB');



%%
% -------------------------------------------------------------------------
% We plot the signal spectrum. 
% In Matlab it is the DFT divided by the number of samples
% -------------------------------------------------------------------------

nsamples = length(sig);
% Number of samples in the input signal.
df = fs/nsamples;
frequency_array = (-nsamples/2:nsamples/2-1)*df;
% Frequency axis.

signal_spectrum = fftshift(fft(sig))/nsamples;
% Spectrum.

fprintf('\n\n%s\n','FFT spectrum analysis')
fprintf('%s','---------------------')

fprintf('\n\n%s%i','The FFT length is equal to the vector length = ', nsamples)

fprintf('\n\n%s\t%3.2f\t%s\n','Expected power of each tone in double-sided spectrum= ',signal_peak_amplitude^2/4,'W');
fprintf('%s\t%3.2f\t%s\n','Expected power of each tone in double-sided spectrum= ',10*log10(signal_peak_amplitude^2/4),'dBW');


fprintf('\n%s\t%3.2f\t%s\n','Retrieved peak power in double-sided fft spectrum= ',max(abs(signal_spectrum).^2),'W');
fprintf('%s\t%3.2f\t%s\n','Retrieved peak power in double-sided fft spectrum= ',10*log10(max(abs(signal_spectrum).^2)),'dBW');


figure('Name','spectrum of noisy sinusoidal signal - FFT length = vector length')
plot(frequency_array,10*log10(abs(signal_spectrum).^2),'b-')
xlabel('frequency (Hz)');
ylabel('power (dBU)');
ylim([-90 0]);




%%
% -------------------------------------------------------------------------
% We plot the spectrum with an FFT length different than the length of the
% signal vector
% -------------------------------------------------------------------------
fft_length_2 = 2^17;
% FFT length

df_2 = fs/fft_length_2;
frequency_array_2 = (-fft_length_2/2:fft_length_2/2-1)*df_2;
% Frequency axis.
signal_spectrum_2 = fftshift(fft(sig,fft_length_2)) / nsamples;
% Spectrum.
% Observe the normalisation by Nsamples (and not by fft_length_2) !


fprintf('\n\n%s%i','The FFT length is equal to ', fft_length_2)

fprintf('\n\n%s\t%3.2f\t%s\n','Retrieved peak power in double-sided fft spectrum= ',max(abs(signal_spectrum_2).^2),'W');
fprintf('%s\t%3.2f\t%s\n','Retrieved peak power in double-sided fft spectrum= ',10*log10(max(abs(signal_spectrum_2).^2)),'dBW');


figure('Name',['spectrum of noisy sinusoidal signal - FFT length = ' num2str(fft_length_2)]);
plot(frequency_array_2,10*log10(abs(signal_spectrum_2).^2),'b-');
xlabel('frequency (Hz)');
ylabel('power (dBU)');
ylim([-90 0]);

% We should observe the effect of zero padding: windowing of the signal by
% a rectangular window => in the spectral domain convolution of the
% spectrum by a sinc function.
% Indeed this is what we observe.





% -------------------------------------------------------------------------
% Use of Welch method: power
% -------------------------------------------------------------------------
block_size = nsamples;
overlap_samples = 0;
fft_length = nsamples;
window_type = 'rectangular';

block_size = 256*8;
overlap_samples = block_size/2;
fft_length = max(256,2^nextpow2(block_size));
window_type = 'hann';

params_window.type = window_type;
params_window.length = block_size;
params_window.symmetry = 'periodic';
window = dsp_window(params_window);
[nenb,window_spectrum,angular_frequency_normalised] = dsp_window_properties(window,fft_length);
% Window property.

[pwelch_power_1,freq_1] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'one_sided','power',fs);
% [pwelch_power_2,freq_2] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'one_sided','power');

fprintf('\n\n%s\n','Welch periodogram: power mode')
fprintf('%s','-----------------------------')

fprintf('\n\n%s\t%3.2f\t%s\n','Peak power in psd_welch single-sided power spectrum= ',max(pwelch_power_1),'W');
fprintf('%s\t%3.2f\t%s\n','Peak power in psd_welch single-sided power spectrum= ',10*log10(max(pwelch_power_1)),'dBW');


figure('Name','power - psd_welch');
h1 = plot(freq_1,10*log10(pwelch_power_1),'b-');
xlabel('frequency (Hz)')
ylabel('power (dBW)')
ylim([-65 3]);
grid on
h2 = hline(10*log10(signal_peak_amplitude^2/2),'r-');
legend([h1 h2],{'Welch power','expected signal power'})

% figure('Name','Power - psd_welch (normalised)');
% plot(freq_2*2,10*log10(pwelch_power_2),'b-');
% xlabel('normalised frequency (\pi rad/sample)')
% ylabel('power (dBW)')
% ylim([-65 3]);
% grid on



%%
% -------------------------------------------------------------------------
% Calculation of the noise PSD
% -------------------------------------------------------------------------
n0 = 2*noise_variance/fs;
% Expected single-sided noise spectral density in W / Hz.

n0n = 2*noise_variance/2/pi;
% Expected two-sided normalised noise spectral density in W/rad/sample.

fprintf('\n\n%s\n','Welch periodogram: power spectral density')
fprintf('%s\n','-----------------------------------------')


fprintf('\n\n%s\t%3.2f\t%s\n','Expected two-sided noise spectral density: ',10*log10(n0/2),'10*log10(W/Hz)')
fprintf('%s\t%3.2f\t%s\n\n','Expected (normalised) two-sided noise spectral density: ',10*log10(n0n/2),'10*log10(W/rad/sample)')

fprintf('\n%s\t%3.2f\t%s\n','Expected one-sided noise spectral density: ',10*log10(n0),'10*log10(W/Hz)')
fprintf('%s\t%3.2f\t%s\n\n','Expected (normalised) one-sided noise spectral density: ',10*log10(n0n),'10*log10(W/rad/sample)')




% -------------------------------------------------------------------------
% Use of Welch method: power spectral density
% -------------------------------------------------------------------------

[psd_welch_psd_3,freq_3] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'one_sided','psd',fs);
[psd_welch_psd_4,freq_4] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

figure('Name','PSD - psd_welch');
h1 = plot(freq_3,10*log10(psd_welch_psd_3),'b-');
xlabel('frequency (Hz)')
ylabel('PSD (10*log10(W/Hz))')
ylim([min(10*log10(psd_welch_psd_3))-5 max(10*log10(psd_welch_psd_3))+5]);
ylim([-65 3]);
grid on
h2 = hline(10*log10(n0),'r-');
legend([h1 h2],{'Welch psd','expected noise psd'})

figure('Name','PSD - psd_welch (normalised)');
h1 = plot(freq_4/pi,10*log10(psd_welch_psd_4),'b-');
xlabel('normalised frequency (\pi rad/sample)')
ylabel('PSD (10*log10(W/rad/sample))')
ylim([min(10*log10(psd_welch_psd_4))-5 max(10*log10(psd_welch_psd_4))+5]);
grid on
h2 = hline(10*log10(n0n),'r-');
legend([h1 h2],{'Welch psd','expected noise psd'})




