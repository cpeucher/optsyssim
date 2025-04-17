% -------------------------------------------------------------------------
% Comparison of average power estimations in the time and frequency domains
% for various types of signals:
%
% 1. Continuous (DC) signal
% 2. Sinusoidal signal (AC + DC)
% 3. White noise signal
% 4. Noisy sinusoidal signal (AC + DC + white noise)
%
% Purpose is to illustrate equality between power estimated in the time and
% frequency domain.
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Clean up
% -------------------------------------------------------------------------
clear all
close all

% -------------------------------------------------------------------------
% Global
% -------------------------------------------------------------------------
nsamples = 2^22;
% Number of samples.

fprintf('\n\n\n%s','Time-frequency power estimations')
fprintf('\n%s','================================')
fprintf('\n%s','================================')

% -------------------------------------------------------------------------
% First, a continuous signal
% -------------------------------------------------------------------------
fprintf('\n\n\n%s','Test of continuous signal')
fprintf('\n%s','=========================')


x = rand(1)*ones(1,nsamples);
% Signal, defined in the time domain.

mean_power_time = sum(abs(x).^2)/nsamples
% Mean power, calculating in the time domain.

X = fft(x)/nsamples;
% Signal spectrum (two-sided).

mean_power_frequency = sum(abs(X).^2)
% Mean power calculated in the frequency domain.

figure('Name','continuous signal')
subplot(2,1,1)
plot(x)
subplot(2,1,2)
plot(10*log10(fftshift(abs(X).^2) + eps))


% -------------------------------------------------------------------------
% Second, a sinusoidal signal
% -------------------------------------------------------------------------

fprintf('\n\n\n%s','Test of sinusoidal signal')
fprintf('\n%s','=========================')

f = 3/128;
xdc = 3;
xac = 0.23;
x = xdc + xac*cos(2*pi*f*[0:1:nsamples-1]);
% Signal, defined in the time domain.
% DC power is 3^2 = 9 (9.54 dB)
% AC power is 0.23^2/2 = 0.02645 (-15.78 dB)



signal_variance = var(x)
% Should correspond to the ac power

signal_ac_power  = sum((x - xdc).^2)/nsamples
% Should be more or less the same as signal_variance, apart from the fact
% that the variance returns the sample variance, 
% i.e. divided by nsamples -  1

signal_dc_power_expected = xdc^2
% Expected DC power.

signal_ac_power_expected = xac^2/2
% Expected 

X = fft(x)/nsamples;
% Signal spectrum (two-sided)
% In the two-sided spectrum:
% - DC line should be at 9.54 dB
% - AC lines should be at -18.78 dB each (total at -15.78 dB)

X1 = X(1:nsamples/2);
X1(2:end) = sqrt(2)*X1(2:end);
% Signal spectrum (one-sided).
% - DC line should be at 9.54 dB
% - AC line should be at -15.78 dB 

figure('Name','sinusoidal signal')
subplot(3,1,1)
plot(x)
subplot(3,1,2)
plot(10*log10(fftshift(abs(X).^2)))
subplot(3,1,3)
plot(10*log10(abs(X1).^2))

mean_power_time = sum(x.^2)/nsamples
% Mean power, calculating in the time domain.
% Should be 9.02645

mean_power_frequency_one_sided = sum(abs(X).^2)
mean_power_frequency_two_sided = sum(abs(X1).^2)
% Hopefully should be the same as above.


% -------------------------------------------------------------------------
% Then a random signal
% -------------------------------------------------------------------------
fprintf('\n\n\n%s','Test of random signal')
fprintf('\n%s','=====================')

sigma2 = 0.1
% Wished noise variance.
xn = sqrt(sigma2)*randn(1,nsamples);

mean_signal_time = mean(xn)

noise_signal_variance = var(xn)
% Check variance of noise realisation.

noise_signal_variance_check = sum(xn.^2 - mean(xn)^2)/(nsamples - 1)
% Check the variance estimator, just for fun.

Xn = fft(xn)/nsamples;
% FFT of noise signal.

figure('Name','random signal')
subplot(2,1,1)
plot(xn)
subplot(2,1,2)
plot(10*log10(fftshift(abs(Xn).^2)))

mean_power_frequency = sum(abs(Xn).^2)


% -------------------------------------------------------------------------
% And finally a noisy sinusoidal signal
% -------------------------------------------------------------------------

fprintf('\n\n\n%s','Test of noisy sinusoidal signal')
fprintf('\n%s','===============================')


% Noise signal.
xt= x + xn;
% Noisy signal, defined in the time domain.
% DC signal power is 3^2 = 9 (9.54 dB)
% AC signal power is 0.23^2/2 = 0.02645 (-15.78 dB)
% Noise power is 0.1 (-10 dB)

noise_signal_variance = var(xn)
% Check variance of noise realisation.

signal_ac_power = var(x)

signal_dc_power = xdc^2

mean_power_time = sum(xt.^2)/nsamples
% Mean power, calculating in the time domain.
% Should be 9.02645

total_signal_variance = var(xt)
% Should be sum of ac power and noise variance

sum_ac_power_noise_var_expected = sigma2 + signal_ac_power_expected
sum_ac_power_noise_var_estimated = noise_signal_variance + signal_ac_power 

Xt = fft(xt)/nsamples;
% FFT of noise signal.

mean_power_frequency = sum(abs(Xt).^2)

