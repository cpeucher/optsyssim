% -------------------------------------------------------------------------
% Test of white Gaussian noise generation and retrieval of the noise
% PSD with the Welch periodogram method
%
% 2022-06-24
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
% Normally distributed random sequence
% -------------------------------------------------------------------------
nsamples = 2^17;
variance_specified = 0.5;
mean_specified = 0;
noise_normal = mean_specified + sqrt(variance_specified)*randn(1,nsamples);
% We generate Gaussian noise with specified mean and variance

x = [-5:0.01:5]*sqrt(variance_specified) + mean_specified;
% Axis for plotting theoretical pdf.

pdf_theoretical = 1/sqrt(2*pi*variance_specified)*exp(-(x - mean_specified).^2/2/variance_specified);
% Theoretical pdf.

variance_normal_check = var(noise_normal);
% Estimate variance.
mean_normal_check = mean(noise_normal);
% Estimate mean.

fprintf('\n\n\n%s\n','Normally distributed random sequence')
fprintf('%s\n\n','====================================')

fprintf('%s%3.6f\n','Mean= ',mean_normal_check);
fprintf('%s%3.6f\n','Expected mean= ',mean_specified);
fprintf('%s%3.6f\n','Variance= ',variance_normal_check);
fprintf('%s%3.6f\n\n','Expected variance= ',variance_specified);

[r_normal,k] = local_xcorr(noise_normal);
% Calculate autocorrelation of the sequence.
% We check the whiteness of the generated noise sequence.


figure('Name','normally distributed random sequence - autocorrelation')
plot(k,r_normal,'b-')
xlabel('lag k')
ylabel('autocorrelation r(k)')


figure('Name','normally distributed random sequence - histogram')
histogram(noise_normal,101,'Normalization','pdf')
hold on
plot(x,pdf_theoretical,'r')
xlabel('n')
ylabel('count')
legend('histogram','theoretical pdf')
% Plot histogram.
% We verify that the variable is normally distributed.

block_size = 100;
overlap_samples = 50;
fft_length = 128;
window_type = 'hann';
[psd_normal,norm_freq] = psd_welch(noise_normal,block_size,overlap_samples,window_type,fft_length,'two_sided','psd');
% Calculate power spectral density using Welch periodogram method.

psd_double_sided_normal = variance_normal_check/2/pi;
% Expected double-sided noise spectral density, from the estimated noise
% variance.

figure('Name','normally distributed random sequence - periodogram')
h1 = plot(norm_freq/pi,psd_normal);
hold on
h2 = hline(psd_double_sided_normal,'r--');
xlabel('normalised angular frequency (\times \pi rad/sample)')
ylabel('power spectral density (power per rad/sample)')
ylim([0 1.1*max(psd_normal)]);
legend([h1,h2],{'periodogram','theoretical psd'},'Location','SouthEast')


noise_fft_normal = fft(noise_normal)/nsamples;
% FFT of normal noise.
noise_fft_normal_tot = sum(abs(noise_fft_normal).^2);
% Total noise power calculated in the frequency domain.

fprintf('%s%3.6f\n\n','Noise power from FFT domain= ',noise_fft_normal_tot);







