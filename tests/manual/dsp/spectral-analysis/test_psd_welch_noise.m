% -------------------------------------------------------------------------
% Test of white noise signal generation and retrieval of the noise
% PSD with the Welch periodogram method
%
% 2021-05-19
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
% Uniformly distributed random sequence
% -------------------------------------------------------------------------
nperbin = 1000;
% Expected number of samples per bin in the histogram.
nbins = 101;
% Number of bins in the histogram.
nsamples = nperbin*nbins;
% Number of samples.

a = -0.5;
b = 0.5;
noise_uniform = a + (b - a)*rand(1,nsamples);
% We generate noise uniformly distributed in the interval [-0.5 0.5].


variance_uniform_check = var(noise_uniform);
% Estimate variance.
mean_uniform_check = mean(noise_uniform);
% Estimate mean.

fprintf('\n\n\n%s\n','Uniformly distributed random sequence')
fprintf('%s\n\n','=====================================')

fprintf('%s%3.6f\n','Mean= ',mean_uniform_check);
fprintf('%s%3.6f\n','Expected mean= ',(a+b)/2);
fprintf('%s%3.6f\n','Variance= ',variance_uniform_check);
fprintf('%s%3.6f\n\n','Expected variance= ',(b-a)^2/12);

[r_uniform,k] = local_xcorr(noise_uniform);
% Calculate autocorrelation of the sequence.
% We check the whiteness of the generated noise sequence.


figure('Name','uniformly distributed random sequence - autocorrelation')
plot(k,r_uniform,'b-')
xlabel('lag k')
ylabel('autocorrelation r(k)')


figure('Name','uniformly distributed random sequence - histogram')
histogram(noise_uniform,101)
xlabel('n')
ylabel('count')
hline(nperbin,'r--')
% Plot histogram.
% We verify that the variable is uniformly distributed.

block_size = 100;
overlap_samples = 50;
fft_length = 128;
window_type = 'hann';
[psd_uniform,norm_freq] = psd_welch(noise_uniform,block_size,overlap_samples,window_type,fft_length,'two_sided','psd');
% Calculate power spectral density using Welch periodogram methof.

psd_double_sided_uniform = variance_uniform_check/2/pi;
% Expected double-sided noise spectral density, from the estimated noise
% variance.

figure('Name','uniformly distributed random sequence - periodogram')
h1 = plot(norm_freq/pi,psd_uniform);
ylim([0 1.1*max(psd_uniform)]);
hold on
h2 = hline(psd_double_sided_uniform,'r--');
xlabel('normalised angular frequency (\times \pi rad/sample)')
ylabel('power spectral density (power per rad/sample)')
legend([h1,h2],{'periodogram','expected'},'Location','SouthEast')



%%
% -------------------------------------------------------------------------
% Normally distributed random sequence
% -------------------------------------------------------------------------
variance_specified = 0.5;
mean_specified = 0;
noise_normal = mean_specified + sqrt(variance_specified)*randn(1,nsamples);
% We generate Gaussian noise with specified mean and variance


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
histogram(noise_normal,101)
xlabel('n')
ylabel('count')
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
legend([h1,h2],{'periodogram','expected'},'Location','SouthEast')


noise_fft_normal = fft(noise_normal)/nsamples;
% FFT of normal noise.
noise_fft_normal_tot = sum(abs(noise_fft_normal).^2);
% Total noise power calculated in the frequency domain.

fprintf('%s%3.6f\n\n','Noise power from FFT domain= ',noise_fft_normal_tot);







