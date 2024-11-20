% -------------------------------------------------------------------------
% Test of Welch periodogram function psd_welch.m
%
% We test our own implementation of Welch periodogram for power spectral
% density (PSD) estimation.
% we reproduce the examples in the pwelch.m help page:
% https://fr.mathworks.com/help/signal/ref/pwelch.html.
%
% 2021-04-20
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean up
% -------------------------------------------------------------------------
clear all
close all

% ------------------------------------------------------------------------- 
% Define path to modules
% ------------------------------------------------------------------------- 
this_path = mfilename('fullpath');% extract path to this m-file
this_path_parts = regexp(this_path,filesep,'split');% split the path to this m-file into parts.
study_drive = char(this_path_parts(1));

opt_library_path = [study_drive '\_cpeu\optsyssim'];
addpath(genpath(opt_library_path));

% ------------------------------------------------------------------------- 
% Dock figures
% ------------------------------------------------------------------------- 
% set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultFigureWindowStyle','normal');



% -------------------------------------------------------------------------
% TEST 1: Welch Estimate Using Default Inputs
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
n = 0:319;
x = cos(pi/4*n)+randn(size(n));

% Normalised angular frequency is pi/4 rad/sample = 0.7854 rad/sample
% Peak amplitude of the signal is 1 U
% Power of the signal is 1/2 U^2.
% In double-sided, the power of each peak should be -6 dBU
% In single-sided, the power of the signal peak should be -3 dBU
% Noise variance is sig^2 = 1
% Single sided spectral noise density: N0 U^2/(rad/sample)
% Noise power = N0*2*2*pi = 1
% Hence N0 = 1/4/pi U^2/(rad/sample) and 10*log10(N0) = 
% 

%----
% Periodogram parameters
%----
block_size = 71;
overlap_samples = 35;
fft_length = 256;
window_type = 'hamming';

%----
% Normalised single-sided PSD
%----
[psd,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

figure('Name','Example 1 - PSD - single-sided')
plot(norm_freq/pi,10*log10(psd));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-12 8]);

%----
% Normalised double-sided PSD
%----
[psd,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'two_sided','psd');


figure('Name','Example 1 - PSD - double-sided')
plot(norm_freq/pi,10*log10(psd));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-12 8]);

%----
% Single-sided power - normalised frequency axis
%----
[power,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','power');

figure('Name','Example 1 - Power - single-sided')
plot(norm_freq/pi,10*log10(power));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('power (dBU)');
grid on;
%ylim([-12 8]);


%----
% Double-sided power - normalised frequency axis
%----
[power,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'two_sided','power');

figure('Name','Example 1 - Power - double-sided')
plot(norm_freq/pi,10*log10(power));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('power (dBU)');
grid on;
%ylim([-12 8]);

%-----------------
% Test result: OK
%-----------------



%%
% -------------------------------------------------------------------------
% TEST 2: Welch Estimate Using Specified Segment Length
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
n = 0:511;
x = cos(pi/3*n)+randn(size(n));

%----
% Periodogram parameters
%----
block_size = 132;
overlap_samples = 66;
fft_length = 256;
window_type = 'hamming';


[psd,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

figure('Name','Test 2 - Welch Estimate Using Specified Segment Length')
plot(norm_freq/pi,10*log10(psd));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-15 15]);

%-----------------
% Test result: OK
%-----------------



%%
% -------------------------------------------------------------------------
% TEST 3: Welch Estimate Specifying Segment Overlap
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
n = 0:319;
x = cos(pi/4*n)+randn(size(n));


%----
% Periodogram parameters
%----
block_size = 100;
overlap_samples = 25;
fft_length = 256;
window_type = 'hamming';


[psd,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

figure('Name','Test 3 - Welch Estimate Specifying Segment Overlap')
plot(norm_freq/pi,10*log10(psd));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-15 10]);

%-----------------
% Test result: OK
%-----------------


%%
% -------------------------------------------------------------------------
% TEST 4: Welch Estimate Using Specified DFT Length
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
n = 0:319;
x = cos(pi/4*n)+randn(size(n));


%----
% Periodogram parameters
%----
block_size = 100;
overlap_samples = 50;
fft_length = 640;
window_type = 'hamming';


[psd,norm_freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

figure('Name','Test 4 - Welch Estimate Using Specified DFT Length')
plot(norm_freq/pi,10*log10(psd));
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-12 8]);

%-----------------
% Test result: OK
%-----------------


%%
% -------------------------------------------------------------------------
% TEST 5: Welch PSD Estimate of Signal with Frequency in Hertz
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
fs = 1000;
t = 0:1/fs:5-1/fs;
x = cos(2*pi*100*t) + randn(size(t));

%----
% Periodogram parameters
%----
block_size = 500;
overlap_samples = 300;
fft_length = 500;
window_type = 'hamming';




[psd,Freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','psd',fs);

figure('Name','Test 5 - Welch PSD Estimate of Signal with Frequency in Hertz - PSD')
plot(Freq,10*log10(psd));
xlabel('frequency (Hz)');
ylabel('PSD (dBU -10 log_{10}(Hz))');
grid on;
ylim([-35 -5]);

[power,Freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'one_sided','power',fs);

figure('Name','Test 5 - Welch PSD Estimate of Signal with Frequency in Hertz - power')
plot(Freq,10*log10(power));
xlabel('frequency (Hz)');
ylabel('power (dBU)');
grid on;
ylim([-35 -5]);

%-----------------
% Test result: OK
%-----------------



%%
% -------------------------------------------------------------------------
% TEST 6: DC-Centered Power Spectrum
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
fs = 1000;
t = 0:1/fs:5-1/fs;
noisevar = 1/4;
x = cos(2*pi*100*t)+sqrt(noisevar)*randn(size(t));

%----
% Periodogram parameters
%----
block_size = 500;
overlap_samples = 300;
fft_length = 500;
window_type = 'hamming';


[power,Freq] = psd_welch(x,block_size,overlap_samples,window_type,fft_length,'two_sided','power',fs);

figure('Name','Test 6 - DC-Centered Power Spectrum')
plot(Freq,10*log10(power));
xlabel('frequency (Hz)');
ylabel('power (dBU)');
grid on;
ylim([-35 -5]);

%-----------------
% Test result: OK
%-----------------




%%
% -------------------------------------------------------------------------
% TEST 7: Welch PSD Estimate of a Multichannel Signal
% -------------------------------------------------------------------------

%----
% Reset random number generator
%----
rng default

%----
% Generate signal
%----
N = 1024;
n = 0:N-1;

w = pi./[2;3;4];
x = cos(w*n)' + randn(length(n),3);

%----
% Periodogram parameters
%----
block_size = N/4;
overlap_samples = 0;
fft_length = max(256,2^nextpow2(length(block_size)));
window_type = 'hamming';
% These correspond to the default parametes of Matlab pwelch(x)

% [psd,normFreq] = psd_welch(x);
[psd1,normFreq1] = psd_welch(x(:,1),block_size,overlap_samples,window_type,fft_length,'one_sided','psd');
[psd2,normFreq2] = psd_welch(x(:,2),block_size,overlap_samples,window_type,fft_length,'one_sided','psd');
[psd3,normFreq3] = psd_welch(x(:,3),block_size,overlap_samples,window_type,fft_length,'one_sided','psd');
% Our implementation of psd_welch does not operate on "multichannel" inputs.
% We need to call it individually for each signal (column of the input
% matrix x).

figure('Name','Test 7 - Welch PSD Estimate of a Multichannel Signal')
plot(normFreq1/pi,10*log10(psd1),'Color',[0    0.4470    0.7410]);
hold on;
plot(normFreq2/pi,10*log10(psd2),'Color', [0.8500    0.3250    0.0980]);
plot(normFreq3/pi,10*log10(psd3),'Color',[0.9290    0.6940    0.1250]);
xlabel('normalised frequency (\times \pi rad/sample)');
ylabel('PSD (dBU -10 log_{10}(rad/sample))');
grid on;
ylim([-15 15]);

%-----------------
% Test result: Failed. But we are unsure of the parameters actually used in
% the Matlab documentation exemple. Changing the number of segments from
% N/8 to N/4 gives closer results.
%-----------------















