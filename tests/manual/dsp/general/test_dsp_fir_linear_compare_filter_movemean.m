% -------------------------------------------------------------------------
% Comparison of moving average performed by the following methods:
% 1. Matlab movmean function
% 2. Matlab filter function
% 3. Our own implementation of a FIR filter: dsp_fir_linear.m
%
% The 3 methods are applied to a sinusoidal noisy signal.
%
% 2021-03-09
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all;
close all;

% -------------------------------------------------------------------------
% Reset random numbers generator
% -------------------------------------------------------------------------
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
nsamples = 2^10;    % Number of samples. 
noise_var = 0.3;    % Noise variance.
block_length = 100;

% -------------------------------------------------------------------------
% Define signal
% -------------------------------------------------------------------------
nperiod = 512;
x = cos(2*pi*[0:1:nsamples - 1]/nperiod);

figure('Name','Noise-free signal')
plot(x)
xlabel('sample')
ylabel('x')
xlim([1 nsamples])


% -------------------------------------------------------------------------
% Add AWGN
% -------------------------------------------------------------------------
xn = x + sqrt(noise_var)*randn(1,nsamples);
% Add white Gaussian noise.   

figure('Name','Signal with AWGN')
plot(xn)
xlabel('sample')
ylabel('x')
xlim([1 nsamples])


% -------------------------------------------------------------------------
% Filtering
% -------------------------------------------------------------------------

% ----
% movmean
% ----

ym1 = movmean(x,[block_length - 1 0],'EndPoints','discard');
% Discards the elements for which a full averaging over block_length cannot
% be done. 
% The length of the vector should be nsamples - block_length + 1.

ym1_base = [NaN*ones(1,block_length -1) ym1];
% To keep the vector consistent with the indices of the samples, we add
% block_length - 1 NaN in front.
% The length of the vector should now be nsamples.

ym2 = movmean(x,[block_length - 1 0],'EndPoints','fill');
% The matlab movmean function actually has an option to substitute
% nonexisting elements with NaN.
% We check that ym1_base & ym2 are strictly identical.

fprintf('\n\n%s\n','Matlab built in movmean function:')
tic
ym3 = movmean(x,[block_length - 1 0]);
% In this case, the first block_length - 1 elements of the output are 
% calculated based on the summation of less than block_length input
% elements.
% We check that ym1_base, ym2, and ym3 are identical for indices greater or
% equal to block_length.
toc

% ----
% filter
% ----

fprintf('\n%s\n','Matlab built in filter function:')
tic
a = 1;
b = ones(1,block_length);
yf = filter(b,a,x)/block_length;
toc

% ----
% dsp_fir_linear
% ----

fprintf('\n%s\n','dsp_fir_linear function:')
tic
z = zeros(1,block_length);
b = ones(1,block_length);
yd = dsp_fir_linear(x,b,z)/block_length;
toc

% -------------------------------------------------------------------------
% Compare filtered waveforms
% -------------------------------------------------------------------------
figure('Name','Compare filtered waveforms')
plot(x,'k-')
hold on
plot(ym1_base,'b--')
plot(yf,'g--' )
plot(yd,'r:');
xlabel('sample')
ylabel('y')
xlim([1 nsamples])
legend('original (noise-free)','movmean','filter','dsp-fir-linear')


% -------------------------------------------------------------------------
% Compensate delay
% -------------------------------------------------------------------------
ym1_base_comp = circshift(ym1_base,[0 -block_length/2]);

gate = zeros(1,nsamples);
gate(300:700) = 1;
[~,imax_x]  = max(x.*gate)
[~,imax_y]  = max(ym1_base_comp.*gate)
% We check that the maximum of the sinusoidal signals overlap.

figure('Name','After delay compensation')
plot(x,'k--')
hold on
plot(ym1_base_comp,'b-')
xlabel('sample')
ylabel('y')
xlim([1 nsamples])
legend('original noise-free sinusoidal signal','after movmean + delay compensation')