% -------------------------------------------------------------------------
% Test of equivalent noise bandwidth calculation
%
% We generate real and complex white Gaussian noise, filter it through an
% electrical low-pass filter and an optical bandpass filter, respectively,
% calulate the equivalent noise bandwidth of the filter, and check that the
% input and output noise power levels are as expected.
%
% 2022-04-28
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Preparation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Clean up
% ------------------------------------------------------------------------- 
clear all
close all

% ------------------------------------------------------------------------- 
% Specify display format
% ------------------------------------------------------------------------- 
format long

% ------------------------------------------------------------------------- 
% Dock figures
% ------------------------------------------------------------------------- 
% set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultFigureWindowStyle','normal');

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility of the
% results.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
% reset(stream);

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Global simulation parameters
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Define global variables
% ------------------------------------------------------------------------- 
global reference_frequency 
global frequency_array 
global time_array
global dt
global df 
global CONSTANT
% global space_grid

% -------------------------------------------------------------------------
% Set global simulation parameters
% -------------------------------------------------------------------------
% reference_frequency = 193.1e12;
% nsamples_per_symbol = 128;
% nsymbols = 128;
% symbol_rate = 40e9;
% 
% nsamples = nsamples_per_symbol*nsymbols;
% sample_rate = nsamples_per_symbol*symbol_rate;
% [time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);        
        
  
reference_frequency = 193.1e12;
df = 10e6;
nsamples = 2^20;
 
sample_rate = nsamples*df;
dt = 1/sample_rate;
time_array = (0:nsamples-1)*dt;
frequency_array = (-nsamples/2:nsamples/2-1)*df;

% -------------------------------------------------------------------------
% Space grid
% -------------------------------------------------------------------------
% xrange = [-15e-6, 15e-6];  
% yrange = [-15e-6, 15e-6]; 
% nxpoints = 2001;
% nypoints = 2001;
% space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
% [space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Startup routines
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% ------------------------------------------------------------------------- 
% Load physical constants
% ------------------------------------------------------------------------- 
CONSTANT = core_load_constants();
% Load essential physical constants.

% ------------------------------------------------------------------------- 
% Start time 
% ------------------------------------------------------------------------- 
start_time = datetime("now");
fprintf('\n\n%s%s\n\n','Simulation started on ',start_time);


%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
fig_interpreter = 'latex';



% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% Now we are ready to implement the system
% -------------------------------------------------------------------------        
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 


% -------------------------------------------------------------------------
% First, real noise
% -------------------------------------------------------------------------


npsd2 = 5;
% Two-sided noise spectral density, in W/Hz.
% Obs: this is the 2-sided noise spectral density.
% One needs to integrate over both negative and positive frequencies to
% obtain the total noise power.
% According to usual conventions, if N0 is the 1-sided power spectral
% density, then npsd2 = N0/2.

noise_power_specified = npsd2*sample_rate;
% Total noise power.


sig = sqrt(noise_power_specified)*randn(1,nsamples);
% Generate white Gaussian noise.

noise_power_in = var(sig);
% Estimate generated noise power.

fprintf('\n\n\n%s\n','Test 1: real noise and ELPF')
fprintf('%s','===========================')

fprintf('\n\n%s\t%3.2f%s\n','Specified noise power:',noise_power_specified,' W');
fprintf('%s\t\t%3.2f%s\n','Actual noise power:',noise_power_in,' W');

fprintf('\n%s\t%3.8f%s\n','Specified 2-sided PSD:',npsd2,' W/Hz');
fprintf('%s\t\t%3.8f%s\n','Actual 2-sided PSD:',noise_power_in/sample_rate,' W/Hz');


% We estimate the noise PSD using Welch periodogram
block_size = 100;
overlap_samples = 50;
fft_length = 128;
window_type = 'hann';
[psd,freq] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'two_sided','psd',sample_rate);
% Calculate power spectral density using Welch periodogram method.

figure('Name','power spectral density')
plot(freq,psd)
hline(npsd2,'r--')
ylim([0 npsd2*1.3])
xlabel('frequency (Hz)')
ylabel('PSD (W/Hz)')

params_elpf.type = 'gaussian';
params_elpf.order = 1;
params_elpf.f3dB = 200e9;
tf = elec_tf_elpf(params_elpf,frequency_array);
% Generate filter transfer function

figure('Name','electrical filter transfer function')
plot(frequency_array,abs(tf).^2)
xlabel('frequency (Hz)')
ylabel('transmission')


tf0 = 1;
enb = calc_enb('lowpass',tf,tf0,df);
% Calculate the equivalent noise bandwidth of the filter.
% We consider the filter to be low-pass.

factor = enb/params_elpf.f3dB;

fprintf('\n%s\t%3.8f%s\n','Filter cut-off frequency:',params_elpf.f3dB,' Hz');
fprintf('%s\t\t%3.8f%s\n','Filter ENB (low-pass):',enb,' Hz');
fprintf('%s\t\t%3.8f\n','Filter ENB/f3dB factor (low-pass):',factor);


sig_out = elec_filter(sig,tf);
% Filter the noise.

noise_power_out = var(sig_out);
% Estimate the noise power at output of filter.

noise_power_out_expected = npsd2*num_int1d_simpson(abs(tf).^2,df);
% The expected noise power is the input noise psd multiplied by the
% integral of the squared modulus of the filter transfer function.
% This is the total power over positive and negative frequencies.

noise_power_out_expected_from_enb = 2*enb*npsd2;
% Expected output noise power from the enb. The factor 2 stems from the
% fact that the noise is integrated over all positive and negative
% frequencies while the low-pass enb is calculated.
% No surprise, this is the same calculation as the one resulting in
% noise_power_out_expected...



fprintf('\n%s\t\t\t%3.8f%s\n','Expected output noise power:',noise_power_out_expected,' W');
fprintf('%s\t%3.8f%s\n','Expected output noise power (from ENB):',noise_power_out_expected_from_enb,' W');
fprintf('%s\t\t\t\t%3.8f%s\n','Actual output noise power:',noise_power_out,' W');



%%
% -------------------------------------------------------------------------
% Second, complex noise
% -------------------------------------------------------------------------


npsd2 = 5;
% Two-sided noise spectral density, in W/Hz.
% Obs: this is the 2-sided noise spectral density.
% One needs to integrate over both negative and positive frequencies to
% obtain the total noise power.
% This is the total noise power (for both quadratures)

noise_power_specified = npsd2*sample_rate;
% Total noise power.

sig = sqrt(noise_power_specified/2)*randn(1,nsamples) + 1j*sqrt(noise_power_specified/2)*randn(1,nsamples);
% Generate white Gaussian noise.
% We split the the noise power onto both quadrature.

noise_power_in_i = var(real(sig));
% Actual noise power: in-phase.
noise_power_in_q = var(imag(sig));
% Actual noise power: quadrature.
noise_power_in_sum = noise_power_in_i + noise_power_in_q;
% Actual noise power: sum of noise power in both quadratures.

noise_power_in_tot = var(sig);
% We check that var applied to a complex variable returns the sum of the
% variances of the real and imaginary parts.




fprintf('\n\n\n%s\n','Test 2: complex noise and OBPF')
fprintf('%s','==============================')

fprintf('\n\n%s\t\t\t\t%3.2f%s\n','Specified noise power:',noise_power_specified,' W');
fprintf('%s\t\t%3.2f%s\n','Actual noise power (in-phase):',noise_power_in_i,' W');
fprintf('%s\t%3.2f%s\n','Actual noise power (quadrature):',noise_power_in_q,' W');
fprintf('%s\t\t\t%3.2f%s\n','Actual noise power (sum):',noise_power_in_sum,' W');
fprintf('%s\t\t\t%3.2f%s\n','Actual noise power (total):',noise_power_in_tot,' W');

fprintf('\n%s\t%3.8f%s\n','Specified 2-sided PSD:',npsd2,' W/Hz');
fprintf('%s\t\t%3.8f%s\n','Actual 2-sided PSD:',noise_power_in_sum/sample_rate,' W/Hz');


% We estimate the noise PSD using Welch periodogram
block_size = 100;
overlap_samples = 50;
fft_length = 128;
window_type = 'hann';
[psd,freq] = psd_welch(sig,block_size,overlap_samples,window_type,fft_length,'two_sided','psd',sample_rate);
% Calculate power spectral density using Welch periodogram method.

figure('Name','power spectral density')
plot(freq,psd)
hline(npsd2,'r--')
ylim([0 npsd2*1.3])
xlabel('frequency (Hz)')
ylabel('PSD (W/Hz)')


params_obpf.type = 'gaussian';
params_obpf.order = 1;
params_obpf.bandwidth = 400e9;
params_obpf.centre_frequency = 0;
tf = opt_tf_obpf(params_obpf,frequency_array);
% Generate filter transfer function

figure('Name','optical filter transfer function')
plot(frequency_array,abs(tf).^2)
xlabel('frequency (Hz)')
ylabel('transmission')


tf0 = 1;
enb = calc_enb('bandpass',tf,tf0,df);
% Calculate the equivalent noise bandwidth of the filter.
% We consider the filter to be low-pass.

factor = enb/params_obpf.bandwidth;

fprintf('\n%s\t%3.8f%s\n','Filter cut-off frequency:',params_obpf.bandwidth,' Hz');
fprintf('%s\t\t%3.8f%s\n','Filter ENB (low-pass):',enb,' Hz');
fprintf('%s\t\t%3.8f\n','Filter ENB/f3dB factor (low-pass):',factor);


sigopt_in.x = sig;
sigopt_in.y = zeros(1,nsamples);

sigopt_out = opt_filter(sigopt_in,tf);
sig_out = sigopt_out.x;
% Filter the noise.

noise_power_out = var(sig_out);
% Estimate the noise power at output of filter.

noise_power_out_expected = npsd2*num_int1d_simpson(abs(tf).^2,df);
% The expected noise power is the input noise psd multiplied by the
% integral of the squared modulus of the filter transfer function.
% This is the total power over positive and negative frequencies.

noise_power_out_expected_from_enb = enb*npsd2;
% Expected output noise power from the enb. 
% No surprise, this is the same calculation as the one resulting in
% noise_power_out_expected...



fprintf('\n%s\t\t\t%3.8f%s\n','Expected output noise power:',noise_power_out_expected,' W');
fprintf('%s\t%3.8f%s\n','Expected output noise power (from ENB):',noise_power_out_expected_from_enb,' W');
fprintf('%s\t\t\t\t%3.8f%s\n','Actual output noise power:',noise_power_out,' W');








        







%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');        
        

% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% Display duration (or wrapup: core_wrapup(start_time,0);)
% ------------------------------------------------------------------------- 
core_display_duration(start_time,datetime("now"));