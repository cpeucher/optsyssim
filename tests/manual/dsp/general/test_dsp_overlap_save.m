% -------------------------------------------------------------------------
% Test of digital filtering by overlap-save method. Function:
% dsp_overlap_save.
%
% We generate a PAM2 signal using digital pulse shaping with a
% raised-cosine impulse.
% The filtering is performed using either a standard FIR filter in the time
% domain, or in the frequency domain using the overlap-save method.
%
% 2021-10-08
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all
format long

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
fig_interpreter = 'latex';




% -------------------------------------------------------------------------
% Global parameters for digital impulse synthesis.
% -------------------------------------------------------------------------
nsymbols = 512*16;
symbol_rate = 32e9;

ts = 1/symbol_rate;
% Symbol duration, in s.

% -------------------------------------------------------------------------
% Global parameters for analog waveform synthesis.
% -------------------------------------------------------------------------
global reference_frequency frequency_array time_array dt df
reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);




% -------------------------------------------------------------------------
% Determine tap coefficients h[n] of the pulse-shaping filter
% -------------------------------------------------------------------------
% This is done by time domain sampling of the raised cosine impulse
% response.

nlengthsymbs = 16;
% Filter depth, in terms of number of symbols.
% 16 means 8 symbols on each side.
% Even number.
nos = 4;
% Oversampling factor. Even number.

ntaps = nlengthsymbs*nos + 1;
% FIR filter length. Ensure it is an odd number.
nn = (ntaps - 1)/2;
% nn so that ntaps = 2*nn + 1
% Integer since ntaps is odd.

tsamp = ((0:ntaps-1) - nn)*ts/nos;
% Sampling time instants.


roll_off = 0.1;
ht = elec_pulse_rc(tsamp,0,ts,roll_off);
% Sample the filter analogue impulse response.

figure('Name','FIR impulse response')
stem(tsamp*symbol_rate,ht,'r')
xlim([-1,1]*nn/nos)
xlabel('time t/T_s')
ylabel('amplitude')


% -------------------------------------------------------------------------
% Generate modulated PAM2 signal with RC pulse shaping.
% -------------------------------------------------------------------------
m = 2;
bit_pattern = dsp_data_binary(nsymbols,m);
% Generate bit pattern.

sig_imp = dsp_upsample(bit_pattern,nos);
% Impulse signal, up-sampled to match the FIR filter.


tic
[sig_fir,~] = dsp_fir_linear(sig_imp,ht,zeros(1,ntaps));
% Apply the filter to the impulse.
toc

tic
sig_os = dsp_overlap_save(sig_imp,ht,ntaps - 1);
toc

figure('Name','PAM2 - digital filter output')
stem(sig_fir,'bo-')
hold on
stem(sig_os,'r*--')
legend('FIR','overlap-save')
xlim([0 500])
xlabel('sample')
ylabel('amplitude')

% -------------------------------------------------------------------------
% Digital-to-analog conversion. To check the generated signal makes
% sense...
% -------------------------------------------------------------------------
sig_zoh = dac_zoh(sig_fir,nsamples/length(sig_fir));
% Basic zero-order hold.

params_elpf.type = 'gaussian';
params_elpf.order = 4;
params_elpf.f3dB = 0.625*symbol_rate;
sig_lpf = elec_elpf(sig_zoh,params_elpf);
% Electrical low-pass filter.

params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'Synthesised PAM2 waveform';
meas_eye(sig_lpf,params_eye);
% Eye diagram.



% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

