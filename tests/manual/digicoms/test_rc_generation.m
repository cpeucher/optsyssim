% -------------------------------------------------------------------------
% Test of digital synthesis of raised-cosine pulses.
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
% Global parameteres
% -------------------------------------------------------------------------
global reference_frequency frequency_array time_array dt df

reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
nsymbols = 128*4;
symbol_rate = 32e9;

% -----------------------------------------------------------------
% Derived parameters
% -----------------------------------------------------------------
nsamples = nsamples_per_symbol*nsymbols;
% Total number of samples.
sample_rate = nsamples_per_symbol*symbol_rate;
% Simulation sample rate.

% -----------------------------------------------------------------
% Create time and frequency axes
% -----------------------------------------------------------------
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);



% -------------------------------------------------------------------------
% Analogue impulse response / spectrum
% -------------------------------------------------------------------------
roll_off = 0.1;
ts = 1/symbol_rate;

pulse = elec_pulse_rc(time_array,time_array(length(time_array)/2),ts,roll_off);

figure('Name','analogue pulse shape')
plot((time_array - time_array(length(time_array)/2))*symbol_rate,pulse,'b')
xlabel('normalised time t/T_s');
ylabel('amplitude')
xlim([-8 8]);


spectrum_ft = num_ft(pulse,dt);
% Calculation of pulse spectrum by Fourier-transforming the pulse shape.

spectrum_direct = calc_rc_spectrum(frequency_array,ts,roll_off);
% Direct calculation of the raised cosine spectrum.


figure('Name','analogue spectrum - log')
plot(frequency_array/1.0e9,20*log10(abs(spectrum_ft)),'b')
hold on
plot(frequency_array/1.0e9,20*log10(abs(spectrum_direct)+eps),'r--')
xlim([-50 50]);
xlabel('frequency (GHz)')
ylabel('|H|^2 (dB)')
legend('Fourier transform','direct calculation')
% We plot the spectra calculated by Fourier transform and through the
% analytical expression.
% We can see the effect of truncation of the impulse reponse in the Fourier
% transform spectrum.
% Increase the duration of the time window to improve the rejection.
% Furtermore, in case roll_off = 0 (sinc pulse) we observe Gibbs
% oscillations. 
% This can be cured by introducing a proper windowing function, e.g.
% params_window.length = nsamples;
% params_window.type = 'hamming';
% params_window.symmetry = 'symmetric';
% ww = dsp_window(params_window);
% spectrum_ft = num_ft(pulse.*ww);


figure('Name','analogue spectrum - lin')
plot(frequency_array/1.0e9,abs(spectrum_ft),'b')
hold on
plot(frequency_array/1.0e9,spectrum_direct,'r--')
xlim([-50 50]);
legend('Fourier transform','direct calculation')
xlabel('frequency (GHz)')
ylabel('|H|')


%%
% -------------------------------------------------------------------------
% FIR shaping filter synthesis - time sampling
% -------------------------------------------------------------------------
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

ht = elec_pulse_rc(tsamp,0,ts,roll_off);
% Sample the filter analogue impulse response.

figure('Name','FIR design - time domain')
plot((time_array - time_array(length(time_array)/2))*symbol_rate,pulse,'b')
hold on
stem(tsamp*symbol_rate,ht,'r')
xlim([-1,1]*nn/nos)
xlabel('time t/T_s')
ylabel('amplitude')
legend('analogue pulse shape','FIR filter coefficient')


% -------------------------------------------------------------------------
% FIR shaping filter synthesis - frequency sampling
% -------------------------------------------------------------------------
fsa = nos*symbol_rate;
% Sampling frequency.


freq_sampling = linspace(0,fsa/2,1000);
spectrum_sampling = calc_rc_spectrum(freq_sampling,ts,roll_off);

spectrum_sampling = spectrum_sampling/spectrum_sampling(1);



hf = dsp_fir_design_frequency_sampling(freq_sampling,spectrum_sampling,fsa,ntaps,1);

figure('Name','Compare time- and frequency-domain FIR coefficients')
stem([-nn:nn],ht/max(ht),'bs')
hold on
stem([-nn:nn],hf/max(hf),'r*')
xlabel('tap index')
ylabel('amplitude (a.u.)')
legend('time sampling','frequency sampling')


% -------------------------------------------------------------------------
% Generation of rc pulse from digital samples
% -------------------------------------------------------------------------
ximpulse = zeros(1,nsymbols*nos);
ximpulse(length(ximpulse)/2) = 1;
% Create an impulse.

[yimpulse,~] = dsp_fir_linear(ximpulse,ht,zeros(1,ntaps));
% Apply the filter to the impulse.

yimpulse = circshift(yimpulse,[1,-nn+1]);
% Compensate the delay of the FIR filter.



figure('Name','FIR pulse synthesis - Digital filter output')
stem([0:length(yimpulse)-1]-length(yimpulse)/2,yimpulse,'bo')
xlabel('sample')
ylabel('amplitude (a.u.)')
xlim([-100 100])

sig_imp = dac_zoh(yimpulse,nsamples/length(yimpulse));
% Basic zero-order hold.

figure('Name','FIR pulse synthesis - After ZOH')
plot(time_array - time_array(nsamples/2),sig_imp,'r')
hold on
plot(time_array - time_array(nsamples/2),pulse,'b')
xlabel('time (s)')
ylabel('amplitude (a.u.)')
xlim([-4 4]*1.0e-10)
legend('ZOH outout','analog waveform')


%%
% -------------------------------------------------------------------------
% Generate PAM2 signal.
% -------------------------------------------------------------------------
% We generate a raised-cosine PAM2 signal using digital filtering of an
% impulse signal, followed by zero-order hold and low-pass filtering.

m = 2;
bit_pattern = generate_binary(nsymbols,m);
% Generate bit pattern.

sig_imp = dsp_upsample(bit_pattern,nos);
% Impulse signal, up-sampled to match the FIR filter.

figure('Name','PAM2 RC generation - impulse signal')
stem(sig_imp,'b-')
xlabel('sample')
ylabel('amplitude (a.u.)')


[sig_fir,~] = dsp_fir_linear(sig_imp,ht,zeros(1,ntaps));
% Apply the filter to the impulse.

sig_os = dsp_overlap_save(sig_imp,ht,ntaps - 1);
% We also filter using the overlap-save method.


figure('Name','PAM2 - digital filter output - direct FIR & OS')
stem(sig_fir,'bo-')
hold on
stem(sig_os,'r*--')
legend('FIR','overlap-save')
xlim([0 500])
xlabel('sample')
ylabel('amplitude (a.u.)')

% We observe the outputs are identical in the direct FIR calculation and
% the calculation by overlap-save method.



sig_zoh = dac_zoh(sig_fir,nsamples/length(sig_fir));
% Basic zero-order hold.


params_elpf.type = 'gaussian';
params_elpf.order = 4;
params_elpf.f3dB = 0.625*symbol_rate;
sig_lpf = elec_elpf(sig_zoh,params_elpf);
% Electrical low-pass filter.

figure('Name','PAM2 RC generation - after FIR + ZOH + LPF')
plot(time_array/1.0e-12,sig_lpf,'b-')
xlabel('time (ps)')
ylabel('amplitude (a.u.)')


params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'PAM2 RC generation - after FIR + ZOH + ELPF';
meas_eye(sig_lpf,params_eye);
% Eye diagram.

% We compare with an ideal RC PAM2 signal generation obtained from the
% convolution of the impulse signal with a theoretical RC pulse shape.

params_mod.pulse_shape = 'rc';
params_mod.roll_off = roll_off;
params_mod.symbol_rate = symbol_rate;
params_mod.normalisation = 'peak';
sig_analog = elec_modulator(bit_pattern,params_mod);
params_eye.name = 'PAM2 RC generation - with ideal analog pulse shape';
meas_eye(sig_analog,params_eye);




%%
% -------------------------------------------------------------------------
% Generate QPSK signal
% -------------------------------------------------------------------------
m = 4;
% Constellation order. M points with log2(M) even (square constellation)
data = generate_binary(nsymbols,m);
% Generate binary data.

[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words.

symbs = diffenc_qpsk(words_dec);
% Differential encoding.


sig_imp = dsp_upsample(symbs,nos);
% Impulse signal.

[sig_fir,~] = dsp_fir_linear(sig_imp,ht,zeros(1,ntaps));
% Apply the filter to the impulse signal.

sig_zoh = dac_zoh(sig_fir,nsamples/length(sig_fir));
% Basic zero-order hold.

sig_lpf = elec_elpf(sig_zoh,params_elpf);
% Electrical low-pass filter.

params_eye.name = 'QPSK RC generation - after FIR + ZOH + ELPF';
meas_eye(sig_lpf(8*nsamples_per_symbol+1:end),params_eye);
% We remove first symbols since they correspond to zero power due to
% initialisation of the FIR. Not good for QPSK.

% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(5)

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

