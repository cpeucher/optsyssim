




% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% amplifiers                                                                amplifiers
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% opt_amplifier
% Basic system model of an optical amplifier
% /src/amplifiers/
% -------------------------------------------------------------------------
params_optamp.mode = 'power';%'pain';%'saturation';
params_optamp.pol = 'x';%'y';%'both';
params_optamp.gain = 20;
params_optamp.output_power = 10;
params_optamp.noise_figure = 4;
params_optamp.add_noise = 1;
sig = opt_amplifier(sig,params_optamp);
% Optical amplifier


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% auxiliary                                                                 auxiliary
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% calc_disp_broadening_gauss_1
% 2nd and 3rd order dispersion-induced rms broadening of 1st order Gaussian pulse
% /src/auxiliary/
% -------------------------------------------------------------------------
params.fwhm = 10e-12;
params.chirp = 0;
[broadening_rms,ld,ldp] = calc_disp_broadening_gauss_1(distance,beta2,beta3,params);
% Pulse rms broadening factor

% -------------------------------------------------------------------------
% calc_disp_broadening_gauss_m
% Second-order dispersion-induced rms broadening of super-Gaussian pulses
% /src/auxiliary/
% -------------------------------------------------------------------------
params.fwhm = 10e-12;
params.chirp = 0;
params.order = 2;
[broadening_rms,ld] = calc_disp_broadening_gauss_m(distance,beta2,params);
% Pulse rms broadening factor

% -------------------------------------------------------------------------
% calc_disp_imdd_small_signal
% Fibre frequency response due to dispersion and chirp with direct detection
% /src/auxiliary/
% -------------------------------------------------------------------------
freq = [0:0.1:20]*1.0e9;                   % frequency range, in Hz
lambda = CONSTANT.c/reference_frequency;   % laser wavelength, in m.
dispersion = 17*80*1e-3;                   % dispersion, in s/m.
params_chirpedtx.alpha = 0;                % laser alpha parameter.
params_chirpedtx.fc = 0;                   % chirp corner frequency, in Hz.
hresp = calc_disp_imdd_small_signal(freq,lambda,dispersion,params_chirpedtx); 
% Fibre + DD small-signal frequency response

% -------------------------------------------------------------------------
% calc_dml_frequency_response
% Calculation of small signal frequency response of a DML
% /src/auxiliary/
% -------------------------------------------------------------------------
params_dml.tau_p = 2.6e-12;          % photon lifetime, in s
params_dml.tau_c = 3.17e-9;          % carrier lifetime, in s
params_dml.n_0 = 2.0e24;             % carrier density at transparency, in 1/m^3
params_dml.sigma_g = 3.34e-20;       % gain cross section, in m^2
params_dml.n_g = 4;                  % group effective index
params_dml.Gamma = 0.2408;           % confinement factor
params_dml.V = 3.6e-17;              % active volume, in m^3
params_dml.epsilon_nl = 2.0e-23;     % gain suppression factor
params_dml.alpha = 6;                % linewidth enhancement factor
params_dml.beta = 1.0e-3;            % spontaneous emission factor
params_dml.eta_0 = 0.2;              % differential quantum efficiency
params_dml.emission_frequency = reference_frequency; % emission frequency, in Hz
ibias = 20e-3;                       % bias current, in A
freq = [0:0.1:20]*1e9;               % frequency range, in Hz
dml_response = calc_dml_frequency_response(ibias,params_dml,freq);
% DML small-signal frequency response

% -------------------------------------------------------------------------
% conv_disp_d_beta
% Conversion of dispersion between D, S, C and beta_n
% /src/auxiliary/
% -------------------------------------------------------------------------
dispersion = 16;              % dispersion, in ps/nm/km
dispersion_slope = 0.058;     % dispersion slope, in ps/nm2/km
dispersion_curvature = 0;     % dispersion curvature, in ps/nm3/km
dispersion_spec_frequency = CONSTANT.c/1550e-9;
beta = conv_disp_d_beta([dispersion dispersion_slope dispersion_curvature],'to_beta','eng','si',dispersion_spec_frequency);
% Convert dispersion to beta

% -------------------------------------------------------------------------
% conv_loss_lin_log
% Conversion of loss per unit length from dB/km to m^-1 and calculation of effective length
% /src/auxiliary/
% -------------------------------------------------------------------------
loss = 0.8;       % Loss in dB/km
length = 50;      % Waveguide length, in m.
[alpha,leff] = conv_loss_lin_log(loss,length);
% Calculation of loss in 1/m and effective length


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% characterization                                                          characterization
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% char_opt_average_power 
% Optical average power meter (PWM)
% /src/characterization/
% -------------------------------------------------------------------------
average_power = char_opt_average_power(sig);
% Calculate signal average power

% -------------------------------------------------------------------------
% char_opt_constellation
% Constellation diagram for optical signals
% /src/characterization/
% -------------------------------------------------------------------------
params_constellation.type = 'scatter';%'transitions';%'scatter_transitions';
params_constellation.pol = 'x';% 'y','both';
params_constellation.nsamples_per_symbol = nsamples_per_symbol;
params_constellation.sampling_time = nsamples_per_symbol/2;
params_constellation.normalisation = 'mean';%'max';
params_constellation.save.emf = 0;
% params_constellation.line_color = 'b';
% params_constellation.line_width = 3;
% params_constellation.marker_type = 'o';
% params_constellation.marker_color = 'r';
% params_constellation.marker_size = 36;
% params_constellation.plot_circular_grid = 1;
% params_constellation.plot_axes = 0;
params_constellation.name = 'Constellation';
char_opt_constellation(sig,params_constellation); 
% Plot optical constellation

% -------------------------------------------------------------------------
% char_opt_energy
% Energy of an optical signal in a specified time interval
% /src/characterization/
% -------------------------------------------------------------------------
ep = char_opt_energy(sig,[time_array(1),time_array(end)]);
% Energy of the signal

% -------------------------------------------------------------------------
% char_opt_peak_power
% Peak power of optical signal in a given time interval
% /src/characterization/
% -------------------------------------------------------------------------
pp = char_opt_peak_power(sig,time_interval);
% Peak power of optical signal

% -------------------------------------------------------------------------
% char_opt_temporal_chirp
% Temporal chirp of an optical signal
% /src/characterization/
% -------------------------------------------------------------------------
chirp = char_opt_temporal_chirp(sig);
% Calculate temporal chirp

% -------------------------------------------------------------------------
% char_pulse_rms
% rms width of an optical pulse
% /src/characterization/
% -------------------------------------------------------------------------
rms = char_pulse_rms(abs(sig.x).^2); 
% Pulse rms duration, in s

% -------------------------------------------------------------------------
% meas_esa
% Electrical spectrum analyser
% /src/characterization/
% -------------------------------------------------------------------------
params_esa.display_interval = [0 frequency_array(end)];
params_esa.resolution_bandwidth = 0;
params_esa.input_impedance = 1;
params_esa.display = 1;
params_esa.save.txt = 0;
params_esa.save.emf = 0;
params_esa.name = 'RF spectrum';
spectrum = meas_esa(sig,params_esa); 
% Electrical spectrum analyser

% -------------------------------------------------------------------------
% mes_eye
% Visualisation of eye diagram of optical or electrical signals
% /src/characterization/
% -------------------------------------------------------------------------
params_eye.pol = 'x';%'y','both';
params_eye.neyes = 2;
params_eye.nsamples_per_symbol = nsamples_per_symbol;
params_eye.save.txt = 0;
params_eye.save.emf = 0;
params_eye.save.jpg = 0;
params_eye.save.vertical_scale_type = 'auto';%'fixed';
params_eye.save.vertical_scale_max = 1.0e-3;
params_eye.save.display_baseline = 1;
params_eye.colour_grade = 0;
params_eye.name = 'Eye diagram';
meas_eye(sig,params_eye);
% Plot eye diagram

% -------------------------------------------------------------------------
% meas_osa
% Optical spectrum analyser
% /src/characterization/
% -------------------------------------------------------------------------
params_osa.pol = 'x';%'y';'both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 0;%12.5e9;
params_osa.sensitivity = -40;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum';
meas_osa(sig,params_osa); 
% Optical spectrum analyser

% -------------------------------------------------------------------------
% meas_rf_power
% Integrate RF power in a specified bandwidth
% /src/characterization/
% -------------------------------------------------------------------------
params_rf_pwm.centre_frequency = 10e9;
params_rf_pwm.bandwidth = 4*df;
params_rf_pwm.input_impedance = 1;
[rf_power,psd] = meas_rf_power(sig,params_rf_pwm);
% Extract RF power in specified bandwidth

% -------------------------------------------------------------------------
% meas_scope
% Oscilloscope for visualisation of optical and electrical signals
% /src/characterization/
% -------------------------------------------------------------------------
params_scope.visualisers = {'amplitude','power','phase','chirp'};
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'Waveform';
meas_scope(sig,params_scope); 
% Scope




% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% core                                                                      core
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% core_create_time_axis
% Create time and frequency axes
% /src/core/
% -------------------------------------------------------------------------
nsamples_per_symbol = 32;
nsymbols = 128;
symbol_rate = 25e9;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);
% Create time and frequency axes

% -------------------------------------------------------------------------
% core_display_duration
% Display time interval between two events
% /src/core/
% -------------------------------------------------------------------------
start_time = datetime("now");
end_time = datetime("now");
core_display_duration(start_time,end_time);
% Display simulation duration

% -------------------------------------------------------------------------
% isoptical
% Test if a signal is optical
% /src/core/
% -------------------------------------------------------------------------
[fflag,message] = isoptical(sig);

% -------------------------------------------------------------------------
% core_load_constants
% Load essential physical constants
% /src/core/
% -------------------------------------------------------------------------
global CONSTANT
CONSTANT = core_load_constants();
% Load essential physical constants



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% digicoms                                                                  digicoms
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% add_awgn
% Addition of complex AWGN to signal samples in the digital domain
% /src/digicoms/
% -------------------------------------------------------------------------
symbs = add_awgn(symbs,esn0_db); 
% Add WGN to symbols

% -------------------------------------------------------------------------
% add_cfo
% Add carrier frequency offset (CFO) to the signal samples
% /src/digicoms/
% -------------------------------------------------------------------------
cfo_normalised = cfo_absolute/symbol_rate;
symbs = add_cfo(symbs,cfo_normalised);
% Add carrier frequency offset

% -------------------------------------------------------------------------
% calc_ber
% Bit-error-ratio calculation by comparison of two binary vectors
% /src/digicoms/
% -------------------------------------------------------------------------
ber = calc_ber(bits,bits_ref);
% BER calculation

% -------------------------------------------------------------------------
% calc_evm
% Error vector magnitude (EVM) calculation
% /src/digicoms/
% -------------------------------------------------------------------------
[evm_max,evm_rms] = calc_evm(constellation,clust);
% EVM calculation

% -------------------------------------------------------------------------
% calc_eye
% Quick eye diagram calculation
% /src/digicoms/
% -------------------------------------------------------------------------
eye = calc_eye(sig,nsamples_per_symbol);
figure('Name','eye diagram')
plot(time_array(1:nsamples_per_symbol),eye,'b');
% Plot eye diagram

% -------------------------------------------------------------------------
% calc_rc_spectrum
% Calculation of raised-cosine pulse spectrum
% /src/digicoms/
% -------------------------------------------------------------------------
roll_off = 0.1;
freq = linspace(0,fs/2,1000);
H = calc_rc_spectrum(freq,symbol_rate,roll_off); 

% -------------------------------------------------------------------------
% calc_ser
% Calculation of symbol error ratio (SER)
% /src/digicoms/
% -------------------------------------------------------------------------
ser = calc_ser(symbs_cx,symbs_ref); 
% SER calculation

% -------------------------------------------------------------------------
% clustering
% Data-aided clustering of received symbols
% /src/digicoms/
% -------------------------------------------------------------------------
symbs_rx = normalise_constellation(symbs_rx,norm_es);
cluster_da = clustering(symbs_rx,words_dec,m);
% Data-aided clustering
symbs_cx = decision_qam_square_hard(symbs_rx,m);
words_dec_rx = demapping(symbs_cx,constellation);
cluster_dd = clustering(symbs_cx,words_dec_rx,m);
% Decision-directed clustering

% -------------------------------------------------------------------------
% conv_bin2dec
% Converts vector of bits (binary) to vector of words (decimal)
% /src/digicoms/
% -------------------------------------------------------------------------
[words_dec,words_bin] = conv_bin2dec(bits_bin,log2(m)); 
% Convert bit string to words

% -------------------------------------------------------------------------
% conv_dec2bin
% Converts vector of words (decimal) to vector of bits (binary)
% /src/digicoms/
% -------------------------------------------------------------------------
bits_bin = conv_dec2bin(words_dec,nob);
% Convert string of words to bits

% -------------------------------------------------------------------------
% decision_med
% Minimum Euclidian distance (MED) decision
% /src/digicoms/
% -------------------------------------------------------------------------
[symbs_cx,words_dec] = decision_med(symbs_rx,constellation);
% Minimum Euclidian distance decision

% -------------------------------------------------------------------------
% decision_qam_square_hard
% Hard decision for square QAM
% /src/digicoms/
% -------------------------------------------------------------------------
symbs_cx = decision_qam_square_hard(symbs_rx,m);
% Hard decision

% -------------------------------------------------------------------------
% decision_qpsk_hard
% Hard decision for QPSK
% /src/digicoms/
% -------------------------------------------------------------------------
es_rx = sum(abs(symbs_rx).^2)/length(symbs_rx);
symbs_rx = symbs_rx*sqrt(2/es_rx);
[symbs_cx,bits,words_dec] = decision_qpsk_hard(symbs_rx,'qpsk_gray');
% Hard decision

% -------------------------------------------------------------------------
% define_constellation
% Define constellations of digital modulation formats
% /src/digicoms/
% -------------------------------------------------------------------------
[constellation,norm_es,norm_emax] = define_constellation(type,m);
% Define constellation

% -------------------------------------------------------------------------
% demapping
% Demapping of digital symbols
% /src/digicoms/
% -------------------------------------------------------------------------
words_dec = demapping(symbs_cx,constellation);
% Demapping

% -------------------------------------------------------------------------
% diffdec_qam
% Differential decoding for square mQAM
% /src/digicoms/
% -------------------------------------------------------------------------
nsymbols = 2^20;
m = 64;
esn0_db = 20;
bits = generate_binary(nsymbols,m);
symbs = diffenc_qam(bits,m);
symbs_rx = add_awgn(symbs,esn0_db);
[bits_rx,symbs_diff] = diffdec_qam(symbs_rx,m);
ber = calc_ber(bits_rx,bits(log2(m)+1:end)) 

% -------------------------------------------------------------------------
% diffdec_qpsk
% Differential decoding for QPSK
% /src/digicoms/
% -------------------------------------------------------------------------
nsymbols = 2^10;
m = 4;
esn0_db = 20;
bits = generate_binary(nsymbols,m);
[words_dec, words_bin] = conv_bin2dec(bits,log2(m));
symbs = diffenc_qpsk(words_dec);
symbs_rx = add_awgn(symbs,esn0_db);
[bits_rx, symbs_diff] = diffdec_qpsk(symbs_rx,'hard');
ber = calc_ber(bits_rx,bits(3:end));

% -------------------------------------------------------------------------
% diffenc_qam
% Generation of differentially encoded symbols for square m-QAM modulation
% /src/digicoms/
% -------------------------------------------------------------------------
symbs = diffenc_qam(bits,m);
% Differential encoding for square QAM

% -------------------------------------------------------------------------
% diffenc_qpsk
% Generation of differentially encoded symbols for QPSK modulation
% /src/digicoms/
% -------------------------------------------------------------------------
symbs = diffenc_qpsk(words_dec);
% Differential encoding for QPSK

% -------------------------------------------------------------------------
% diffenc_qpsk_nextquadrant
% Determine quadrant index for QPSK differential encoding
% /src/digicoms/
% -------------------------------------------------------------------------
iquad(iword) = diffenc_qpsk_nextquadrant(iquad(iword - 1),words_dec(iword - 1));
% Determine quadrant index based on 2-bit word value and previous quadrant index

% -------------------------------------------------------------------------
% extract_snr_constellation
% Extract the SNR from data-aided clusters of a constellation
% /src/digicoms/
% -------------------------------------------------------------------------
snr_db = extract_snr_constellation(cluster_da);
% Extract SNR for each data-aided cluster of the constellation

% -------------------------------------------------------------------------
% generate_binary
% Generation of uniformly distributed binary random data
% /src/digicoms/
% -------------------------------------------------------------------------
bits_bin = generate_binary(nsymbols,m);
% Generation of uniformly distributed binary random data

% -------------------------------------------------------------------------
% mapping
% Mapping decimal words to complex symbols
% /src/digicoms/
% -------------------------------------------------------------------------
symbs = mapping(words_dec,constellation);
% Mapping

% -------------------------------------------------------------------------
% normalise_constellation
% Normalisation of received constellation using energy per symbol
% /src/digicoms/
% -------------------------------------------------------------------------
symbs_rx = normalise_constellation(symbs_rx,norm_es);
% Normalise constellation

% -------------------------------------------------------------------------
% plot_constellation
% Quickly plot digital constellation diagram
% /src/digicoms/
% -------------------------------------------------------------------------
constellation_type = 'plain';%'heat';'cluster';
constellation_name = 'received constellation';
plot_constellation(symbs,constellation_type,constellation_name,[-5:1:5]);
% Plot constellation

% -------------------------------------------------------------------------
% plot_constellation_mapping
% Plot original digital constellation with mapping
% /src/digicoms/
% -------------------------------------------------------------------------
params.name = ['Constellation for '];
params.labels.dec = 1; 
params.labels.bin = 1;
params.axes = 0;
params.save = 0;
plot_constellation_mapping(constellation,params); 
% Plot original constellation

% -------------------------------------------------------------------------
% plot_constellation_mapping
% Reference SEP versus Es/N0 curves over AWGN for various constellations 
% /src/digicoms/
% -------------------------------------------------------------------------
[pse, status] = sep_awgn_reference(esn0_db,mod,m);
% Theoretical symbol-error-probability versus SNR

% -------------------------------------------------------------------------
% serial_parallel
% Serial-to-parallel conversion
% /src/digicoms/
% -------------------------------------------------------------------------
streams_out = serial_parallel(bits_bin,nstreams); 
% Serial-to-parallel conversion





% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp                                                                       dsp
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp/adc-dac                                                               dsp/adc-dac
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% dac_zoh
% Basic zero-order hold interpolation
% /src/dsp/adc-dac/
% -------------------------------------------------------------------------
sig = dac_zoh(samples,ndac);
% Zero-order hold


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp/digital-coherent                                                      dsp/digital-coherent
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% cfo_diffphase
% CFO estimation for M-PSK by differential phase algorithm (Leven)
% /src/dsp/digital-coherent/
% -------------------------------------------------------------------------
tk = [0:1:length(symbs) - 1]/symbol_rate;
mpower = 4;
sample_delay = 1;
block_length = 1000;
[cfo_estimate, cfo_estimation_delay] = cfo_diffphase(symbs,mpower,sample_delay,block_length);
symbs_comp = dsp_delay(symbs,cfo_estimation_delay).*exp(-1j*2*pi*cfo_estimate.*dsp_delay(tk,cfo_estimation_delay)*symbol_rate);
symbs_comp = symbs_comp([2*cfo_estimation_delay + 1:end]);

nblocks = floor(length(symbs)/block_length);
cfo_estimate = cfo_diffphase(symbs,mpower,sample_delay,block_length,'block');
symbs_comp = symbs(1:nblocks*block_length).*exp(-1j*2*pi*cfo_estimate.*tk(1:nblocks*block_length)*symbol_rate);

% -------------------------------------------------------------------------
% cpr_bps
% Carrier phase estimation using the blind phase search (BPS) algorithm
% /src/dsp/digital-coherent/
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% iq_gsop
% Gram-Schmidt orthogonalization procedure for IQ imbalance compensation
% /src/dsp/digital-coherent/
% -------------------------------------------------------------------------
[samps_cmp_r,samps_cmp_i] = iq_gsop(real(samps),imag(samps));
samps_cmp = samps_cmp_r + 1i*samps_cmp_i; 
samps_cmp = normalise_constellation(samps_cmp,norm_es);

% -------------------------------------------------------------------------
% iq_lop
% Löwdin orthogonalization procedure for IQ imbalance compensation
% /src/dsp/digital-coherent/
% -------------------------------------------------------------------------
[samps_cmp_r,samps_cmp_i] = iq_lop(real(samps),imag(samps));
samps_cmp = samps_cmp_r + 1i*samps_cmp_i; 
samps_cmp = normalise_constellation(samps_cmp,norm_es);


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp/general                                                               dsp/general
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------


% -------------------------------------------------------------------------
% dsp_delay
% Digital delay
% /src/dsp/general/
% -------------------------------------------------------------------------
a = linspace(1,10,10);
b = dsp_delay(a,5);
c = dsp_delay(a,5,[1 2 3 4 5]);
aa = linspace(1,50,50);
aa = reshape(aa,5,10);
bb = dsp_delay(aa,3);
cc = dsp_delay(aa,3,[9 10 11]);
dd = dsp_delay(aa,3,[9 10 11;12 13 14;15 16 17;18 19 20;21 22 334]);
% Digital delay

% -------------------------------------------------------------------------
% dsp_farrow
% Implementation of Farrow structure, e.g. for interpolation
% /src/dsp/general/
% -------------------------------------------------------------------------
% Linear interpolator (2 basepoints - 1 sample delay):
nfilters = 2;
ntaps = 2;
c0 = [0, 1];
c1 = [1, -1];
c = [c0;c1].';

% Parabolic interpolator (3 basepoints - 1 sample delay):
nfilters = 3;
ntaps = 3;
c0 = [0,   1,  0];
c1 = [0.5, 0,  -0.5];
c2 = [0.5, -1, 0.5];
c = [c0;c1;c2].';

% Cubic interpolator (4 basepoints - 2 sample delay): 
nfilters = 4;
ntaps = 4;
c0 = [0,    0,     1,    0];
c1 = [-1/6, 1,     -1/2, -1/3];
c2 = [0,   1/2,    -1,   1/2];
c3 = [1/6, -1/2,   1/2,  -1/6]; 
c = [c0;c1;c2;c3].';

% Parabolic interpolator (4 basepoints - 2 sample delay)
ntaps = 4;
nfilters = 3;
alpha = 0.5;
c = [0,  -alpha,     alpha;
    0,   alpha+1,  -alpha;
    1,   alpha-1   -alpha;
    0,  -alpha,     alpha];

y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
% Farrow interpolator


% -------------------------------------------------------------------------
% dsp_fir_design_frequency_sampling
% FIR filter coefficients determination by frequency sampling
% /src/dsp/general/
% -------------------------------------------------------------------------
roll_off = 0.1; 
% Raised-cosine roll-off factor
ts = 1/symbol_rate; 
% Symbol rate, in baud
nos = 4;
% Oversampling factor. Even number.
fsa = nos*symbol_rate;
% Sampling frequency
nlengthsymbs = 16;
% Filter depth, in terms of number of symbols
% 16 means 8 symbols on each side
% Even number
ntaps = nlengthsymbs*nos + 1;
% FIR filter length. Ensure it is an odd number
freq_sampling = linspace(0,fsa/2,1000);
spectrum_sampling = calc_rc_spectrum(freq_sampling,ts,roll_off);
% Analog spectrum of raise-cosine filter
h = dsp_fir_design_frequency_sampling(freq_sampling,spectrum_sampling,fsa,ntaps,1);
% Tap coefficients 

% -------------------------------------------------------------------------
% dsp_fir_frequency_response
% Calculate the frequency response of a FIR filter from its impulse response / tap coefficients
% /src/dsp/general/
% -------------------------------------------------------------------------
[wh,H] = dsp_fir_frequency_response(h);
% Filter frequency response

% -------------------------------------------------------------------------
% dsp_fir_linear
% Linear buffer implementation of a finite impulse response (FIR) filter.
% /src/dsp/general/
% -------------------------------------------------------------------------
[y,z] = dsp_fir_linear(x,b,z);
% FIR filter

% -------------------------------------------------------------------------
% dsp_overlap_save
% Digital filtering in the frequency domain using overlap-and-save method
% /src/dsp/general/
% -------------------------------------------------------------------------
y = dsp_overlap_save(x,h,L);
% Overlap-and-save fitering


% -------------------------------------------------------------------------
% dsp_upsample
% Upsample by an integer factor by adding L - 1 zeros between samples
% /src/dsp/general/
% -------------------------------------------------------------------------
y = dsp_upsample(x,L);
% Upsampling



% -------------------------------------------------------------------------
% dsp_window
% Window functions for DSP applications
% /src/dsp/general/
% -------------------------------------------------------------------------
params_window.type = 'kaiser';%'hamming';'hann';'bartlett';'rectangular';
params_window.length = 15;
params_window.symmetry = 'periodic';%'symmetric';
params_window.beta = 2.4*pi;% For Kaiser window.
w = dsp_window(params_window); 
% Calculate window function

% -------------------------------------------------------------------------
% dsp_window_properties
% Calculate some properties of DSP windows
% /src/dsp/general/
% -------------------------------------------------------------------------
[nenb,spectrum,norm_ang_freq] = dsp_window_properties(w,fft_length);
% Properties of window function






% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp/spectral-analysis                                                     % dsp/spectral-analysis
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% harmonics_first_nyquist
% Calculate frequencies of aliases of harmonics in 1st Nyquist zone
% /src/dsp/spectral-analysis
% -------------------------------------------------------------------------
[freq,cf,cfs] = harmonics_first_nyquist(f,fs,nmax);
% Calculate frequencies of aliases of harmonics in 1st Nyquist zone

% -------------------------------------------------------------------------
% psd_welch
% Power spectral density estimation using the Welch periodogram method
% /src/dsp/spectral-analysis
% -------------------------------------------------------------------------
block_size = 100;
overlap_samples = 50;
fft_length = 640;
window_type = 'hamming';
[psd,norm_freq] = psd_welch(data,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');

% -------------------------------------------------------------------------
% spectrum_one_sided
% Quick calculation of one-sided spectrum of a signal
% /src/dsp/spectral-analysis
% -------------------------------------------------------------------------
[f,X1] = spectrum_one_sided(x,fs);
[f,X1] = spectrum_one_sided(x,fs,2^10);
[f,X1] = spectrum_one_sided(x,fs,2^10,1);
[f,X1] = spectrum_one_sided(x,fs,[],1); 
% One-sided spectrum





% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% electrical                                                                electrical
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% calc_enb
% Equivalent noise bandwidth of optical or electrical filter
% /src/electrical/
% -------------------------------------------------------------------------
tf = elec_tf_elpf(params_elpf,frequency_array);
enb = calc_enb_tf('lowpass',tf,max(abs(tf)),df);
% Equivalent noise bandwidth

% -------------------------------------------------------------------------
% elec_coder_manchester
% Manchester encoder without rise time
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_coder_manchester(bit_pattern,nsamples_per_symbol); 
% Manchester encoder

% -------------------------------------------------------------------------
% elec_coder_nrz
% Electrical NRZ encoder without rise time
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_coder_nrz(bit_pattern,nsamples_per_symbol); 
% NRZ encoder without rise time

% -------------------------------------------------------------------------
% elec_coder_rz
% Electrical RZ encoder 
% /src/electrical/
% -------------------------------------------------------------------------
params_rz.duty_cycle = 0.5;
params_rz.rise_time = 1/symbol_rate/8;
params_rz.normalisation = 1;
[sig,params_rz] = elec_coder_rz(bit_pattern,nsamples_per_symbol,params_rz); 
% RZ signal

% -------------------------------------------------------------------------
% elec_dc_block
% Electrical DC block
% /src/electrical/
% -------------------------------------------------------------------------
[sig, dclevel] = elec_dc_block(sig);
% Electrical DC block

% -------------------------------------------------------------------------
% elec_bias
% Applies bias and set peak-to-peak value to electrical signal
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_bias(sig,dc,pp);
% Set bias and peak-to-peak value


% -------------------------------------------------------------------------
% elec_elpf
% Electrical low-pass filter
% /src/electrical/
% -------------------------------------------------------------------------
params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
params_elpf.order = 4;
params_elpf.f3dB = 0.75*symbol_rate;
params_elpf.roll_off = 0.6;            % for 'raised_cosine' and 'root_raised_cosine'
params_elpf.symbol_rate = symbol_rate; % for 'raised_cosine' and 'root_raised_cosine'
sig = elec_elpf(sig,params_elpf); 
% Electrical low-pass filter
% For backward compatibility only.

% -------------------------------------------------------------------------
% elec_filter
% Electrical filter
% /src/electrical/
% -------------------------------------------------------------------------
params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
params_elpf.order = 4;
params_elpf.f3dB = 0.75*symbol_rate;
params_elpf.roll_off = 0.15;
params_elpf.symbol_rate = symbol_rate;
params_elpf.hold_time = 1/symbol_rate/2;
params_elpf.samples_per_symbol = nsamples_per_symbol;
tf = elec_tf_elpf(params_elpf,frequency_array);
sig = elec_filter(sig,tf);
% Electrical filter

% -------------------------------------------------------------------------
% elec_modulator
% Electrical modulator using convolution of impulses with pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
params_elecmod.pulse_shape = 'rc';%'rrc';'sinc';'sech';'gaussian';
params_elecmod.roll_off = 1;
params_elecmod.symbol_rate = symbol_rate;
params_elecmod.fwhm = 1/symbol_rate/6;
params_elecmod.order = 1;
params_elecmod.normalisation = 'peak';
sig = elec_modulator(symbols,params_elecmod); 
% Electrical modulator

% -------------------------------------------------------------------------
% elec_pulse_gauss
% Electrical Gaussian pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_gauss(time_array,pulse_position,pulse_duration,pulse_order); 
% Gaussian pulse

% -------------------------------------------------------------------------
% elec_pulse_rc
% Raised-cosine pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_rc(time_array,pulse_position,1/symbol_rate,alpha);
% Raised-cosine pulse

% -------------------------------------------------------------------------
% elec_pulse_rectangular
% Rectangular pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_rectangular(time_array,pulse_position,pulse_duration);
% Rectangular pulse

% -------------------------------------------------------------------------
% elec_pulse_rrc
% Root-raised-cosine pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_rrc(time_array,time_array(nsamples/2),1/symbol_rate,alpha); 
% Root-raised-cosine pulse

% -------------------------------------------------------------------------
% elec_pulse_sech
% Electrical sech pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_sech(time_array,pulse_position,pulse_duration); 
% Electrical hyperbolic secant pulse

% -------------------------------------------------------------------------
% elec_pulse_sinc
% sinc pulse shape
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_sinc(time_array,pulse_position,1/symbol_rate); 
% sinc pulse shape

% -------------------------------------------------------------------------
% elec_pulse_sequence_nrz
% Electrical NRZ pulse sequence generation
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% NRZ sequence generation

% -------------------------------------------------------------------------
% elec_resync
% Resynchronise electrical signal to a reference signal
% /src/electrical/
% -------------------------------------------------------------------------
[sig, delay_samples] = elec_resync(sig,sig_ref);
% Resynchronise electrical signal

% -------------------------------------------------------------------------
% elec_rise_time
% Apply rise time to a rectangular electrical signal
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_rise_time(sig,rise_time);
% Rise time

% -------------------------------------------------------------------------
% elec_sinusoidal
% Electrical sinusoidal signal generation
% /src/electrical/
% -------------------------------------------------------------------------
params_rf.frequency = symbol_rate;
params_rf.phase = 0;
params_rf.vpp = 1.0;
params_rf.vdc = 0;
sig = elec_sinusoidal(params_rf); 
% Electrical sinusoidal signal generation

% -------------------------------------------------------------------------
% elec_tf_elpf
% Transfer functions of electrical low-pass filters
% /src/electrical/
% -------------------------------------------------------------------------
params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
params_elpf.order = 4;
params_elpf.f3dB = 0.75*symbol_rate;
params_elpf.roll_off = 0.15;
params_elpf.symbol_rate = symbol_rate;
params_elpf.hold_time = 1/symbol_rate/2;
params_elpf.samples_per_symbol = nsamples_per_symbol;
tf = elec_tf_elpf(params_elpf,frequency_array);
% Transfer functions of electrical low-pass filters

% -------------------------------------------------------------------------
% elec_tf_matched
% Transfer function of matched electrical filter
% /src/electrical/
% -------------------------------------------------------------------------
tf = elec_tf_matched(sig_pulse,dt);
% Transfer function of matched electrical filter
sig = elec_filter(sig,tf);



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% fibres                                                                    fibres
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------


% -------------------------------------------------------------------------
% opt_nlse_scalar_basic
% Basic scalar nonlinear Schroedinger equation
% /src/fibres/
% -------------------------------------------------------------------------
params_fibre.nonlinear_coefficient = 1.3e-3;% nonlinear coefficient, in 1/W/m
params_fibre.dispersion = 16;% dispersion, in ps/nm/km
params_fibre.dispersion_slope = 0.058;% dispersion slope, in ps/nm2/km
params_fibre.dispersion_curvature = 0;% dispersion curvature, in ps/nm3/km
params_fibre.dispersion_spec_frequency = reference_frequency;
params_fibre.loss = 0.2;% loss, in dB/km
params_fibre.length = 10e3;% fibre length, in m
numparams_fibre.max_step_size = 1;% maximum step size, in m
numparams_fibre.max_phase_shift = 1e-3;% maximum nonlinear phase shift, in radians
params_fibre.beta_coefficients = dispersion_conv_d_beta([params_fibre.dispersion params_fibre.dispersion_slope params_fibre.dispersion_curvature],'to_beta','eng','si',params_fibre.dispersion_spec_frequency);
params_fibre.loss_alpha = conv_loss_lin_log(params_fibre.loss);
sig = opt_nlse_scalar_basic(sig,params_fibre,numparams_fibre);
% Basic scalar nonlinear Schroedinger equation


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% filters                                                                   filters
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% extract_dispersion_from_tf
% Group delay and dispersion from optical filter complex transfer function
% /src/filters/
% -------------------------------------------------------------------------
% save_dispersion.status = 0;
% save_dispersion.file_name = 'filter_tf.dat';
% [phase,fgd,gd,fgvd,gvd] = extract_dispersion_from_tf(frequency_array,tf,'si',save_dispersion);
% Extract filter GD and GVD

% -------------------------------------------------------------------------
% opt_filter
% Optical filter
% /src/filters/
% -------------------------------------------------------------------------
sig = opt_filter(sig,tf);
% Optical filter

% -------------------------------------------------------------------------
% opt_mrr
% Micro-ring resonator filter in add-drop configuration
% /src/filters/
% -------------------------------------------------------------------------
params_mrr.field_round_trip_loss = 0.96;
params_mrr.power_coupling_1 = 0.63;
params_mrr.power_coupling_2 = params_mrr.power_coupling_1;
params_mrr.centre_frequency = 0;
params_mrr.fsr = 100e9;
params_mrr.visualiser_status = 0;
params_mrr.save_tf.status = 0;
params_mrr.save_tf.file_name = 'mrr_tf.dat';
sig_add = opt_no_sig;
[sig_through,sig_drop] = opt_mrr(sig,sig_add,params_mrr); 
% Add-drop MRR 

% -------------------------------------------------------------------------
% opt_mzdi
% Mach-Zehnder delay interferometer
% /src/filters/
% -------------------------------------------------------------------------
params_mzdi.delay_target = 1/symbol_rate;
params_mzdi.coupling_ratio_in = 0.5;
params_mzdi.coupling_ratio_out = 0.5;
params_mzdi.mode = 'tuned';%'general';
params_mzdi.phase_shift = 0;
[sig21,sig22] = opt_mzdi(sig,opt_nosig,params_mzdi); 
% Mach-Zehnder delay interferometer

% -------------------------------------------------------------------------
% opt_tf_fbg
% Transfer function of fibre Bragg grating
% /src/filters/
% -------------------------------------------------------------------------
params_fbg.length = 0.03;                 % Grating length, in m.
params_fbg.centre_frequency = 193.1e12;   % Centre frequency of the grating, in Hz.
params_fbg.n0 = 1.45;                     % Effective index of the fibre.
params_fbg.dn = 1.0e-4;                   % Refractive index change.
params_fbg.m = 0;                         % Control of average index change.
params_fbg.chirp_rate = 0;                % Laser linear chirp rate, in 1/m.
params_fbg.apodisation.type = 'uniform';     
params_fbg.apodisation.eta = 0;
params_fbg.apodisation.fwhm = 0;
params_fbg.apodisation.profile = 0;
numparams_fbg.nsections = 1000;           % Number of uniform grating sections.
params_fbg.phase_shift = zeros(1,numparams.nsections);% Position of phase shifts.
freq = frequency_array + reference_frequency;
tf = opt_tf_fbg(freq,params_fbg,numparams_fbg);
% Transfer function of fibre Bragg grating

% -------------------------------------------------------------------------
% opt_tf_fp
% Transfer function of Fabry-Perot optical filter
% /src/filters/
% -------------------------------------------------------------------------
params_fp.mode = 'transmission';%'reflection';
params_fp.finesse = 100;
params_fp.fsr = 100e9;
params_fp.centre_frequency = 0;
tf = opt_tf_fp(frequency_array,params_fp);
% Transfer function of Fabry-Perot optical filter

% -------------------------------------------------------------------------
% opt_tf_gt
% Transfer function of Gires-Tournois interferometer
% /src/filters/
% -------------------------------------------------------------------------
params_gt.reflectivity = 0.9;
params_gt.fsr = 100e9;
params_gt.centre_frequency = 0;
tf = opt_tf_gt(frequency_array,params_gt);
% Transfer function of Gires-Tournois interferometer

% -------------------------------------------------------------------------
% opt_tf_matched
% Transfer function of optical matched filter
% /src/filters/
% -------------------------------------------------------------------------
tf = opt_tf_matched(pulse,freq);
% Transfer function of optical matched filter

% -------------------------------------------------------------------------
% opt_tf_mrr
% Transfer function of add-drop micro-ring resonator filter
% /src/filters/
% -------------------------------------------------------------------------
params_mrr.field_round_trip_loss = 0.96;
params_mrr.power_coupling_1 = 0.63;
params_mrr.power_coupling_2 = params_mrr.power_coupling_1;
params_mrr.centre_frequency = 0;
params_mrr.fsr = 100e9;
vis_mrr.visualiser_status = 1;
vis_mrr.save_tf.status = 0;
vis_mrr.save_tf.file_name = 'mrr_tf.dat';
tf = opt_tf_mrr(frequency_array,params_mrr,vis_mrr); 
% Transfer function of add-drop micro-ring resonator filter

% -------------------------------------------------------------------------
% opt_tf_mrr_crow
% Coupled microring resonator (CROW) scattering matrix 
% /src/filters/
-------------------------------------------------------------------------
params_crow.nrings = 3;
params_crow.power_coupling = [0.45 0.5 0.3 0.45]';
params_crow.radius = [10e-6 20e-6 10e-6];
params_crow.phase_shift = [0 0 0];
params_crow.neff = [2.54 2.54 2.54]';
params_crow.loss_log = [1e-2 1e-2 1e-2]';
M = opt_tf_mrr_crow(frequency_array,params_crow);
tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';
% Transfer function of coupled MRRs

% -------------------------------------------------------------------------
% opt_tf_obpf
% Transfer functions of some standard optical bandpass filters
% /src/filters/
% -------------------------------------------------------------------------
params_obpf.type = 'gaussian';%'none','rectangular_ideal','rectangular';
params_obpf.centre_frequency = 0;
params_obpf.bandwidth = 40e9;
params_obpf.order = 4;
params_obpf.attenuation_in_band = -1;
params_obpf.attenuation_out_band = -20;
tf = opt_tf_obpf(params_obpf,frequency_array);
% Optical bandpass filter transfer function


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% logical                                                                   logical
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% logical_adapt_binary_sequence
% Adapt binary sequence to simulation time window
% /src/logical/
% -------------------------------------------------------------------------
seq = logical_adapt_binary_sequence(seq,nsymbols); 
% Adapt binary sequence to simulation time window

% -------------------------------------------------------------------------
% logical_alternate
% Generation of alternate binary sequence
% /src/logical/
% -------------------------------------------------------------------------
seq = logical_alternate(alternate0,nsymbols);
% Generation of alternate binary sequence

% -------------------------------------------------------------------------
% logical_differential_encoder_binary
% Differential encoder for DPSK or duobinary signal generation
% /src/logical/
% -------------------------------------------------------------------------
bit_pattern_diffenc = logical_differential_encoder_binary(bit_pattern);
% Differential encoder

% -------------------------------------------------------------------------
% logical_prbs
% Generation of pseudo-random binary sequences
% /src/logical/
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';%'de_bruijn';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate PRBS

% Other generating polynomials:
% params_prbs.order = 10;
% params_prbs.poly = [10 3 0];%[10 7 0];
% params_prbs.seed = [0 1 1 0 1 1 1 0 0 1];

% params_prbs.order = 11;
% params_prbs.poly = [11 9 0];%[11 8 5 2 0];
% params_prbs.seed = [0 1 1 0 1 1 1 0 0 1 0];


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% modes                                                                     modes
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% beam_gauss
% Scalar field of a fundamental Gaussian beam
% /src/modes/
% -------------------------------------------------------------------------
position = [0 0];
lambda = 1550e-9;
k = 2*pi/lambda;
w0 = 10e-6;
z = 0;
power = 1e-3;
u = sqrt(power)*beam_gauss(position,k,w0,z);
% Create Gaussian beam

% -------------------------------------------------------------------------
% create_space_grid
% Create 2D space grid in the transverse plane
% /src/modes/
% -------------------------------------------------------------------------
xrange = [-15e-6, 15e-6];  
yrange = [-15e-6, 15e-6]; 
nxpoints = 2001;
nypoints = 2001;
space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
[space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 

% -------------------------------------------------------------------------
% fibre_si_b
% Calculate normalised propagation constant for circular-core step-index fibre
% /src/modes/
% -------------------------------------------------------------------------
params_fibre.a = 4.7e-6;
params_fibre.n1 = 1.4628;
params_fibre.n2 = 1.4600;
lambda = 1550e-9;
mode_type = 'LP';
mode_l = 0;
V = 2*pi*params_fibre.a*sqrt(params_fibre.n1^2 - params_fibre.n2^2)/lambda;
[b, nmodes] = fibre_si_b(mode_type,mode_l,V,params_fibre);

% -------------------------------------------------------------------------
% fibre_si_lp_field
% Calculate scalar field of the LPlm mode of a step-index fibre
% /src/modes/
% -------------------------------------------------------------------------
xrange = [-1, 1]*50e-6;  
yrange = [-1, 1]*50e-6; 
nxpoints = 2001;
nypoints = 2001;
global space_grid
space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
[space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 

params_fibre.n1 =  1.45;
params_fibre.n2 =  1.4361;
params_fibre.a = 25e-6;
mode_type = 'LP';
mode_l = 1;             % Index l of mode LPlm
mode_m = 1;             % Index m of mode LPlm. Check that mode_m <= nmodes
lambda = 850e-9; 
V = 2*pi*params_fibre.a*sqrt(params_fibre.n1^2 - params_fibre.n2^2)/lambda;
[b,nmodes] = fibre_si_b('LP',mode_l,V,params_fibre);
% Normalised propagation constant
field = fibre_si_lp_field(params_fibre,mode_l,b(mode_m),V,'even'); 
% Mode field

visparams.limit_radius = Inf;
visparams.show_core_limit = 0;
visparams.show_core_linewidth = 1;
visparams.save = 0;
visparams.colormap = 'jet';%'hot';
visparams.name = ['Mode field distribution for ' mode_type num2str(mode_l) num2str(mode_m)];
fibre_si_plot_mode(field,params_fibre,visparams);

% -------------------------------------------------------------------------
% fibre_si_plot_mode
% Plot field or intensity distribution (irradiance pattern) of fibre mode
% /src/modes/
% -------------------------------------------------------------------------
visparams.limit_radius = Inf;
visparams.show_core_limit = 0;
visparams.save = 0;
visparams.colormap = 'jet';%'hot';
visparams.name = 'Mode field distribution';
fibre_si_plot_mode(field,params_fibre,visparams);

% -------------------------------------------------------------------------
% plot_irradiance
% Plot intensity distribution
% /src/modes/
% -------------------------------------------------------------------------
visparams.display = '2d';% '2d_slicex','2d_slicey','2d_slicexy'
visparams.name = 'Intensity distribution for xxx beam';
visparams.xslice = 0;
visparams.yslice = 10e-6;
hfig = plot_irradiance(I,visparams);


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% modulators                                                                modulators
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% mod_iq
% IQ optical modulator (single-polarization)
% /src/modulators/
% -------------------------------------------------------------------------
params_iq.vpi_mzm = 1;
params_iq.vpi_pm = 1;
params_iq.bias_i1 = params_iq.vpi_mzm;
params_iq.bias_i2 = 0;
params_iq.bias_q1 = params_iq.vpi_mzm;
params_iq.bias_q2 = 0;
params_iq.split_i = 0.5;
params_iq.split_q = 0.5;
params_iq.loss_i = 0;
params_iq.loss_q = 0;
v_ps = params_iq.vpi_pm/2;
sig_i1 = params_iq.vpi_mzm*(nrz_data_sig_i - 0.5);
sig_i2 = -params_iq.vpi_mzmM*(nrz_data_sig_i - 0.5);
sig_q1 = params_iq.vpi_mzm*(nrz_data_sig_q - 0.5);
sig_q2 = -params_iq.vpi_mzm*(nrz_data_sig_q - 0.5);
sig = mod_iq(sig,sig_i1,sig_i2,sig_q1,sig_q2,v_ps,params_iq);
% IQ modulator


% -------------------------------------------------------------------------
% mod_linear
% Linear intensity modulator
% /src/modulators/
% -------------------------------------------------------------------------
extinction_ratio_db = 30;
sig = mod_linear(sig,nrz_data_sig,extinction_ratio_db); 
% Linear intensity modulator


% -------------------------------------------------------------------------
% mod_mzm
% Mach-Zehnder modulator (MZM)
% /src/modulators/
% -------------------------------------------------------------------------
vpi = 1.0;
driving_signal_1 = vpi/2*(nrz_data_sig-0.5);
driving_signal_2 = -vpi/2*(nrz_data_sig-0.5);
bias_1 = 1.5*vpi;
bias_2 = 0;
split_in = 0.5;
split_out = 0.5;
loss = 0;
sig = mod_mzm(sig,driving_signal_1,driving_signal_2,bias_1,bias_2,vpi,split_in,split_out,loss);
% Mach-Zehnder modulator


% -------------------------------------------------------------------------
% mod_mzm_pulse_carver
% Mach-Zehnder modulator-based pulse carver 
% /src/modulators/
% -------------------------------------------------------------------------
optical_clock_frequency = symbol_rate;
optical_clock_duty_cycle = 'rz33';%'rz50';%'rz67';
sig = mod_mzm_pulse_carver(sig,optical_clock_frequency,optical_clock_duty_cycle);
% Mach-Zehnder modulator-based pulse carver


% -------------------------------------------------------------------------
% mod_pm
% Phase modulator (PM)
% /src/modulators/
% -------------------------------------------------------------------------
vpi = 1.0;
loss = 0;
sig = mod_pm(sig,nrz_data_sig,vpi,loss); 
% Phase modulator

% -------------------------------------------------------------------------
% mod_polm
% Polarisation modulator
% /src/modulators/
% -------------------------------------------------------------------------
vpi = 1.0;
loss = 0;
sig = mod_polm(sig,nrz_data_sig,vpi,loss); 
% Polarisation modulator



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical                                                                 numerical
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical/differentiation                                                 numerical/differentiation
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% num_diff
% pproximation of first or second-order derivatives by central finite differences
% /src/numerical/differentiation/
% -------------------------------------------------------------------------
[xn,dy] = num_diff(order,x,y);

% -------------------------------------------------------------------------
% num_diff_1d_pb
% Differentiate evenly spaced 1D data with periodic boundary conditions
% /src/numerical/differentiation/
% -------------------------------------------------------------------------
dfunc = num_diff_1d_pb(func,h);
% Differentiation



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical/functions                                                       numerical/functions
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% func_sinc
% Function sinc(x) = sin(x)/x
% /src/numerical/functions/
% -------------------------------------------------------------------------
y = func_sinc(x);

% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical/integration                                                     % numerical/integration
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% num_int1d_simpson
% Numerical integration using Simpson rule
% /src/numerical/integration/
% -------------------------------------------------------------------------
integral = num_int1d_simpson(func,dx);

% -------------------------------------------------------------------------
% num_int2d_simpson
% Double integral of function defined on Cartesian grid by Simpson method
% /src/numerical/integration/
% -------------------------------------------------------------------------


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical/matrices                                                        numerical/matrices
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% prod_mm
% Product of time or frequency dependent 2x2 matrices
% /src/numerical/matrices/
% -------------------------------------------------------------------------
C = prod_mm(A,B);
% Matrix-matrix product

% -------------------------------------------------------------------------
% prod_mv
% Product of time or frequency dependent 2-element vector and 2x2 matrix
% /src/numerical/matrices
% This function was initially created to manipulate Jones vector and is 
% identical to jones_prod_mv under /src/polarization
% We just created this function since its scope of use is much larger while
% keeping a self-contained polarization package.
% We are aware this is very bad practice.
% -------------------------------------------------------------------------
B = prod_mv(M,A);
% Vector-matrix product


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% numerical/transforms                                                      numerical/transforms
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% num_ft
% Numerical estimation of the continuous-time Fourier transform based on DFT
% /src/numerical/transforms/
% -------------------------------------------------------------------------
F = num_ft(f,dt);
% For time interval [0 T]
F = num_ft(f,dt,0);
% For time interval [-T/2 T/2]
F = num_ft(f,dt,(T1 + T2)/2);
% For time interval [T1 T2]
% Fourier transform

% -------------------------------------------------------------------------
% num_ift
% Numerical estimation of the inverse continuous Fourier transform based on IDFT
% /src/numerical/transforms/
% -------------------------------------------------------------------------
f = num_ift(F,dt);
% Inverse Fourier transform






% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% passive                                                                   passive
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% opt_att
% Optical attenuator 
% /src/passive/
% -------------------------------------------------------------------------
params_att.mode = 'attenuation';%'power';
params_att.target = 10;  % attenuation in dB or total output power in dBm
sig = opt_att(sig,params_att);
% Optical attenuator


% -------------------------------------------------------------------------
% opt_combiner_y_junction
% Optical Y-junction combiner with 50% combining power ratio
% /src/passive/
% -------------------------------------------------------------------------
sig = opt_combiner_y_junction(sig1,sig2);
% Optical Y-junction combiner with 50% combining power ratio



% -------------------------------------------------------------------------
% opt_coupler_2x2
% 2x2 optical directional coupler
% /src/passive/
% -------------------------------------------------------------------------
[sig21,sig22] = opt_coupler_2x2(sig11,sig12,'lin',0.5); 
% 2x2 optical directional coupler


% -------------------------------------------------------------------------
% opt_delay
% Optical delay
% /src/passive/
% -------------------------------------------------------------------------
include_phase_shift = 0;%1;
delay_target = 100e-12;
[sig,delay_actual] = opt_delay(sig,include_phase_shift,delay_target); 
% Optical delay


% -------------------------------------------------------------------------
% opt_dispersion
% Linear and lossless dispersive element
% /src/passive/
% -------------------------------------------------------------------------
params_dispersion.dispersion = 17;           % Dispersion in ps/nm/km
params_dispersion.dispersion_slope = 0.058;  % Dispersion slope in ps/nm^2/km
params_dispersion.dispersion_curvature = 0;  % Dispersion curvature, in ps/nm^3/km
params_dispersion.dispersion_spec_frequency = reference_frequency;
z = 1e3; 
sig = opt_dispersion(sig,params_dispersion,z); 
% Linear and lossless dispersive element


% -------------------------------------------------------------------------
%  opt_hybrid90_2x4
% 90 degree 2x4 optical hybrid
% /src/passive/
% -------------------------------------------------------------------------
[sigout1,sigout2,sigout3,sigout4] = opt_hybrid90_2x4(sig,sig_lo); 
% 90 degree 2x4 optical hybrid

% -------------------------------------------------------------------------
% opt_otdm_mux
% OTDM multiplexer
% /src/passive/
% -------------------------------------------------------------------------
otdm_mux_multiplication_factor = 4;
otdm_mux_mode = 'serial';% 'serial_random';'parallel';
[sig,otdm_mux_delays] = opt_otdm_mux(sig,symbol_rate,...
      otdm_mux_multiplication_factor,params_prbs.order,otdm_mux_mode);
% OTDM multiplexer

% -------------------------------------------------------------------------
% opt_phase_shift
% Optical phase-shifter
% /src/passive/
% -------------------------------------------------------------------------
 sig = opt_phase_shift(sig,phase_shift); 
 
 
 
% -------------------------------------------------------------------------
% opt_splitter_y_junction
% Optical Y-junction splitter with 50% combining power ratio
% /src/passive/
% -------------------------------------------------------------------------
[sig1,sig2] = opt_splitter_y_junction(sig);
% Y-junction splitter with 50% power splitting ratio



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% performance                                                               performance
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
 
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% performance/ber                                                           performance/ber
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
 
% -------------------------------------------------------------------------
% ber_gauss
% BER estimation using Gaussian approximation
% /src/performance/ber/
% -------------------------------------------------------------------------
params_bergauss.ignore_bits_start = 0;
params_bergauss.distribution = 'gauss';
params_bergauss.threshold_mode = 'optimum';%'optimum_search';%'fixed';
params_bergauss.threshold = 0.5;
params_bergauss.sample_mode = 'optimum';%'fixed';
params_bergauss.sample_index = [ ];
[BER,threshold,sample] = ber_gauss(sig,bit_pattern,params_bergauss); 
% BER calculation with Gaussian approximation


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% performance/eo                                                            performance/eo
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% eops
% Calculation of eye opening penalty at optimum sampling time
% /src/performance/eo/
% -------------------------------------------------------------------------
[eop,eo_norm,opt_sampling] = eops(abs(sig.x).^2,bit_pattern,params_eop.norm_power,-1);
% Eye opening calculation



% -------------------------------------------------------------------------
% eopw
% Calculation of eye opening penalty
% /src/performance/eo/
% -------------------------------------------------------------------------
params_eop.norm_power = 2*char_opt_average_power(sig);
sig = abs(sig.x).^2;
params_eop.method = 'eow';%'cr';
params_eop.bit_pattern = bit_pattern;
params_eop.eo_b2b = -1;%eyeop.eo_norm;
params_eop.display_eye = 1;
params_eop.eye_width_samples = 12;
params_eop.bits_to_ignore_start = 0;
params_eop.eye_display_name = 'Eye diagram';
eop = eopw(sig,params_eop); 
%Eye opening calculation



% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% polarisation                                                              polarisation
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% calc_dop 
% Calculate degree of polarization (DOP) from Stokes parameters
% /src/polarization/
% -------------------------------------------------------------------------
S = char_opt_stokes(sig);
Sav = calc_stokes_average(S,nsamples);
dop = calc_dop(Sav);
% Calculate degree of polarization


% -------------------------------------------------------------------------
% char_opt_stokes
% Characterization of instantaneous Stokes parameters
% /src/polarization/
% -------------------------------------------------------------------------
S = char_opt_stokes(sig);
% Calculate instantaneous Stokes parameters
S = S./S(1,:);
% Normalization by S0

% -------------------------------------------------------------------------
% jones_polarizer_x
% Jones matrix of polarizer // x
% /src/polarization/
% -------------------------------------------------------------------------
P = jones_polarizer_x(loss,extinction_ratio); 
% Jones matrix of polarizer // x

% -------------------------------------------------------------------------
% jones_polarizer_y
% Jones matrix of polarizer // y
% /src/polarization/
% -------------------------------------------------------------------------
P = jones_polarizer_y(loss,extinction_ratio); 
% Jones matrix of polarizer // x


% -------------------------------------------------------------------------
% plot_poincare_sphere
% Plot Poincare sphere
% /libs/poincare/
% -------------------------------------------------------------------------
S = char_opt_stokes(sig);
hfig = plot_poincare_sphere([file_name_core_figure '_poincare']);
plot3(S(2,:)./S(1,:),S(3,:)./S(1,:),S(4,:)./S(1,:),'LineStyle','-','LineWidth',3,'Color','r','Marker','none','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r')
% Trajectory
plot3(S(2)/S(1),S(3)/S(1),S(4)/S(1),'Marker','o','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
% SOP


% -------------------------------------------------------------------------
% pol_dgd
% Differential group delay (DGD)
% /src/polarization/
% -------------------------------------------------------------------------
dgd_target = 50e-12;
include_phase_shift = 0;
[sig,dgd_actual] = pol_dgd(sig,include_phase_shift,dgd_target);
% Differential group delay


% -------------------------------------------------------------------------
% pol_pbs 
% Polarization beam splitter (PBS)
% /src/polarization/
% -------------------------------------------------------------------------
pbs_il = 0;         % polarizer insertion loss, in dB
pbs_er = Inf;       % polarization extinction ratio, in dB
pbs_angle = 0;      % polarizer angle, in rad
[sig_x, sig_y] = pol_pbs(sig,pbs_angle,pbs_il,pbs_er); 
% Polarization beam splitter


% -------------------------------------------------------------------------
% pol_pdl
% Polarization dependent loss
% /src/polarization/
% -------------------------------------------------------------------------
il = 0;       % Insertion loss, in dB (positive number)
pdl = 3;      % Polarization dependent loss, in dB (positive number)
sig = pol_pdl(sig,insertion_loss,pdl); 
% Polarization dependent loss


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% receivers                                                                 receivers
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------


% -------------------------------------------------------------------------
% rx_coherent
% Front-end for polarisation- and phase-diversity coherent receiver
% /src/receiver/
% -------------------------------------------------------------------------
params_coh.pbs_sig.angle = 0;
params_coh.pbs_sig.loss = 0;
params_coh.pbs_sig.extinction_ratio = Inf;
params_coh.pbs_lo.angle = -pi/4;
params_coh.pbs_lo.loss = 0;
params_coh.pbs_lo.extinction_ratio = Inf;
params_coh.pd.responsivity = 1;
params_coh.pd.thermal_noise_density = 0*10e-12;
params_coh.pd.shot_noise = 'off';
params_coh.pd.dark_current = 0;
params_coh.elpf.type = 'bessel';
params_coh.elpf.order = 4;
params_coh.elpf.f3dB = 0.7*symbol_rate;
[sig_x_i,sig_x_q,sig_y_i,sig_y_q] = rx_coherent(sig,sig_lo,params_coh);
% Coherent front-end

% -------------------------------------------------------------------------
% rx_dd
% Receiver incl. optical filtering and direct and interferometric detection
% /src/receiver/
% -------------------------------------------------------------------------
params_rx.type = 'dd';%'id';
params_rx.obpf.type = 'gaussian';
params_rx.obpf.order = 1;
params_rx.obpf.bandwidth = 4*symbol_rate;
params_rx.obpf.centre_frequency = 0;

% Parameters for direct-detection receiver:
params_rx.elpf.type = 'bessel';
params_rx.elpf.order = 4;
params_rx.elpf.f3dB = 0.75*symbol_rate;
params_rx.pd.responsivity = 1;
params_rx.pd.shot_noise = 0;
params_rx.pd.thermal_noise_density = 0;
params_rx.pd.dark_current = 0;

% Parameters for interferometric-detection receiver:
params_rx.mzdi.input_port = 'upper';%'lower';
params_rx.mzdi.mode = 'tuned';%'general';
params_rx.mzdi.delay = 1/symbol_rate;
params_rx.mzdi.phase_shift = 0;

params_rx.upper.pd.responsivity = 1;
params_rx.upper.pd.shot_noise = 0;
params_rx.upper.pd.thermal_noise_density = 0;
params_rx.upper.pd.dark_current = 0;
params_rx.lower.pd.responsivity = 1;
params_rx.lower.pd.shot_noise = 0;
params_rx.lower.pd.thermal_noise_density = 0;
params_rx.lower.pd.dark_current = 0;

params_rx.upper.elpf.type = 'bessel';
params_rx.upper.elpf.order = 4;
params_rx.upper.elpf.f3dB = 0.75*symbol_rate;
params_rx.lower.elpf.type = 'bessel';
params_rx.lower.elpf.order = 4;
params_rx.lower.elpf.f3dB = 0.75*symbol_rate;
params_rx.bd.detection_mode = 'balanced';%'single_ended_upper';%'single_ended_lower';
params_rx.bd.electrical_delay = 0;
params_rx.bd.polarity = 1;%-1.

sig = rx_dd(sig,params_rx);

% For MZDI input_port = 'upper' and mode = 'Tuned'
% 'single_ended_upper' -> compare to original data sequence and ignore 
% 1st bit. The optical signal is AMI.
% 'single_ended_lower' -> compare to inverted original data sequence and 
% ignore 1st bit. The optical signal is DB.



% -------------------------------------------------------------------------
% rx_pin
% PIN receiver with shot and thermal noise and low-pass filtering
% /src/receiver/
% -------------------------------------------------------------------------
params_pin.pd.responsivity = 1;
params_pin.pd.thermal_noise_density = 10e-12;
params_pin.pd.shot_noise = 0;
params_pin.pd.dark_current = 0;
params_pin.elpf.type = 'none';%'bessel';'gaussian';'rc';'rectangular';
params_pin.elpf.order = 4;
params_pin.elpf.f3dB = 0.7*symbol_rate;
sig = rx_pin(sig,params_pin);
% PIN receiver

% -------------------------------------------------------------------------
% rx_resynchronise
% Retiming of an electrical or optical signal based on intensity cross-correlation ("clock recovery")
% /src/receiver/
% -------------------------------------------------------------------------
sig_type = 'opt';%'elec';
[sig,imax] = rx_resynchronise(sig,seq,sig_type);
% Resync signal

% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% sources                                                                   sources
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% opt_broadband_time
% Optical broadband source (time domain generation)
% /src/sources/
% -------------------------------------------------------------------------
power_dbm = 0;  % power of each frequency component, in dBm
sig = opt_broadband_time(power_dbm); 
% Optical broadband source

% -------------------------------------------------------------------------
% opt_comb
% Arbitrary frequency comb generation
% /src/sources/
% -------------------------------------------------------------------------
params_comb.nlines = 13;
% Number of lines in the comb.
params_comb.centre_frequency = reference_frequency;
% Centre frequency of the comb, in Hz.
params_comb.line_spacing = 100e9;
% Frequency spacing of the lines in the comb, in Hz.
params_comb.frequency = params_comb.centre_frequency - (params_comb.nlines - 1)/2*params_comb.line_spacing+(0:1:params_comb.nlines -1 )*params_comb.line_spacing;
% Frequencies of the lines in the comb, in Hz.
params_comb.power = 1e-3*randn(1,params_comb.nlines);
% Power of the lines in the comb, in W.
params_comb.phase = zeros(1,params_comb.nlines);
% Phases of the lines in the comb, in rad.
params_comb.linewidth = 1e6*rand(1,params_comb.nlines);
% Linewidths of the lines in the comb, in Hz.
params_comb.jones.x = ones(1,params_comb.nlines);
params_comb.jones.y = zeros(1,params_comb.nlines);
% Jones vectors of the lines in the comb. Should be normalised.
[sig,params_comb.line_indices] = opt_comb(params_comb);
% Frequency comb

% -------------------------------------------------------------------------
% opt_laser_cw
% Continuous wave (CW) laser
% /src/sources/
% -------------------------------------------------------------------------
params_cw.power = 1.0e-3;
params_cw.linewidth = 0;
params_cw.emission_frequency = 193.1e12;
sig = opt_laser_cw(params_cw);
% CW laser

% -------------------------------------------------------------------------
% opt_laser_dml
% Directly modulated single-mode laser using rate equations
% /src/sources/
% -------------------------------------------------------------------------
params_dml.tau_p = 2.6e-12;          % photon lifetime, in s
params_dml.tau_c = 3.17e-9;          % carrier lifetime, in s
params_dml.n_0 = 2.0e24;             % carrier density at transparency, in 1/m^3
params_dml.sigma_g = 3.34e-20;       % gain cross section, in m^2
params_dml.n_g = 4;                  % group effective index
params_dml.Gamma = 0.2408;           % confinement factor
params_dml.V = 3.6e-17;              % active volume, in m^3
params_dml.epsilon_nl = 2.0e-23;     % gain suppression factor, in m^3
params_dml.alpha = 6;                % linewidth enhancement factor
params_dml.beta = 1.0e-3;            % spontaneous emission factor
params_dml.eta_0 = 0.2;              % differential quantum efficiency
params_dml.emission_frequency = reference_frequency; % emission frequency, in Hz
numparams_dml.ode_solver_options = odeset('RelTol',1e-8);% ODE solver parameters
numparams_dml.npass = 2;             % number of iterations of the ODE solver
numparams_dml.check_convergence = 0;  % display convergence of densities.
sig = opt_laser_dml(sig,params_dml,numparams_dml); 
% Directly-modulated laser

% -------------------------------------------------------------------------
% opt_nosig
% Dummy optical signal generation
% /src/sources/
% -------------------------------------------------------------------------
sig = opt_nosig();
% Dummy optical signal generation

% -------------------------------------------------------------------------
% opt_pulse_gauss
% Gaussian pulse
% /src/sources/
% -------------------------------------------------------------------------
pulse = opt_pulse_gauss(time_array,order,peak_power,position,duration,chirp);
% Gaussian pulse

% -------------------------------------------------------------------------
% opt_pulse_sech
% Hyperbolic secant pulse
% /src/sources/
% -------------------------------------------------------------------------
pulse = opt_pulse_sech(time_array,peak_power,position,duration,chirp);
% Hyperbolic secant pulse

% -------------------------------------------------------------------------
% opt_source_pulse 
% Optical pulse source
% /src/sources/
% -------------------------------------------------------------------------
params_pulse_train.type = 'gaussian';%'sech';
params_pulse_train.order = 1;
params_pulse_train.emission_frequency = reference_frequency + 100e9;
params_pulse_train.peak_power = 1e-3;
params_pulse_train.fwhm = 50e-12;
params_pulse_train.chirp = 5;
sig = opt_source_pulse(seq,params_pulse_train); 
% Optical pulse source


% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% transmitters                                                              transmitters
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% tx
% General optical transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.type = 'ook_nrz';%'ook_rz33', 'ook_rz50', 'ook_rz67',
%                       'ook_db_nrz_delay_add', 'ook_db_nrz_low_pass'
%                       'dpsk_nrz_mzm', 'dpsk_nrz_pm'
%                       'dpsk_rz33_mzm', 'dpsk_rz50_mzm', 'dpsk_rz67_mzm'
%                       'qpsk_nrz_pm_pm', 'qpsk_nrz_mzm_pm', 'qpsk_nrz_iq'
%                       'pam4_nrz_mzm'
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.bit_pattern_1 = bit_pattern_1;
params_tx.bit_pattern_2 = bit_pattern_2;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx(params_tx);
% General optical transmitter

% -------------------------------------------------------------------------
% tx_dpsk_nrz_mzm
% NRZ-DPSK transmitter using a Mach-Zehnder modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_dpsk_nrz_mzm(params_tx); 
% NRZ-DPSK transmitter using a Mach-Zehnder modulator

% -------------------------------------------------------------------------
% tx_dpsk_nrz_pm
% NRZ-DPSK transmitter using a phase modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_dpsk_nrz_pm(params_tx); 
% NRZ-DPSK transmitter using a phase modulator

% -------------------------------------------------------------------------
% tx_dpsk_rz33_mzm
% 33% RZ-DPSK transmitter using a Mach-Zehnder modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_dpsk_rz33_mzm(params_tx); 
% 33% RZ-DPSK transmitter using a Mach-Zehnder modulator

% -------------------------------------------------------------------------
% tx_dpsk_rz50_mzm
% 50% RZ-DPSK transmitter using a Mach-Zehnder modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_dpsk_rx50_mzm(params_tx); 
% 50% RZ-DPSK transmitter using a Mach-Zehnder modulator

% -------------------------------------------------------------------------
% tx_dpsk_rz67_mzm
% 67% RZ-DPSK transmitter using a Mach-Zehnder modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = logical_differential_encoder_binary(bit_pattern);
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_dpsk_rx67_mzm(params_tx); 
% 67% RZ-DPSK transmitter using a Mach-Zehnder modulator

% -------------------------------------------------------------------------
% tx_duobinary_nrz_delay_add
% Delay-and-add NRZ duobinary transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_duobinary_nrz_delay_add(params_tx);
% Delay-and-add NRZ duobinary transmitter

% -------------------------------------------------------------------------
% tx_duobinary_nrz_low_pass
% Low-pass filter NRZ duobinary transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_duobinary_nrz_low_pass(params_tx); 
% Low-pass filter NRZ duobinary transmitter

% -------------------------------------------------------------------------
% tx_fsk
% Ideal optical binary FSK transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.type = 'fsk';%'cpfsk';'msk';'gmsk';
params_tx.symbol_rate = symbol_rate;
params_tx.emission_frequency = reference_frequency;
params_tx.tone_spacing = 2*params_tx.symbol_rate; % For FSK and CPFSK only.
params_tx.rise_time = 1/params_tx.symbol_rate/4; 
params_tx.bt = 0.3;                               % For GMSK only.
params_tx.power = 1.0e-3;
sig = tx_fsk(bit_pattern,params_tx);
% Ideal optical binary FSK transmitter

% -------------------------------------------------------------------------
% tx_laser_chirped
% Black-box laser with frequency chirp
% /src/transmitters/
% -------------------------------------------------------------------------
params_chirped_tx.emission_frequency = reference_frequency;   % laser emission frequency, in Hz
params_chirped_tx.power_average = 1.0e-3;                     % signal average power, in W
params_chirped_tx.extinction_ratio = 10;                      % signal extinction ratio, in dB
params_chirped_tx.alpha = 2.5;                                % alpha parameter
params_chirped_tx.kappa = 12e12;                              % adiabatic chirp, in Hz/W
[sig,freq_data] = tx_laser_chirped(nrz_data_sig,params_chirped_tx); 
% Black-box laser with frequency chirp

% -------------------------------------------------------------------------
% tx_ook_nrz
% NRZ-OOK transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_ook_nrz(params_tx);
% NRZ-OOK transmitter

% -------------------------------------------------------------------------
% tx_ook_rz33
% 33% RZ-OOK transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_ook_rz33(params_tx); 
% 33% RZ-OOK transmitter

% -------------------------------------------------------------------------
% tx_ook_rz50
% 50% RZ-OOK transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_ook_rz50(params_tx); 
% 50% RZ-OOK transmitter

% -------------------------------------------------------------------------
% tx_ook_rz67
% 67% RZ-OOK transmitter
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_ook_rz67(params_tx); 
% 67% RZ-OOK transmitter

% -------------------------------------------------------------------------
% tx_pam4_mzm
% PAM4 transmitter using a Mach-Zehnder modulator (chirp-free, infinite ER)
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.bit_pattern_1 = bit_pattern_1;
params_tx.bit_pattern_2 = bit_pattern_2;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_pam4_mzm(params); 
% PAM4 transmitter

% -------------------------------------------------------------------------
% tx_qam16_nrz_iq
% NRZ-16QAM transmitter using ideal optical IQ modulator
% /src/transmitters/
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1e-3;
params_tx.bit_pattern.b1 = round(rand(1,nsymbols));
params_tx.bit_pattern.b2 = round(rand(1,nsymbols));
params_tx.bit_pattern.b3 = round(rand(1,nsymbols));
params_tx.bit_pattern.b4 = round(rand(1,nsymbols));
params_tx.rise_time = 1/symbol_rate/4;
params_tx.alpha = 0.64;
params_tx.beta = 1;
[sig,ui,uq] = tx_qam16_nrz_iq(params_tx); 
% NRZ-16QAM transmitter using ideal optical IQ modulator

% -------------------------------------------------------------------------
% tx_qpsk_nrz_iq
% NRZ-QPSK transmitter using optical IQ modulator
% /src/transmitters/
% -------------------------------------------------------------------------
u = round(rand(1,nsymbols));
v = round(rand(1,nsymbols));
[params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('parallel',u,v);
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_qpsk_nrz_iq(params_tx); 
% NRZ-QPSK transmitter using optical IQ modulator

% -------------------------------------------------------------------------
% tx_qpsk_nrz_mzm_pm
% NRZ-QPSK transmitter using a Mach-Zehnder modulator followed by a phase modulator
% /src/transmitters/
% -------------------------------------------------------------------------
u = round(rand(1,nsymbols));
v = round(rand(1,nsymbols));
[params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('serial',u,v);
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_qpsk_nrz_mzm_pm(params_tx);
% NRZ-QPSK transmitter using a Mach-Zehnder modulator followed by a phase modulator

% -------------------------------------------------------------------------
% tx_qpsk_nrz_pm_pm
% NRZ-QPSK transmitter using two phase modulators
% /src/transmitters/
% -------------------------------------------------------------------------
u = round(rand(1,nsymbols));
v = round(rand(1,nsymbols));
[params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('serial',u,v);
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_qpsk_nrz_pm_pm(params_tx); 
% NRZ-QPSK transmitter using two phase modulators

% -------------------------------------------------------------------------
% tx_qpsk_rz50_iq
% 50% RZ-DQPSK transmitter using optical IQ modulator
% /src/transmitters/
% -------------------------------------------------------------------------
u = round(rand(1,nsymbols));
v = round(rand(1,nsymbols));
[params_tx.bit_pattern_1,params_tx.bit_pattern_2,~,~,~] = logical_differential_encoder_dqpsk('parallel',u,v);
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.rise_time = 1/symbol_rate/4;
sig = tx_qpsk_rz50_iq(params_tx);
% 50% RZ-DQPSK transmitter using optical IQ modulator
