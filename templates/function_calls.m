




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
% Plot optical constellation.

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
% mes_esa
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
% dsp                                                                       dsp
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------
% dsp/general                                                               dsp/general
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

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
% dsp/digicoms                                                              dsp/digicoms
% --------------------------------------------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------------------------------------------

% -------------------------------------------------------------------------
% dsp_constellation_define
% Define constellations of digital modulation formats
% /src/dsp/digicoms/
% -------------------------------------------------------------------------
[constellation,norm_es,norm_emax] = dsp_constellation_define(type,m);
% Define constellation

% -------------------------------------------------------------------------
% dsp_constellation_plot
% Quickly plot digital constellation diagram
% /src/dsp/digicoms/
% -------------------------------------------------------------------------
constellation_type = 'plain';%'heat';'cluster';
constellation_name = 'received constellation';
dsp_constellation_plot(symbs,constellation_type,constellation_name,[-5:1:5]);
% Plot constellation

% -------------------------------------------------------------------------
% dsp_conv_bin2dec
% Converts vector of bits (binary) to vector of words (decimal)
% /src/dsp/digicoms/
% -------------------------------------------------------------------------
[words_dec,words_bin] = dsp_conv_bin2dec(bits_bin,log2(m)); 
% Convert bit string to words

% -------------------------------------------------------------------------
% dsp_data_binary
% Generation of uniformly distributed binary random data
% /src/dsp/digicoms/
% -------------------------------------------------------------------------
bits_bin = dsp_data_binary(nsymbols,m);
% Generation of uniformly distributed binary random data

% -------------------------------------------------------------------------
% dsp_mapping
% Mapping decimal words to complex symbols
% /src/dsp/digicoms/
% -------------------------------------------------------------------------
symbs = dsp_mapping(words_dec,constellation);
% Mapping


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
% /src/auxiliary/
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
params_tx.laser_power = 1.0e-3;
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
