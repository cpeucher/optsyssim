




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% amplifiers
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% characterization
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% char_opt_average_power 
% Optical average power meter (PWM)
% /src/characterization/
% -------------------------------------------------------------------------
average_power = char_opt_average_power(sig);
% Calculate signal average power


% -------------------------------------------------------------------------
% char_opt_temporal_chirp
% Temporal chirp of an optical signal
% /src/characterization/
% -------------------------------------------------------------------------
chirp = char_opt_temporal_chirp(sig);
% Calculate temporal chirp





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% core
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% electrical
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% elec_coder_nrz
% Electrical NRZ encoder without rise time
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_coder_nrz(bit_pattern,nsamples_per_symbol); 
% NRZ encoder without rise time

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
% elec_pulse_sequence_nrz
% Electrical NRZ pulse sequence generation
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% NRZ sequence generation

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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% fibres
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


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
params_fibre.loss_alpha = loss_conv_lin_log(params_fibre.loss);
sig = opt_nlse_scalar_basic(sig,params_fibre,numparams_fibre);
% Basic scalar nonlinear Schroedinger equation





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% modulators
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% numerical
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% num_diff_1d_pb
% Differentiate evenly spaced 1D data with periodic boundary conditions
% /src/numerical/differentiation/
% -------------------------------------------------------------------------
dfunc = num_diff_1d_pb(func,h)
% Differentiation

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




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% passive
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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
 
 
 
 


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% polarization
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% receivers
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% sources
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% transmitters
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

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
