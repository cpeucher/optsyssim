
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% auxiliary
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% prod_mm
% Product of time or frequency dependent 2x2 matrices
% /src/auxiliary/
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
% elec_pulse_sequence_nrz
% Electrical NRZ pulse sequence generation
% /src/electrical/
% -------------------------------------------------------------------------
sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% NRZ sequence generation





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






