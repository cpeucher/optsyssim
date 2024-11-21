% -------------------------------------------------------------------------
% Test of coherent front end function rx_coherent
%
% 2021-11-27
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
% results
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));

%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.interpreter = 'latex';
fig.font_size = 18;

mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));


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
reference_frequency = 193.1e12;
nsamples_per_symbol = 32;
nsymbols = 128;
symbol_rate = 10e9;

nsamples = nsamples_per_symbol*nsymbols;
sample_rate = nsamples_per_symbol*symbol_rate;
[time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);        
        
  
% reference_frequency = 193.1e12;
% df = 10e6;
% nsamples = 2^14;
% 
% sample_rate = nsamples*df;
% dt = 1/sample_rate;
% time_array = (0:nsamples-1)*dt;
% frequency_array = (-nsamples/2:nsamples/2-1)*df;

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
% Load essential physical constants

% ------------------------------------------------------------------------- 
% Start time 
% ------------------------------------------------------------------------- 
start_time = datetime("now");
fprintf('\n\n%s%s\n\n','Simulation started on ',start_time);





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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% PRBS
% -------------------------------------------------------------------------
params_prbs_1.type = 'shift_register';
params_prbs_1.order = 7;
params_prbs_1.poly = [7 6 0];%[7 3 0];[7 1 0];

params_prbs_2.type = 'shift_register';
params_prbs_2.order = 7;
params_prbs_2.poly = [7 3 0];%[7 1 0];[7 6 0];


% -------------------------------------------------------------------------
% Electrical signal 
% -------------------------------------------------------------------------
rise_time = 1/symbol_rate/4;  % NRZ signal rise time, in s.

% -------------------------------------------------------------------------
% QPSK transmitter
% -------------------------------------------------------------------------
params_tx.emission_frequency = reference_frequency;
params_tx.linewidth = 0;
params_tx.power = 1.0e-3;
params_tx.rise_time = rise_time;

% -------------------------------------------------------------------------
% Waveform display
% -------------------------------------------------------------------------
nsymbols_display = 8;





        
        
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------

params_prbs_1.seed = [1 1 1 0 1 1 0];
bit_pattern_1 = logical_prbs(params_prbs_1);
bit_pattern_1 = logical_adapt_binary_sequence(bit_pattern_1,nsymbols);

params_prbs_1.seed = [0 0 0 1 1 1 0];
bit_pattern_3 = logical_prbs(params_prbs_1);
bit_pattern_3 = logical_adapt_binary_sequence(bit_pattern_3,nsymbols);

params_prbs_2.seed = [0 1 0 0 1 0 0];
bit_pattern_2 = logical_prbs(params_prbs_2);
bit_pattern_2 = logical_adapt_binary_sequence(bit_pattern_2,nsymbols);

params_prbs_2.seed = [1 0 1 0 1 1 0];
bit_pattern_4 = logical_prbs(params_prbs_2);
bit_pattern_4 = logical_adapt_binary_sequence(bit_pattern_4,nsymbols);
% Generate random sequences and adapt lengths to simulation time window

nrz_data_sig_1 = elec_pulse_sequence_nrz(bit_pattern_1,rise_time);
nrz_data_sig_2 = elec_pulse_sequence_nrz(bit_pattern_2,rise_time);
nrz_data_sig_3 = elec_pulse_sequence_nrz(bit_pattern_3,rise_time);
nrz_data_sig_4 = elec_pulse_sequence_nrz(bit_pattern_4,rise_time);
% Encode the binary bit pattern onto NRZ signals


figure('Name','Original bit pattern')
subplot(4,1,1)
plot(time_array/1.0e-12,nrz_data_sig_1,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
ylim([0 1.1]);
title('bit pattern 1')
subplot(4,1,2)
plot(time_array/1.0e-12,nrz_data_sig_2,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
ylim([0 1.1]);
title('bit pattern 2')
subplot(4,1,3)
plot(time_array/1.0e-12,nrz_data_sig_3,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
ylim([0 1.1]);
title('bit pattern 3')
subplot(4,1,4)
plot(time_array/1.0e-12,nrz_data_sig_4,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
ylim([0 1.1]);
title('bit pattern 4')
% Check original patterns



params_tx.bit_pattern_1 = bit_pattern_2;
params_tx.bit_pattern_2 = bit_pattern_1;
% Patterns to apply to the IQ modulator
% Observe that we invert the patterns to indeed recover bit_pattern_1 on
% the I output and bit_pattern_2 on the Q output of the coherent front-end.
sig_x = tx_qpsk_nrz_iq(params_tx);
% QPSK transmitter for -x polarisation

params_tx.bit_pattern_1 = bit_pattern_4;
params_tx.bit_pattern_2 = bit_pattern_3;
sig_y = tx_qpsk_nrz_iq(params_tx);
% QPSK transmitter for -y polarisation
% The signal is, for the time being, polarised along -x.

sig_sig.x = sig_x.x;
sig_sig.y = sig_y.x; 
% Combine the two polarisations

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
params_eye.name = 'Eye diagram - transmitter output';
meas_eye(sig_x,params_eye);
% Eye diagram


params_constellation.type = 'scatter_transitions';%'transitions';%'scatter_transitions';
params_constellation.pol = 'x';%'y','both';
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
char_opt_constellation(sig_x,params_constellation);
% Constellation visualisation

params_cw.power = 1.0e-3;
params_cw.linewidth = 0;
params_cw.emission_frequency = params_tx.emission_frequency;
sig_lo = opt_laser_cw(params_cw);
% Local oscillator CW laser

params_coherent_front_end.pbs_sig.angle = 0;
params_coherent_front_end.pbs_sig.loss = 0;
params_coherent_front_end.pbs_sig.extinction_ratio = Inf;
params_coherent_front_end.pbs_lo.angle = -pi/4;
params_coherent_front_end.pbs_lo.loss = 0;
params_coherent_front_end.pbs_lo.extinction_ratio = Inf;
params_coherent_front_end.pd.responsivity = 1;
params_coherent_front_end.pd.thermal_noise_density = 0*10e-12;
params_coherent_front_end.pd.shot_noise = 0;
params_coherent_front_end.pd.dark_current = 0;
params_coherent_front_end.elpf.type = 'none';%'bessel';%'none','bessel','gaussian','rc', 'rectangular';
params_coherent_front_end.elpf.order = 4;
params_coherent_front_end.elpf.f3dB = 0.7*symbol_rate;
[sig_x_i,sig_x_q,sig_y_i,sig_y_q] = rx_coherent(sig_sig,sig_lo,params_coherent_front_end);



figure('Name','Recovered bit pattern')
subplot(4,1,1)
plot(time_array/1.0e-12,sig_x_i,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
% ylim([0 1.1]);
title('in-phase (-x pol)')
subplot(4,1,2)
plot(time_array/1.0e-12,sig_x_q,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
% ylim([0 1.1]);
title('quadrature (-x pol)')
subplot(4,1,3)
plot(time_array/1.0e-12,sig_y_i,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
% ylim([0 1.1]);
title('in-phase (-y pol)')
subplot(4,1,4)
plot(time_array/1.0e-12,sig_y_q,'b')
xlabel('time (ps)')
ylabel('voltage (a.u.)')
xlim([0 nsymbols_display/symbol_rate/1.0e-12]);
% ylim([0 1.1]);
title('quadrature (-y pol)')
% Check received bit patterns









% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(2)

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
% Display simulation duration
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
core_display_duration(start_time,datetime("now"));
