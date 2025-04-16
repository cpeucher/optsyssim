% -------------------------------------------------------------------------
% Test of 16-QAM signal generation using an IQ modulator driven with 2
% PAM-4 electrical signals (i.e. this "historical" implementation does not
% use a DAC). 
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
nsamples_per_symbol = 64;
nsymbols = 2048;
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
% Generation of NRZ QAM16 
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
% Check bit patterns
% -------------------------------------------------------------------------
rise_time = 1/symbol_rate/4;  
% NRZ signal rise time, in s

nrz_data_sig_1 = elec_pulse_sequence_nrz(params_tx.bit_pattern.b1,rise_time);
nrz_data_sig_2 = elec_pulse_sequence_nrz(params_tx.bit_pattern.b2,rise_time);
nrz_data_sig_3 = elec_pulse_sequence_nrz(params_tx.bit_pattern.b3,rise_time);
nrz_data_sig_4 = elec_pulse_sequence_nrz(params_tx.bit_pattern.b4,rise_time);

figure('Name','Original bit patterns')
subplot(4,1,1)
plot(time_array*symbol_rate,nrz_data_sig_1)
hold on
stem([0:1:nsymbols - 1] + 0.5,nrz_data_sig_1(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('pattern 1')
subplot(4,1,2)
plot(time_array*symbol_rate,nrz_data_sig_2)
hold on
stem([0:1:nsymbols - 1] + 0.5,nrz_data_sig_2(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('pattern 2')
subplot(4,1,3)
plot(time_array*symbol_rate,nrz_data_sig_3)
hold on
stem([0:1:nsymbols - 1] + 0.5,nrz_data_sig_3(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('pattern 3')
subplot(4,1,4)
plot(time_array*symbol_rate,nrz_data_sig_4)
hold on
stem([0:1:nsymbols - 1] + 0.5,nrz_data_sig_4(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('pattern 4')
xlabel('bit number')

% -------------------------------------------------------------------------
% Check IQ modulator driving signals
% -------------------------------------------------------------------------
figure('Name','IQ modulator Driving signals')
subplot(2,1,1)
plot(time_array*symbol_rate,ui)
hold on
stem([0:1:nsymbols - 1] + 0.5,ui(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('u_i')
subplot(2,1,2)
plot(time_array*symbol_rate,uq)
hold on
stem([0:1:nsymbols - 1] + 0.5,uq(nsamples_per_symbol/2:nsamples_per_symbol:end));
xlim([0 16])
xticks([0.5:1:16])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'})
ylabel('u_q')
xlabel('symbol number')



% -------------------------------------------------------------------------
% Optical spectrum
% -------------------------------------------------------------------------
params_osa.pol = 'x';%'y';'both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 1e9;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'Optical spectrum - QAM16 NRZ';
meas_osa(sig,params_osa); 
% Optical spectrum analyser.


% -------------------------------------------------------------------------
% Eye diagram of driving signal
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
params_eye.name = 'Eye diagram - uI';
meas_eye(ui,params_eye);
% Plot eye diagram

% -------------------------------------------------------------------------
% Intensity eye diagram
% -------------------------------------------------------------------------
params_eye.name = 'Eye diagram - QAM16 NRZ';
meas_eye(abs(sig.x).^2,params_eye);

% -------------------------------------------------------------------------
% Constellation
% -------------------------------------------------------------------------
params_constellation.type = 'scatter_transitions';%'transitions';%'scatter';
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
params_constellation.plot_circular_grid = 0;
params_constellation.plot_axes = 0;
params_constellation.name = 'Constellation - QAM16 NRZ';
char_opt_constellation(sig,params_constellation); 
% Plot optical constellation






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
