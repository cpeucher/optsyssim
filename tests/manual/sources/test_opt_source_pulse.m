% -------------------------------------------------------------------------
% Test of opt_source_pulse function
%
% We generate different types of optical pulse trains using the
% opt_source_pulse function
% 1. Single pulse
% 2. Periodic pulse train
% 3. Alternate pulse train
% 4. Intensity binary-modulated pulse train
% 5. 16-QAM modulated pulse train
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
nsamples_per_symbol = 128;
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
% Test 1: generation of isolated optical pulse
% -------------------------------------------------------------------------
symbs = zeros(1,nsymbols);
symbs(nsymbols/2) = 1;


params_pulse_train.pulse_shape = 'gaussian';%'sech';
params_pulse_train.order = 1;
params_pulse_train.emission_frequency = reference_frequency + 100e9;
params_pulse_train.peak_power = 1e-3;
params_pulse_train.fwhm = 100e-12;
params_pulse_train.chirp = 0;
sig = opt_source_pulse(symbs,params_pulse_train); 


params_scope.visualisers = {'power'};%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'scope - isolated pulse';
meas_scope(sig.x,params_scope);
% Oscilloscope

params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'osa - isolated pulse';
meas_osa(sig,params_osa);
% Optical spectrum analyser.



% -------------------------------------------------------------------------
% Test 2: generation of unmodulated pulse train.
% -------------------------------------------------------------------------
symbs = ones(1,nsymbols);

params_pulse_train.pulse_shape = 'sech';
params_pulse_train.emission_frequency = reference_frequency;
params_pulse_train.fwhm = 10e-12;

sig = opt_source_pulse(symbs,params_pulse_train); 

params_scope.visualisers = {'power'};%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.ascii = 0;
params_scope.name = 'scope - unmodulated pulse train';
meas_scope(sig.x,params_scope);
% Oscilloscope


params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 0;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'osa  - unmodulated pulse train';
meas_osa(sig,params_osa);
% Optical spectrum analyser.


% -------------------------------------------------------------------------
% Test 3: generation of alternate pulse train.
% -------------------------------------------------------------------------
symbs = logical_alternate('alternate0',nsymbols);

params_pulse_train.pulse_shape = 'sech';
params_pulse_train.emission_frequency = reference_frequency - 300e9;
params_pulse_train.fwhm = 10e-12;

sig = opt_source_pulse(symbs,params_pulse_train); 

params_scope.visualisers = {'power'};%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'scope - alternate pulse train';
meas_scope(sig.x,params_scope);
% Oscilloscope


params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 0;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'osa  - alternate pulse train';
meas_osa(sig,params_osa);
% Optical spectrum analyser.




% -------------------------------------------------------------------------
% Test 4: generation of binary intensity modulated pulse train.
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
symbs = logical_prbs(params_prbs);
symbs = logical_adapt_binary_sequence(symbs,nsymbols);

params_pulse_train.pulse_shape = 'sech';
params_pulse_train.emission_frequency = reference_frequency - 300e9;
params_pulse_train.fwhm = 10e-12;

sig = opt_source_pulse(symbs,params_pulse_train); 

params_scope.visualisers = {'power'};%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'scope - binary intensity modulated pulse train';
meas_scope(sig.x,params_scope);
% Oscilloscope


params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 0;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'osa  - binary intensity modulated pulse train';
meas_osa(sig,params_osa);
% Optical spectrum analyser.


% -------------------------------------------------------------------------
% Test 5: generation of 16-QAM modulated pulse train.
% -------------------------------------------------------------------------
m = 16;
% Constellation order.
[constellation,norm_es,norm_emax] = define_constellation('qam16_gray',m);
% Define the constellation.
data = generate_binary(nsymbols,m);
% Generate binary data.
[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 4-bit words.
symbs = mapping(words_dec,constellation);
% Mapping.

constellation_type = 'plain';%'heat';'cluster';
constellation_name = 'generated digital constellation';
plot_constellation(symbs,constellation_type,constellation_name,[-5:1:5]);


params_pulse_train.pulse_shape = 'sech';
params_pulse_train.emission_frequency = reference_frequency - 0*300e9;
params_pulse_train.fwhm = 10e-12;

sig = opt_source_pulse(symbs,params_pulse_train); 

params_scope.visualisers = {'power'};%'power_phase';'power_chirp';'power_phase_chirp';
params_scope.display_interval = [0 time_array(end)];
params_scope.save.ascii = 0;
params_scope.name = 'scope - binary intensity modulated pulse train';
meas_scope(sig.x,params_scope);
% Oscilloscope


params_osa.pol = 'x';%'y','both';
params_osa.display_interval = [frequency_array(1) frequency_array(end)];
params_osa.resolution_bandwidth = 0;
params_osa.sensitivity = -100;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'osa  - binary intensity modulated pulse train';
meas_osa(sig,params_osa);
% Optical spectrum analyser.


samps = sqrt(1000)*sig.x(nsamples_per_symbol/2:nsamples_per_symbol:end);
% Sample the generated waveform.

constellation_name = 'recovered digital constellation';
plot_constellation(samps,constellation_type,constellation_name,[-5:1:5]);

symbs_error = sum(abs(samps - symbs).^2);        






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
