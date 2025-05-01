% -------------------------------------------------------------------------
% Test of eye opening calculation using eops and eopw functions
%
% A signal severely affected by inter-symbole interference (but where the
% eye remains open) is generated using the opt_laser_dml function.
% Its eye opening is then characterized using:
% 1. The single sampling time eye opening determination function eops
% 2. The wide sampling time eye opening determination function eopw, where
% the sampling windows is set to a single sample
% 3. The wide sampling time eye opening determination function eopw, where
% the sampling windows is set to a range of samples
%
% 1. and 2. results in the same eye opening values. The returned values of
% the optimum sampling time are different due to the fact that retiming is 
% performed in 2, but not in 1. It is however shown that the
% sampling instant is the same with respect to the eye diagram.
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


params_prbs.type = 'shift_register';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window. 


rise_time = 1/symbol_rate/4;  % NRZ signal rise time, in s.
nrz_data_sig = elec_pulse_sequence_nrz(bit_pattern,rise_time);
% Encode the binary bit pattern onto an NRZ signal.


idc = 30.0e-3;
% DC level of the DML driving signal, in A
ipp = 20.0e-3;
% Peak-to-peak value of the DML driving signal, in A.

sig_drive = idc + ipp*(nrz_data_sig - 0.5);
% Driver + bias tee


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
sig = opt_laser_dml(sig_drive,params_dml,numparams_dml); 
% DML


params_osa.pol = 'x';%'y';'both';
params_osa.display_interval = [-6*symbol_rate 6*symbol_rate];
params_osa.resolution_bandwidth = 12.5e9;
params_osa.sensitivity = -60;
params_osa.display = 1;
params_osa.save.txt = 0;
params_osa.save.emf = 0;
params_osa.save.jpg = 0;
params_osa.name = 'DML output spectrum';
meas_osa(sig,params_osa); 
% Optical spectrum analyser

params_scope.visualisers = {'power'};
params_scope.display_interval = [0 time_array(end)];
params_scope.save.emf = 0;
params_scope.name = 'DML output waveform';
meas_scope(sig.x,params_scope); 
% Scope

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
params_eye.name = 'DML output eye diagram';
meas_eye(sig,params_eye);
% Plot eye diagram



% -------------------------------------------------------------------------
% Extraction of the eye opening at the optimum sampling time using legacy
% function eops
% -------------------------------------------------------------------------


% params_eop.norm_power = 2*char_opt_average_power(sig);
params_eop.norm_power = 1;
% Normalisation power set to 1 to return the absolute eye opening

% First with the "single-sample" eye opening function eops
fprintf('\n\n\n%s','Use of eops function')
fprintf('\n%s',    '====================')


tic
[eop,eo_norm,opt_sampling] = eops(abs(sig.x).^2,bit_pattern,params_eop.norm_power,-1)
toc

vline((opt_sampling-1)*dt/1.0e-12,'r-')
% Indicate the optimum sampling time on the eye diagram.


% -------------------------------------------------------------------------
% Extraction of the eye opening at the optimum sampling time using eopw
% restricted to 1 sample
% -------------------------------------------------------------------------

% Next with the "wide" eye opening function eopw with 1 sample width
params_eop.method = 'eow';%'cr';
params_eop.bit_pattern = bit_pattern;
params_eop.eo_b2b = -1;%eo_norm;
params_eop.display_eye = 1;
params_eop.eye_width_samples = 1;
params_eop.bits_to_ignore_start = 0;
params_eop.eye_display_name = 'Eye diagram with eye opening';

fprintf('\n\n\n%s','Use of eopw function with 1 sample baseline')
fprintf('\n%s',    '===========================================')


tic
eopw_res = eopw(abs(sig.x).^2,params_eop)
toc


% Conclusion:
% Same eye opening. Different optimum sampling indices are returned. But
% this is linked to the fact that the signal is retimed in eopw, but not in
% eops
% The location of the optimum is clearly the same on the eye diagram,
% We furthermore have checked that the difference between the two returned
% values of optimum sample index correspond to the index shift in the
% resynchronise function used in eopw.



%%
% -------------------------------------------------------------------------
% Finally we look for the eye opening defined by considering the maximum
% height of a rectangle of a specified width that fits withing the eye
% opening.
% -------------------------------------------------------------------------

params_eop.eye_width_samples = 18;
params_eop.bits_to_ignore_start = 0;
params_eop.eye_display_name = 'Eye diagram with eye opening';

fprintf('\n\n\n%s%i%s','Use of eopw function with ',params_eop.eye_width_samples,' sample(s) baseline')
fprintf('\n%s','==============================================')

tic
eopw_res = eopw(abs(sig.x).^2,params_eop)
toc








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
