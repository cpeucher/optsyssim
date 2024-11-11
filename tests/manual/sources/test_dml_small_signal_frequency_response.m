% -------------------------------------------------------------------------
% We compare the small-signal frequency response of a directly modulated
% laser obtained numerically through resolution of the rate equations to
% the theoretical response obtained analytically by linearisation of the
% rate equations.
%
% To this end we apply sinusoidal signals of varying frequencies with a 
% small modulation index to the rate equation model of the directly
% modulated laser and characterise the resulting peak-to-peak optical
% power.
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
do_debug = 0;
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
nsymbols = 32;
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
% DML parameters
% -------------------------------------------------------------------------
params_dml.tau_p = 2.6e-12;          % photon lifetime, in s
params_dml.tau_c = 3.17e-9;          % carrier lifetime, in s
params_dml.n_0 = 2.0e24;             % carrier density at transparency, in 1/m^3
params_dml.sigma_g = 3.34e-20;       % gain cross section, in m^2
params_dml.n_g = 4;                  % group effective index
params_dml.Gamma = 0.2408;           % confinement factor
params_dml.V = 3.6e-17;              % active volume, in m^3
params_dml.epsilon_nl = 0*2.0e-23;   % gain suppression factor, in m^3
params_dml.alpha = 6;                % linewidth enhancement factor
params_dml.beta = 1.0e-3;            % spontaneous emission factor
params_dml.eta_0 = 0.2;              % differential quantum efficiency
params_dml.emission_frequency = reference_frequency; % emission frequency, in Hz

% Numerical parameters:
numparams_dml.ode_solver_options = odeset('RelTol',1e-8);% ODE solver parameters
numparams_dml.npass = 2;             % number of iterations of the ODE solver
numparams_dml.check_convergence = 1;  % display convergence of densities.


% -------------------------------------------------------------------------
% Electrical driving signal
% -------------------------------------------------------------------------
ibias = 50e-3;
% Bias current.

ipp = 1.0e-3;
% Peak-to-peak current.





        
        
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Determination of rf frequencies for small signal frequency response
% characterisation
% -------------------------------------------------------------------------
nperiods_min = 4;
% Minimum number of signal periods over time window. This will decide the
% lowest possible modulation frequency.
nsamples_per_period_max = nsamples/nperiods_min;
% Maximimum number of samples per period.
mod_period_max = nsamples*dt/nperiods_min;
% Maximum allowable period of the modulating signal.
mod_freq_min = 1/mod_period_max;
% Corresponding minimum frequency.

fprintf(1,'\n%s\t%3.2f\t%s','Minimum modulation frequency:',mod_freq_min/1.0e9,'GHz');

nsamples_per_period_min = 32;
% Minimum number of samples per period of the modulating signal.
mod_period_min = nsamples_per_period_min*dt;
% Corresponding minimum period of the modulating signal.
mod_freq_max = 1/mod_period_min;
% Corresponding maximum frequency of the modulating signal.
nperiods_max = nsamples*dt/mod_period_min;
% Maximum number of periods over the the time window. 

fprintf(1,'\n%s\t%3.2f\t%s\n\n','Maximum modulation frequency:',mod_freq_max/1.0e9,'GHz');


mod_freq_numerical = [nperiods_min:1:nperiods_max]/dt/nsamples;
% Modulation frequencies, in Hz.

% index     frequency
% 5         2.5 GHz
% 13        5 GHz
% 21        7.5 GHz
% 29        10 GHz
% 45        15 GHz
% 61        20 GHz
% 77        25 GHz


if do_debug
    mod_freq_numerical = mod_freq_numerical(29);
    % Single frequency for debug purpose (waveform visualisation).
else
    mod_freq_numerical = mod_freq_numerical(mod_freq_numerical < 20e9);
    % Restrict the modulation frequencies to 20 GHz.    
end










%%
% -------------------------------------------------------------------------
% Numerical determination of small-signal frequency response
% -------------------------------------------------------------------------
delta_power = zeros(1,length(mod_freq_numerical));
% Pre-initialisation.

for ifreq = 1:5:length(mod_freq_numerical) 
    % Loop over modulation frequencies.
    
    sig_drive = ibias + 0.5*ipp*cos(2*pi*mod_freq_numerical(ifreq)*time_array);
    % Sinusoidal driving signal applied to the DML.
    
    if do_debug
        % Plot driving signal.        
        figure('Name','Driving signal')
        plot(time_array/1.0e-12,sig_drive/1.0e-3);
        xlabel('time (ps)');
        ylabel('current (mA)');        
    end
    
    
    sig = opt_laser_dml(sig_drive,params_dml,numparams_dml); 
    % DML.
    
    sig_power = abs(sig.x).^2;   
    % Emitted power.    
    
    if do_debug  
        % Plot emitted power
        figure('Name','Emitted power')
        plot(time_array/1.0e-12,sig_power/1.0e-3)
        xlabel('time (ps)')
        ylabel('emitted power (mW)')        
    end
    
    delta_power(ifreq) = max(sig_power) - min(sig_power);
    % Peak-to-peak optical power.
        
end
% End of loop over modulation frequencies.



%%
% -------------------------------------------------------------------------
% Analytical small-signal frequency response
% -------------------------------------------------------------------------
mod_freq_analytical = logspace(0,11,1000);
% Modulation frequency range.
dml_response = calc_dml_frequency_response(ibias,params_dml,mod_freq_analytical);
% Analytical calculation of the small signal frequency response.

%%
% -------------------------------------------------------------------------
% Compare analytical and numerical small-signal frequency responses
% -------------------------------------------------------------------------
fig_name = [file_name_core_figure '_'];
hfig = figure('Name',fig_name);
plot(mod_freq_analytical/1.0e9,dml_response.HdB,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(mod_freq_numerical/1.0e9,10*log10(delta_power/delta_power(1)),'Color','r','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('frequency (GHz)','Interpreter',fig.interpreter)
ylabel('$|H|$ (dB)','Interpreter',fig.interpreter)
legend('analytical','numerical','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([0 20])
% ylim([-0.5 0.05])
if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end   





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
