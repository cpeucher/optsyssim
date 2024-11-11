% -------------------------------------------------------------------------
% We calculate the static emitted power versus applied current curve of a
% semiconductor laser and compare it to the theoretical curve.
% 
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
nsamples_per_symbol = 32;
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
ibias_min = 0e-3;       % minimum bias current, in A
ibias_max = 60e-3;      % maximum bias current, in A
ibias_step = 0.5e-3;    % bias current step, in A     




        
        
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------

if do_debug
    ibias = 25e-3;
    numparams_dml.check_convergence = 1;
else
    ibias = [ibias_min:ibias_step:ibias_max];
    numparams_dml.check_convergence = 0;
end

emitted_power_numerical = zeros(1,length(ibias));
emitted_power_analytical = zeros(1,length(ibias));
% Preinitialization.

% -------------------------------------------------------------------------
% Numerical determination of emitted power versus bias current curve
% -------------------------------------------------------------------------
for ii = 1:length(ibias)
    % Loop over bias current.
    
    sig_drive = ibias(ii)*ones(1,nsamples);
    % Constant driving signal.
    
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
    
    emitted_power_numerical(ii) = mean(sig_power);
    % Numerical emitted power.    
    
    dml_response = calc_dml_frequency_response(ibias(ii),params_dml,logspace(0,11,1000));
    % Analytical calculation of laser response.
    
    emitted_power_analytical(ii) = dml_response.Pout;
    % Analytical emitted power.
    
end
% End of loop over bias current.
%%
% -------------------------------------------------------------------------
% Analytical  response
% -------------------------------------------------------------------------



%%
% -------------------------------------------------------------------------
% Compare analytical and numerical power versus current curves
% -------------------------------------------------------------------------
fig_name = [file_name_core_figure '_'];
hfig = figure('Name',fig_name);
plot(ibias/1.0e-3,emitted_power_analytical/1.0e-3,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(ibias/1.0e-3,emitted_power_numerical/1.0e-3,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('bias current (mA)','Interpreter',fig.interpreter)
ylabel('emitted power (mW)','Interpreter',fig.interpreter)
legend('analytical','numerical','Location','NorthWest','Box','off','Interpreter',fig.interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([0 60])
ylim([0 5])
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
