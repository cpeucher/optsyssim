% -------------------------------------------------------------------------
% Test of pol_pbs function
% We input linearly polarized light to a polarization beam splitter (PBS).
% The angle of the PBS is rotated. 
% We look at the evolution of the power at the output of each port.
% The insertion loss (IL) and polarization extinction ratio (PER) can be
% varied.
% Two modes of operation are tested:
% 1. A static mode where the angle of the PBS has a fixed value and
% multiple simulations are conducted in order to assess the power levels at
% the output ports as a function of the PBS angle.
% 2. A dynamic mode where the PBS angle is made to vary continuously over 
% the time window of the simulation. The power levels ar the output ports
% of the PBS can then be monitored over a single simulation run.
%
% Christophe Peucheret (christophe.peucheret@univ-rennes1.fr)
% 2024-08-01
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
% results.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'test','fig');
file_name_core_data = strrep(mfilename,'test','data');
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');


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
% global CONSTANT
% global space_grid

% -------------------------------------------------------------------------
% Set global simulation parameters
% -------------------------------------------------------------------------
reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
nsymbols = 128;
symbol_rate = 40e9;

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
% Load essential physical constants.

% ------------------------------------------------------------------------- 
% Start time 
% ------------------------------------------------------------------------- 
start_time = clock;
fprintf('\n\n%s%s\n\n','Simulation started on ',datestr(start_time));


%--------------------------------------------------------------------------
% Switches
%--------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;
fig.interpreter = 'latex';



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
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------



params_cw.power = 1.0e-3;
params_cw.linewidth = 0;
params_cw.emission_frequency = 193.1e12;
[sig,params_cw.actual_emission_frequency] = opt_laser_cw(params_cw);
% CW laser

power_in = char_opt_average_power(sig);
% Check laser power = power at input of polarizer

nangles = 20;
% Number of PBS angles

theta_static = 2*pi/nangles*(0:nangles - 1);
% Polarizer angles

pbs_il = 3;
% Insertion loss, in dB
pbs_er = 20;
% Polarization extinction ratio, in dB

power_static_i = zeros(1,nangles);
power_static_j = zeros(1,nangles);
% Pre-allocate arrays in which the power at the output ports of the PBS
% will be stored.

% First we use a static implementation of the PBS
for iangle = 1:nangles
    % Loop over PBS angles
    
    [sig_i,sig_j] = pol_pbs(sig,theta_static(iangle),pbs_il,pbs_er);
    % Polarization beam splitter
    
    power_static_i(iangle) = char_opt_average_power(sig_i);
    power_static_j(iangle) = char_opt_average_power(sig_j);
    
end

power_static_tot = power_static_i + power_static_j;

% Second, we use a dynamic implementation of the PBS
theta_dynamic = 2*pi/nsamples*(0:nsamples - 1);

[sig_i,sig_j] = pol_pbs(sig,theta_dynamic,pbs_il,pbs_er);
% Polarization beam splitter

power_dynamic_i = abs(sig_i.x).^2 + abs(sig_i.y).^2;
power_dynamic_j = abs(sig_j.x).^2 + abs(sig_j.y).^2;

power_dynamic_tot = power_dynamic_i + power_dynamic_j;


fig_name = [file_name_core_figure '_transmission_lin'];
hfig = figure('Name',fig_name);
plot(theta_dynamic/pi,power_dynamic_i/power_in,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(theta_static/pi,power_static_i/power_in,'Color','b','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_dynamic/pi,power_dynamic_j/power_in,'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_static/pi,power_static_j/power_in,'Color','r','LineStyle','none','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_dynamic/pi,power_dynamic_tot/power_in,'Color','g','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_static/pi,power_static_tot/power_in,'Color','g','LineStyle','none','Marker','^','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('PBS angle ($\times\pi$ rad) ','Interpreter',fig.interpreter)
ylabel('transmission','Interpreter',fig.interpreter)
legend('port $\mathbf{i}$ - dyn.','port $\mathbf{i}$ - stat.', 'port $\mathbf{j}$ - dyn.', 'port $\mathbf{j}$ - stat.','total - dyn.','total - stat.','Location','North','Box','off','Interpreter',fig.interpreter,'FontSize',15,'NumColumns',3)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 2])
ylim([-0.05 1.4])

hfig.Position = [680 678 560*1.5 420];






fig_name = [file_name_core_figure '_transmission_log'];
hfig = figure('Name',fig_name);
plot(theta_dynamic/pi,10*log10(power_dynamic_i/power_in),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(theta_static/pi,10*log10(power_static_i/power_in),'Color','b','LineStyle','none','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_dynamic/pi,10*log10(power_dynamic_j/power_in),'Color','r','LineStyle','--','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_static/pi,10*log10(power_static_j/power_in),'Color','r','LineStyle','none','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_dynamic/pi,10*log10(power_dynamic_tot/power_in),'Color','g','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta_static/pi,10*log10(power_static_tot/power_in),'Color','g','LineStyle','none','Marker','^','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('PBS angle ($\times\pi$ rad) ','Interpreter',fig.interpreter)
ylabel('transmission (dB)','Interpreter',fig.interpreter)
legend('port $\mathbf{i}$ - dyn.','port $\mathbf{i}$ - stat.', 'port $\mathbf{j}$ - dyn.', 'port $\mathbf{j}$ - stat.','total - dyn.','total - stat.','Location','North','Box','off','Interpreter',fig.interpreter,'FontSize',15,'NumColumns',3)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 2])
ylim([-40 15])

hfig.Position = [680 678 560*1.5 420];









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
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% Display duration (or wrapup: core_wrapup(start_time,0);)
% ------------------------------------------------------------------------- 
core_display_duration(start_time,clock);
% ------------------------------------------------------------------------- 
% End of core file
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 
% ------------------------------------------------------------------------- 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Guang Tong Xin Xi Tong Fang Zhen
% C. Peucheret (christophe.peucheret@univ-rennes.fr) 2009-20xx
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------