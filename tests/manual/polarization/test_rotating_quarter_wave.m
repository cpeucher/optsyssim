% -------------------------------------------------------------------------
% We reproduce the Stokes parameter measurement method based on a rotating
% quarter-wave plate described in e.g.
% B. Schaefer, et al., “Measuring the Stokes polarization parameters,” 
% Am. J. Phys. 75, 163 (2007). doi: 10.1119/1.2386162
% and earlier work.
%
% Here we implement a proof-of-concept simulation based on an ideal setup:
% - The initial orientation of the slow axis of the quarter-wave plate is
% known. It is taken equal to zero at the time origin.
% - The retardation of the quarter-wave plate is exactly pi/2
% - The orientation of the linear polarizer is exactly along the x axis of
% the canonical basis
% In reality, some calibration procedure would be needed, as reported in a
% number of references in the literature.
%
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
% Choice of colormap
% -------------------------------------------------------------------------
line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));

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
global CONSTANT
% global space_grid

% -------------------------------------------------------------------------
% Set global simulation parameters
% -------------------------------------------------------------------------
reference_frequency = 193.1e12;
nsamples_per_symbol = 128;
nsymbols = 32;
symbol_rate = 25e9;

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

%%
% -------------------------------------------------------------------------
% Signal generation
% -------------------------------------------------------------------------
alpha = pi/4;
delta = pi/2;
% Right-handed circular polarization

alpha = pi/8;
delta = pi/6;
% Somehow arbitrary polarization

sig.x = cos(alpha)*ones(1,nsamples);
sig.y = sin(alpha)*ones(1,nsamples).*exp(1j*delta);
% Generate signal with desired state-of-polarization


fprintf('\n\n%s\n','Expected normalised Stokes parameters:')

S_exp = [1 cos(2*alpha) sin(2*alpha)*cos(delta) sin(2*alpha)*sin(delta)].'
% Expected normalised Stokes vector.

%%
% -------------------------------------------------------------------------
% Characterization of Stokes vector using standard method with linear and
% circular polarizers
% -------------------------------------------------------------------------
S = char_opt_stokes(sig);
% Instantaneous Stokes vector


fprintf('\n%s\n','Parameters retrieved from conventional Jones to Stokes conversion:')

S_std = S(:,1)/S(1,1)
% Actual normalised Stokes vector



% -------------------------------------------------------------------------
% Characterization of Stokes vector using the rotating quarter waveplate
% method
% -------------------------------------------------------------------------
nperiods = 16;
% Number of periods of the detected signal.
% Choose an even number.
% The number of rotation of the quarter wave plate is nperiods/2.

theta0 = 0;
% Angle of the slow-axis of the quarter-wave plate with respect to the -x
% axis of the canonical basis at the origin.
% In this proof-of-concept exemple, recovery of the Stokes parameters 
% from the FFT of the detected signal only works when theta0 = 0.
% Otherwise some calibration needs to be implemented. 
% See e.g. M. J. Romerein et al., “Calibration method using a single 
% retarder to simultaneously measure polarization and fully characterize 
% a polarimeter over a broad range of wavelengths,” Applied Optics 50,5382
% (2011). doi: 10.1364/AO.50.005382
% and numerous other references.


theta = nperiods*pi/nsamples*(0:nsamples - 1) - theta0;
% Wave plate angle. 

fm = nperiods*sample_rate/2/nsamples;
% Rotation frequency of quarter-wave plate

sig = pol_retarder(sig,theta,pi/2);
% Rotating quarter wave plate



S_wp = char_opt_stokes(sig);
hfig = plot_poincare_sphere([file_name_core_figure '_qwp_op']);
plot3(S_wp(2,:)./S_wp(1,:),S_wp(3,:)./S_wp(1,:),S_wp(4,:)./S_wp(1,:),'LineStyle','-','LineWidth',3,'Color','r','Marker','none','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r')
% We plot the evolution of the state-of-polarization after the quarter
% waveplate on the Poincare sphere.
plot3(S_exp(2)/S_exp(1),S_exp(3)/S_exp(1),S_exp(4)/S_exp(1),'Marker','s','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
% Indicate the SOP one wishes to determine on the Poincare sphere


sig = pol_polarizer(sig,0,0,Inf);
% Linear horizontal polariser

S_pol = char_opt_stokes(sig);
plot3(S_pol(2,:)./S_pol(1,:),S_pol(3,:)./S_pol(1,:),S_pol(4,:)./S_pol(1,:),'Marker','o','MarkerSize',15,'MarkerEdgeColor','b','MarkerFaceColor','b');
% Display state-of-polarization after linear polarizer (!) on Poincare
% sphere.


fig_name = [file_name_core_figure '_stokes'];
hfig = figure('Name',fig_name);
plot(theta/pi,S_wp(1,:)./S_wp(1,:),'Color',line_color(1,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
hold on
plot(theta/pi,S_wp(2,:)./S_wp(1,:),'Color',line_color(2,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta/pi,S_wp(3,:)./S_wp(1,:),'Color',line_color(3,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
plot(theta/pi,S_wp(4,:)./S_wp(1,:),'Color',line_color(4,:),'LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('quarter wave plate angle ($\times\pi$ rad)','Interpreter',fig.interpreter)
ylabel('normalized Stokes parameter','Interpreter',fig.interpreter)
legend('$S_0$','$S_1$','$S_2$','$S_3$','Location','SouthWest','Box','off','Interpreter',fig.interpreter,'FontSize',15,'NumColumns',2)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([0 2])
ylim([-1.2 1.2])


id = abs(sig.x).^2;
% Signal detection
% The linear horizontal polarizer being ideal (infinite polarization
% extinction ratio),we do not expect (and verify...) any power along the -y
% axis.


fig_name = [file_name_core_figure '_detected_time'];
hfig = figure('Name',fig_name);
plot(theta/pi,id,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('quarter wave plate angle ($\times\pi$ rad)','Interpreter',fig.interpreter)
ylabel('detected signal (a.u.)','Interpreter',fig.interpreter)
%legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])


%%
% -------------------------------------------------------------------------
% We recover the Stokes parameters from the spectrum (FFT) of the detected 
% signal.
% This assumes that theta0 = 0, i.e. at the origin the slow axis of the
% quarter-wave plate is aligned with -x
% -------------------------------------------------------------------------

id_fft = fft(id)/nsamples;
% Calculate spectrum

fig_name = [file_name_core_figure '_detected_freq'];
hfig = figure('Name',fig_name);
semilogy(frequency_array/fm,fftshift(abs(id_fft)),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel('frequency ($\times f_m$ Hz)','Interpreter',fig.interpreter)
ylabel('magnitude','Interpreter',fig.interpreter)
%legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([0 8])
% ylim([-0.5 0.05])


Id_dc = id_fft(1);
% Value of the frequency component at DC. Should be real.
Id_2fm = id_fft(nperiods + 1);
% Value of the frequency component at 2*fm, where fm is the rotation
% frequency of the wave plate. Should be imaginary.
Id_4fm = id_fft(2*nperiods + 1);
% Value of the frequency component at 4*fm. Should be complex.

A_spec = 2*Id_dc;
% Half the amplitude of the cosntant term in the expansion of id
B_spec = 4*1i*Id_2fm;
% Half the amplitude of the sin(2*theta) term in the expansion of id
C_spec = 4*real(Id_4fm);
% Half the amplitude of the cos(4*theta) term in the expansion of id
D_spec = -4*imag(Id_4fm);
% Half the amplitude of the sin(4*theta) term in the expansion of id

% From which we can recover the Stokes parameters:
S_spec = [A_spec - C_spec; 2*C_spec; 2*D_spec; -B_spec];
% Stokes parameters
% Observe that S3 = -B
% There is a minus sign in the relation between S3 and B, unlike what is
% written in the paper by Schaefer et al. (2007). 

fprintf('\n%s\n','Stokes parameters retrieved from rotating quarter wave plate method:')
fprintf('%s\n','Spectral analysis of detected signal')
fprintf('%s\n','Valid only when theta_0 = 0')

S_spec = S_spec/S_spec(1)
% Normalized Stokes parameters.

%%
% -------------------------------------------------------------------------
% We can also recover the Stokes parameter through relations (20) in
% (Schaefer et al., 2007).
% -------------------------------------------------------------------------

A = 2*sum(id)/nsamples;
B = 4*sum(id.*sin(2*theta))/nsamples;
C = 4*sum(id.*cos(4*theta))/nsamples;
D = 4*sum(id.*sin(4*theta))/nsamples;

S_qwp = [A - C; 2*C; 2*D; -B];
% Stokes parameters.

fprintf('\n%s\n','Stokes parameters retrieved from rotating quarter wave plate method:')
fprintf('%s\n','Time-domain "sine/cosine transform" analysis')

dop_qwp = calc_dop(S_qwp)
% Degree of polarization

S_qwp = S_qwp/S_qwp(1)
% Normalized Stokes parameters




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