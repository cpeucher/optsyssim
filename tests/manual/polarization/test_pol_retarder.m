% -------------------------------------------------------------------------
% Test of pol_retarder function
% We look at the transformation of given input states-of-polarization by
% wave plates (also known as retarders).
%
% Christophe Peucheret (christophe.peucheret@univ-rennes1.fr)
% 2024-08-02
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
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
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
nsamples_per_symbol = 32;
nsymbols = 128;
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



% -------------------------------------------------------------------------
% Test 1
% Tranformation of linear polarization at +45 degrees by a quarter wave
% plate whose slow axis is along -x
% We expect a left-hand circular polarization (LHCP)
% -------------------------------------------------------------------------
fprintf('\n\n%s','Test 1')
fprintf('\n%s\n','======')
fprintf('\n%s\t%s','PLATE:','Quarter wave plate with slow axis along -x')
fprintf('\n%s\t\t%s','IN:','Linear polarization at +45 degrees')
fprintf('\n%s\t%s\n','OUT:','Left-hand circular polarization')


alpha = pi/4;
delta = 0;
sig_in.x = cos(alpha)*ones(1,nsamples);
sig_in.y = sin(alpha)*ones(1,nsamples).*exp(1j*delta);
% Generate linearly polarized signal at 45 degrees.

expSin = [1 0 1 0].'
% Expected normalised Stokes vector at input

calcSin = [1 cos(2*alpha) sin(2*alpha)*cos(delta) sin(2*alpha)*sin(delta)].'
% Calculated normalised Stokes vector at input

S = char_opt_stokes(sig_in);
actSin = S(:,1)/S(1,1)
% Actual normalised Stokes vector at input

sig_out = pol_retarder(sig_in,0,pi/2);
% Quarter wave plate with slow axis along x


expSout = [1 0 0 -1].'
% Expected Stokes vector at output
% Corrsponds to LHCP.

S = char_opt_stokes(sig_out);
actSout = S(:,1)/S(1,1)
% Actual normalised Stokes vector at output


% Representation of the transformation on the Poincare sphere: 
sig_out = pol_retarder(sig_in,0,pi*(0:1:nsamples-1)/nsamples/2);
% Quarter wave plate with slow axis along x
S = char_opt_stokes(sig_out);

fig_name = [file_name_core_figure '_quarter_L_to_LHC'];
hfig = plot_poincare_sphere(fig_name);
plot3(S(2,:)./S(1,:),S(3,:)./S(1,:),S(4,:)./S(1,:),'LineStyle','-','LineWidth',3,'Color','r','Marker','none','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r')

% We can also obtain right-handed circular polarization by rotating the
% axis of the quarter wave plate:
sig_out = pol_retarder(sig_in,pi/2,pi*(0:1:nsamples-1)/nsamples/2);
% Quarter wave plate with slow axis along x
S = char_opt_stokes(sig_out);

fig_name = [file_name_core_figure '_quarter_L_to_RHC'];
hfig = plot_poincare_sphere(fig_name);
plot3(S(2,:)./S(1,:),S(3,:)./S(1,:),S(4,:)./S(1,:),'LineStyle','-','LineWidth',3,'Color','r','Marker','none','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r')



%%
% -------------------------------------------------------------------------
% Test 2
% Tranformation of right-hand circular polarization by a quarter
% wave plate whose slow axis is along -x
% We expect a linear polarization at +45 degrees
% -------------------------------------------------------------------------
fprintf('\n\n%s','Test 2')
fprintf('\n%s\n','======')
fprintf('\n%s\t%s','PLATE:','Quarter wave plate with slow axis along -x')
fprintf('\n%s\t\t%s','IN:','Right-hand circular polarization')
fprintf('\n%s\t%s\n','OUT:','Linear polarization at +45 degrees')

alpha = pi/4;
delta = pi/2;
sig_in.x = cos(alpha)*ones(1,nsamples);
sig_in.y = sin(alpha)*ones(1,nsamples).*exp(1j*delta);
% Generate right-hand circular polarization (RHCP)

expSin = [1 0 0 1].'
% Expected normalised Stokes vector at input

calcSin = [1 cos(2*alpha) sin(2*alpha)*cos(delta) sin(2*alpha)*sin(delta)].'
% Calculated normalised Stokes vector at input

S = char_opt_stokes(sig_in);
actSin = S(:,1)/S(1,1)
% Actual normalised Stokes vector at input

sig_out = pol_retarder(sig_in,0,pi/2);
% Quarter wave plate with slow axis along x


expSout = [1 0 1 0].'
% Expected Stokes vector at output
% Corresponds to linear at +45 degrees

S = char_opt_stokes(sig_out);
actSout = S(:,1)/S(1,1)
% Actual normalised Stokes vector at output



% -------------------------------------------------------------------------
% Test 3
% Tranformation of linear polarization at 45 degrees by a half wave plate 
% whose slow axis is along -x
% We expect a linear polarization at +135 degrees
% -------------------------------------------------------------------------
fprintf('\n\n%s','Test 3')
fprintf('\n%s\n','======')
fprintf('\n%s\t%s','PLATE:','Half wave plate with slow axis along -x')
fprintf('\n%s\t\t%s','IN:','Linear polarization at +45 degrees')
fprintf('\n%s\t%s\n','OUT:','Linear polarization at +135 degrees')


alpha = pi/4;
delta = 0;
sig_in.x = cos(alpha)*ones(1,nsamples);
sig_in.y = sin(alpha)*ones(1,nsamples).*exp(1j*delta);
% Generate linearly polarized signal at 45 degrees.

expSin = [1 0 1 0].'
% Expected normalised Stokes vector at input

calcSin = [1 cos(2*alpha) sin(2*alpha)*cos(delta) sin(2*alpha)*sin(delta)].'
% Calculated normalised Stokes vector at input

S = char_opt_stokes(sig_in);
actSin = S(:,1)/S(1,1)
% Actual normalised Stokes vector at input

sig_out = pol_retarder(sig_in,0,pi);
% Quarter wave plate with slow axis along x


expSout = [1 0 -1 0].'
% Expected Stokes vector at output
% Corrsponds to LHCP.

S = char_opt_stokes(sig_out);
actSout = S(:,1)/S(1,1)
% Actual normalised Stokes vector at output


%%
% -------------------------------------------------------------------------
% Test 4
% Tranformation of right-hand circular polarization by a half wave plate 
% whose slow axis is along -x
% We expect a left-hand circular polarization
% -------------------------------------------------------------------------
fprintf('\n\n%s','Test 4')
fprintf('\n%s\n','======')
fprintf('\n%s\t%s','PLATE:','Half wave plate with slow axis along -x')
fprintf('\n%s\t\t%s','IN:','Right-hand circular polarization')
fprintf('\n%s\t%s\n','OUT:','Left-hand circular polarization')

alpha = pi/4;
delta = pi/2;
sig_in.x = cos(alpha)*ones(1,nsamples);
sig_in.y = sin(alpha)*ones(1,nsamples).*exp(1j*delta);
% Generate right-hand circular polarization (RHCP)

expSin = [1 0 0 1].'
% Expected normalised Stokes vector at input

calcSin = [1 cos(2*alpha) sin(2*alpha)*cos(delta) sin(2*alpha)*sin(delta)].'
% Calculated normalised Stokes vector at input

S = char_opt_stokes(sig_in);
actSin = S(:,1)/S(1,1)
% Actual normalised Stokes vector at input

sig_out = pol_retarder(sig_in,0,pi);
% Quarter wave plate with slow axis along x


expSout = [1 0 0 -1].'
% Expected Stokes vector at output
% Corresponds to linear at +45 degrees.

S = char_opt_stokes(sig_out);
actSout = S(:,1)/S(1,1)
% Actual normalised Stokes vector at output










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