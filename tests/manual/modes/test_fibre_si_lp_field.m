% -------------------------------------------------------------------------
% Test of fibre_si_lp_field function
%
% Representation of the transverse field and power distribution of the LPlm
% modes of a weakly-guiding, circular-core, step-index optical fibre.
%
% This includes:
% - Determining the normalised propagation constant via fibre_si_b, which
%   in turns uses the fibre_si_dispersion_equation function.
% - Calculating the field via fibre_si_lp_field 
% - Plotting the field / power distribution via fibre_si_plot_mode
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all
format long

% -------------------------------------------------------------------------
% Strings for filenames
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'run','fig');
file_name_core_data = strrep(mfilename,'run','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));

% ------------------------------------------------------------------------- 
% Reinitialise the random number generator for reproducibility.
% ------------------------------------------------------------------------- 
stream = RandStream.getGlobalStream;
reset(stream);

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
% 
% -------------------------------------------------------------------------

xrange = [-1, 1]*65e-6;  
yrange = [-1, 1]*65e-6; 
nxpoints = 2001;
nypoints = 2001;
global space_grid
space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
[space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 

params_fibre.n1 =  1.45;
params_fibre.n2 =  1.4361;
params_fibre.a = 25e-6;
mode_type = 'LP';
mode_l = 1;             % Index l of mode LPlm
mode_m = 1;             % Index m of mode LPlm. Check that mode_m <= nmodes
lambda = 850e-9; 
V = 2*pi*params_fibre.a*sqrt(params_fibre.n1^2 - params_fibre.n2^2)/lambda;
[b,nmodes] = fibre_si_b('LP',mode_l,V,params_fibre);
% Normalised propagation constant
field = fibre_si_lp_field(params_fibre,mode_l,b(mode_m),V,'even'); 
% Mode field

visparams.limit_radius = Inf;
visparams.show_core_limit = 1;
visparams.show_core_linewidth = 5;
visparams.save = 0;
visparams.colormap = 'jet';%'hot';
visparams.name = ['Mode field distribution for ' mode_type num2str(mode_l) num2str(mode_m)];
fibre_si_plot_mode(field,params_fibre,visparams);


visparams.limit_radius = 62.5e-6;
visparams.show_core_limit = 1;
visparams.show_core_linewidth = 5;
visparams.save = 1;
visparams.colormap = 'hot';
visparams.name = ['Mode power distribution for ' mode_type num2str(mode_l) num2str(mode_m)];
fibre_si_plot_mode(abs(field).^2,params_fibre,visparams);


% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
% alignfigs(2)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

