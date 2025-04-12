% -------------------------------------------------------------------------
% Test the definition of standard digital constellations:
%
% Each defined constellation is plotted with indication of the
% corresponding word (binary and decimal representations) mapped to each
% point of the constellation in the complex plane.
%
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
time_stamp = datestr(datetime('now','TimeZone','Z'),'yyyymmddThhMMSSZ');

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
% Set and plot constellations
% -------------------------------------------------------------------------
params_plot_constellation.labels.dec = 1;
params_plot_constellation.labels.bin = 1;
params_plot_constellation.axes = 0;
params_plot_constellation.save = 0;


% ----
% BPSK
% ----
constellation_order = 2;
constellation_type = 'bpsk';
[constellation_bpsk,norm_es_bpsk,norm_emax_bpsk] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_bpsk,params_plot_constellation);

% ----
% PAM4, natural coding
% ----
constellation_order = 4;
constellation_type = 'pam4_natural';
[constellation_pam4_natural,norm_es_pam4_natural,norm_emax_pam4_natural] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_pam4_natural,params_plot_constellation);

% ----
% PAM4, Gray coding
% ----
constellation_order = 4;
constellation_type = 'pam4_gray';
[constellation_pam4_gray,norm_es_pam4_gray,norm_emax_pam4_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_pam4_gray,params_plot_constellation);

% ----
% QPSK, natural coding
% ----
constellation_order = 4;
constellation_type = 'qpsk_natural';
[constellation_qpsk_natural,norm_es_qpsk_natural,norm_emax_qpsk_natural] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qpsk_natural,params_plot_constellation);

% ----
% QPSK, Gray coding
% ----
constellation_order = 4;
constellation_type = 'qpsk_gray';
[constellation_qpsk_gray,norm_es_qpsk_gray,norm_emax_qpsk_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qpsk_gray,params_plot_constellation);

% ----
% PSK8, Gray coding
% ----
constellation_order = 8;
constellation_type = 'psk8_gray';
[constellation_psk8_gray,norm_es_psk8_gray,norm_emax_psk8_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_psk8_gray,params_plot_constellation);

% ----
% QAM8, rectangular, Gray coding
% ----
constellation_order = 8;
constellation_type = 'qam8_rect_gray';
[constellation_qam8_rect_gray,norm_es_qam8_rect_gray,norm_emax_qam8_rect_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam8_rect_gray,params_plot_constellation);

% ----
% QAM8, square, Gray coding
% ----
constellation_order = 8;
constellation_type = 'qam8_square_gray';
[constellation_qam8_square_gray,norm_es_qam8_square_gray,norm_emax_qam8_square_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam8_square_gray,params_plot_constellation);

% ----
% QAM8, circular
% ----
constellation_order = 8;
constellation_type = 'qam8_circ';
[constellation_qam8_circ,norm_es_qam8_circ,norm_emax_qam8_circ] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam8_circ,params_plot_constellation);

% ----
% QAM8, star
% ----
constellation_order = 8;
constellation_type = 'qam8_star';
[constellation_qam8_star,norm_es_qam8_star,norm_emax_qam8_star] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam8_star,params_plot_constellation);


% ----
% QAM16, natural coding
% ----
constellation_order = 16;
constellation_type = 'qam16_natural';
[constellation_qam8_natural,norm_es_qam8_natural,norm_emax_qam8_natural] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam8_natural,params_plot_constellation);


% ----
% QAM16, Gray coding
% ----
constellation_order = 16;
constellation_type = 'qam16_gray';
[constellation_qam16_gray,norm_es_qam16_gray,norm_emax_qam16_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam16_gray,params_plot_constellation);

% ----
% QAM16, Quadrant coding
% ----
constellation_order = 16;
constellation_type = 'qam16_quadrant';
[constellation_qam16_quadrant,norm_es_qam16_quadrant,norm_emax_qam16_quadrant] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam16_quadrant,params_plot_constellation);

% ----
% QAM32, rectangular, Gray mapping
% ----
constellation_order = 32;
constellation_type = 'qam32_rect_gray';
[constellation_qam32_rect_gray,norm_es_qam32_rect_gray,norm_emax_qam32_rect_gray] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam32_rect_gray,params_plot_constellation);

% ----
% QAM32, cross
% ----
constellation_order = 32;
constellation_type = 'qam32_cross';
[constellation_qam32_cross,norm_es_qam32_cross,norm_emax_qam32_cross] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam32_cross,params_plot_constellation);

% ----
% QAM64, Gray mapping
% ----
constellation_order = 64;
constellation_type = 'qam64_gray';
[constellation_qam64_cross,norm_es_qam64_cross,norm_emax_qam64_cross] = define_constellation(constellation_type,constellation_order);
params_plot_constellation.name = constellation_type;
plot_constellation_mapping(constellation_qam64_cross,params_plot_constellation);



% -------------------------------------------------------------------------
% Figure alignment
% -------------------------------------------------------------------------
% Requires alignfigs function
% Available at https://github.com/nickhale/alignfigs
alignfigs(4)


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

