% -------------------------------------------------------------------------
% Test of the ber_gauss.m function.
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Received power range
% -------------------------------------------------------------------------
prx_min = -30;
prx_max = -15;
prx_step = 1;
% dBm

% -------------------------------------------------------------------------
% Receiver parameters
% -------------------------------------------------------------------------
params_pin.pd.responsivity = 1;
params_pin.pd.thermal_noise_density = 15e-12;
params_pin.pd.shot_noise = 0;
params_pin.pd.dark_current = 0;
params_pin.elpf.type = 'bessel';
params_pin.elpf.order = 4;
params_pin.elpf.f3dB = 0.7*symbol_rate;

% -------------------------------------------------------------------------
% PRBS
% -------------------------------------------------------------------------
params_prbs.type = 'shift_register';
params_prbs.order = 7;
params_prbs.poly = [7 6 0];%[7 3 0];[7 1 0];
params_prbs.seed = [1 1 1 0 1 1 0];
        
        
% -------------------------------------------------------------------------        
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Here we go
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% -------------------------------------------------------------------------

bit_pattern = logical_prbs(params_prbs);
bit_pattern = logical_adapt_binary_sequence(bit_pattern,nsymbols);
% Generate random sequence and adapt length to simulation time window


params_tx.type = 'ook_nrz';
params_tx.emission_frequency = reference_frequency;
params_tx.power = 1.0e-3;
params_tx.linewidth = 0;
params_tx.bit_pattern = bit_pattern;
params_tx.rise_time = 1/symbol_rate/4;
sig_tx = tx(params_tx);
% Generate OOK NRZ signal


prx_dbm = [prx_min:prx_step:prx_max];
% Received power vector


for ipower = 1:length(prx_dbm)    
    
    params_att.mode = 'power';
    params_att.target = prx_dbm(ipower); 
    sig = opt_att(sig_tx,params_att);    
    % Optical Attenuator
    
    pav = char_opt_average_power(sig);    
    prx_check_dbm = 10*log10(pav/1.0e-3);
    % Received power check   
    
    sig = rx_pin(sig,params_pin);
    % PIN receiver
    
    if do_debug
        params_eye.pol = 'x';%'y','both';
        params_eye.neyes = 2;
        params_eye.samples_per_symbol = nsamples_per_symbol;
        params_eye.save.ascii = 0;
        params_eye.save.emf = 0;
        params_eye.save.jpg = 0;
        params_eye.save.vertical_scale_type = 'auto';%'fixed';
        params_eye.save.vertical_scale_max = 1.0e-3;
        params_eye.save.display_baseline = 1;
        params_eye.colour_grade = 0;
        params_eye.name = ['eye diagram for prx =' num2str(prx_dbm(ipower)) ' dBm'];
        meas_eye(sig,params_eye);
        % Eye diagram
    end
        
    params_bergauss.ignore_bits_start = 0;
    params_bergauss.distribution = 'gauss';
    params_bergauss.threshold_mode = 'optimum';%'optimum_search';%'fixed';
    params_bergauss.threshold = [ ];
    params_bergauss.sample_mode = 'optimum';%'fixed';
    params_bergauss.sample_index = nsamples_per_symbol/2;
    [ber_gauss_ber(ipower),ber_gauss_threshold(ipower),ber_gauss_sample(ipower)] = ber_gauss(sig,bit_pattern,params_bergauss);


end
% End of loop over received power.








%%
% -------------------------------------------------------------------------
% Analytical calculations of the BER according to Gaussian noise theory
% -------------------------------------------------------------------------
    hel = elec_tf_elpf(params_pin.elpf,frequency_array);
    hel0 = elec_tf_elpf(params_pin.elpf,0);
    enb = calc_enb('lowpass',hel,hel0,df);

%enb = calc_enb('elec',params_pin.elpf);
% Noise equivalent bandwidth of the receiver.

sig2ther = params_pin.pd.thermal_noise_density^2*enb;
% Thermal noise variance

Q = params_pin.pd.responsivity*2*10.^(prx_dbm/10)*1.0e-3/(2*sqrt(sig2ther));
% Q-factor

ber_theory  = 0.5*erfc(Q/sqrt(2));
% Conversion from Q-factor to BER.


% -------------------------------------------------------------------------
% Retrieve Monte-Carlo error counting results from file
% -------------------------------------------------------------------------
% data_ec = io_read_text_file('./data/data_test_ber_gauss_mc_ber.dat','%f %f %f %f %f %f %f %f %f %f %f %f',1);
% prx_ec = data_ec(:,1);
% ber_ec = data_ec(:,3);


%%
% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

fig_name = [file_name_core_figure '_compare_methods'];
hfig = figure('Name',fig_name);
h1 = plot(prx_dbm,log10(ber_gauss_ber),'Color','b','LineStyle','-','Marker','s','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot(prx_dbm,log10(ber_theory),'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
% h3 = plot(prx_ec,log10(ber_ec),'Color','k','LineStyle','-','Marker','o','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hline(-9,'r--');
xlabel('average received power (dBm)','Interpreter',fig.interpreter)
ylabel('$\log_{10}(\textrm{BER})$','Interpreter',fig.interpreter)
% legend([h1 h3 h2],{'Gaussian estimation' 'error counting' 'Gaussian theory'},'Location','West','Box','off','Interpreter',fig.interpreter,'FontSize',15)
legend([h1 h2],{'Gaussian estimation' 'Gaussian theory'},'Location','West','Box','off','Interpreter',fig.interpreter,'FontSize',15)

% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([-30 -21])
% ylim([-0.5 0.05])

% hfig.PaperUnits = 'centimeters';
% hfig.PaperType = 'a4';
% hfig.PaperOrientation = 'landscape';
% hfig.PaperPosition = [5 5 17 10];

if do_print
    print(fig_name,'-dmeta');
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
