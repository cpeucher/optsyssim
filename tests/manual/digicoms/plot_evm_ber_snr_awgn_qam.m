% -------------------------------------------------------------------------
% Plot relations between EVM, BER, SER and SNR over an AWGN channel for 
% QPSK, 16-QAM and 64-QAM modulation.
% Different normalisations are considered for EVM calculations.
%
% Run:
% test_evm_ber_snr_awgn_qpsk.m
% test_evm_ber_snr_awgn_qam16.m
% test_evm_ber_snr_awgn_qam64.m
% to generate the necessary data files.
%
% Otherwise some sample data files are provided in
% /data/digicoms/
% under the names of: 
% data_evm_ber_snr_awgn_qpsk.mat
% data_evm_ber_snr_awgn_qam16.mat
% data_evm_ber_snr_awgn_qam64.mat
% 
% In particular we attempt to replicate the relation between BER and EVM
% (for the AWGN channel) in
% R. Schmogrow et al., "Error vector magnitude as a performance measure for 
% advanced modulation formats," IEEE Photonics Technology Letters 24, 61 
% (2012) [DOI: 10.1109/LPT.2011.2172405].

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
do_add_figsize_to_filename = 0;
margin_figure = 0;

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.font_name = 'times';
fig.font_size = 15;
fig.interpreter = 'latex';


mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));



% -------------------------------------------------------------------------
% LOAD DATA
% -------------------------------------------------------------------------
data_path = '../../../data/digicoms/';


data_file_name = 'data_evm_ber_snr_awgn_qpsk.mat';
data_qam4 = load([data_path data_file_name]);
% Load Matlab mat file

data_file_name = 'data_evm_ber_snr_awgn_qam16.mat';
data_qam16 = load([data_path data_file_name]);
% Load Matlab mat file

data_file_name = 'data_evm_ber_snr_awgn_qam64.mat';
data_qam64 = load([data_path data_file_name]);
% Load Matlab mat file


%%
% -------------------------------------------------------------------------
% EVM vs SNRmax
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_evmmax_vs_snr_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
plot(data_qam4.esn0_db,100*data_qam4.evm_max_da,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
plot(data_qam16.esn0_db,100*data_qam16.evm_max_da,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(data_qam64.esn0_db,100*data_qam64.evm_max_da,'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

hleg = legend('QPSK','16-QAM','64-QAM');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s/N_0$ (dB)');
y1 = ylabel('$\textrm{EVM}_\textrm{max}$ (\%)');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;

ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
% xlim([0 57]);
% ax.XTick = [0:10:60];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([0 55]);
ax.YTick = [0:10:50];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end




%%
% -------------------------------------------------------------------------
% EVM vs SNRrms
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_evmrms_vs_snr_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
plot(data_qam4.esn0_db,100*data_qam4.evm_rms_da,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
plot(data_qam16.esn0_db,100*data_qam16.evm_rms_da,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(data_qam64.esn0_db,100*data_qam64.evm_rms_da,'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot([0:0.1:30],sqrt(1./10.^([0:0.1:30]/10)),'LineWidth',1.5,'Color','k','LineStyle','none','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');


hleg = legend('QPSK','16-QAM','64-QAM');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s/N_0$ (dB)');
y1 = ylabel('$\textrm{EVM}_\textrm{rms}$ (\%)');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;

ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
% xlim([0 57]);
% ax.XTick = [0:10:60];


ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([0 55]);
ax.YTick = [0:10:50];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end






%%
% -------------------------------------------------------------------------
% BER vs EVM
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ber_vs_evmmax_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
plot(100*data_qam4.evm_max_da,log10(data_qam4.ber),'LineWidth',1.5,'Color','r','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
plot(100*data_qam16.evm_max_da,log10(data_qam16.ber),'LineWidth',1.5,'Color','b','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(100*data_qam64.evm_max_da,log10(data_qam64.ber),'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

hleg = legend('QPSK','16-QAM','64-QAM');
set(hleg,'Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$\textrm{EVM}_\textrm{max}$ (\%)');
y1 = ylabel('$10\log_{10}\left(\textrm{BER}\right)$');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;


ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([0 57]);
ax.XTick = [0:10:60];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([-10 0]);
ax.YTick = [-10:2:0];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end



%%
% -------------------------------------------------------------------------
% SER vs EVM
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ser_vs_evmmax_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
plot(100*data_qam4.evm_max_da,log10(data_qam4.ser),'LineWidth',1.5,'Color','r','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
plot(100*data_qam16.evm_max_da,log10(data_qam16.ser),'LineWidth',1.5,'Color','b','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(100*data_qam64.evm_max_da,log10(data_qam64.ser),'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

hleg = legend('QPSK','16-QAM','64-QAM');
set(hleg,'Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$\textrm{EVM}_\textrm{max}$ (\%)');
y1 = ylabel('$10\log_{10}\left(\textrm{SER}\right)$');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;

ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([0 57]);
ax.XTick = [0:10:60];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([-10 1]);
ax.YTick = [-10:2:0];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end


%%
% -------------------------------------------------------------------------
% SER vs SNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ser_vs_snr_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
h1 = semilogy(data_qam4.esn0_db,data_qam4.ser,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = semilogy(data_qam4.esn0_db,data_qam4.pse,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h3 = semilogy(data_qam16.esn0_db,data_qam16.ser,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h4 = semilogy(data_qam16.esn0_db,data_qam16.pse,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h5 = semilogy(data_qam64.esn0_db,data_qam64.ser,'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h6 = semilogy(data_qam64.esn0_db,data_qam64.pse,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

hleg = legend([h1 h3 h5],{'QPSK','16-QAM','64-QAM'});
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s / N_0$ (dB)');
y1 = ylabel('SER');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;

ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([0 30]);
ax.XTick = [0:5:30];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-8 5]);
% ax.YTick = [ ];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end


%%
% -------------------------------------------------------------------------
% SER vs SNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ber_vs_snr_qpsk_16qam_64qam';
hfig = figure('Name',fig_name);
h1 = semilogy(data_qam4.esn0_db,data_qam4.ber,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = semilogy(data_qam4.esn0_db,data_qam4.pse/2,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h3 = semilogy(data_qam16.esn0_db,data_qam16.ber,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h4 = semilogy(data_qam16.esn0_db,data_qam16.pse/4,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h5 = semilogy(data_qam64.esn0_db,data_qam64.ber,'LineWidth',1.5,'Color','g','LineStyle','none','Marker','^','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
h6 = semilogy(data_qam64.esn0_db,data_qam64.pse/6,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');



hleg = legend([h1 h3 h5],{'QPSK','16-QAM','64-QAM'});
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s / N_0$ (dB)');
y1 = ylabel('BER');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;


ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([0 30]);
ax.XTick = [0:5:30];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-8 1]);
% ax.YTick = [ ];

grid on;

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

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


% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');
