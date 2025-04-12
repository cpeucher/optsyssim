% -------------------------------------------------------------------------
% Relations between EVM, BEP and SNR over an AWGN channel.
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
fig.font_size = 15;
fig.font_name = 'arial';


mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));






% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^20;
% Number of symbols. Should be power of 2.

% -------------------------------------------------------------------------
% Define constellation
% -------------------------------------------------------------------------
m = 64;
% Constellation order. m points with log2(m) even (square constellation)

[constellation,norm_es,norm_emax] = define_constellation('qam64_gray',m);
% Define the constellation for Gray mapping



% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
data = generate_binary(nsymbols,m);
% Generate binary data

[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words



% -------------------------------------------------------------------------
% Mapping
% -------------------------------------------------------------------------
symbs = mapping(words_dec,constellation);
% Mapping

symbs_ref = symbs;
% Save reference signals


% -------------------------------------------------------------------------
% Eb/N0 range
% -------------------------------------------------------------------------
esn0_db = [1:0.5:20];

% -------------------------------------------------------------------------
% Theoretical symbol-error probabilities
% -------------------------------------------------------------------------
pse = sep_awgn_reference(esn0_db,'qam_square',m);


% -------------------------------------------------------------------------
% BER vs SNR
% -------------------------------------------------------------------------
for isnr = 1:length(esn0_db)
    % Loop over SNR    
    
    symbs_rx = add_awgn(symbs,esn0_db(isnr));
    % Add white Gaussian noise  
    
    symbs_rx = normalise_constellation(symbs_rx,norm_es);
    % Normalise constellation   
    
    
    cluster_da = clustering(symbs_rx,words_dec,m);
    % Data-aided clustering
    
     
    symbs_cx = decision_qam_square_hard(symbs_rx,m);
    % Hard decision
    
    words_dec_rx = demapping(symbs_cx,constellation);
    % Demapping. Returns the decimal representation of the symbols.    
    
    cluster_dd = clustering(symbs_cx,words_dec_rx,m);
    % Decision-directed clustering
    
    
    [evm_max_da(isnr), evm_rms_da(isnr)] = calc_evm(constellation,cluster_da);
    % EVM calculation based on data-aided clustering
    
    [evm_max_dd(isnr), evm_rms_dd(isnr)] = calc_evm(constellation,cluster_da);
    % EVM calculation based on decision-directed clustering    
    
    snr_db(isnr) = extract_snr_constellation(cluster_da);
    % Estimate the SNR based on the constellation.   
    
    
    bits_bin = conv_dec2bin(words_dec_rx,2);
    % Convert to binary representation.
    
    ser(isnr) = calc_ser(symbs_cx,symbs_ref);
    % SER calculation.    
    
    ber(isnr) = calc_ber(bits_bin,data);
    % BER calculation.
    
end   


symbol_rate = 20e9; 
% Symbol rate, in bps
noise_bandwidth = 12.5e9; 
% Noise reference bandwidth, in Hz
osnr_db_20gbd = esn0_db + 10*log10(symbol_rate/noise_bandwidth/2);
% OSNR at the specified symbol rate

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
result_file_name = ['data_evm_ber_snr_awgn_qam64'  '.mat'];
save(result_file_name,'esn0_db','snr_db','ber','ser','evm_max_da','evm_max_dd','evm_rms_da','evm_rms_dd','pse');


%%
%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------




% -------------------------------------------------------------------------
% SER vs SNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ser_vs_snr_qam64';
hfig = figure('Name',fig_name);
semilogy(esn0_db,ser,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
semilogy(esn0_db,pse,'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');


hleg = legend('error counting','theory');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s / N_0$ (dB)');
% yyaxis left
y1 = ylabel('SER / $P_s$');
% yyaxis right
% y2 = ylabel(' ');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;


ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([4.5 19]);
ax.XTick = [5:2:19];

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-7 1]);
% ax.YTick = [ ];

grid on;

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'r';
% ylim([0 1.1]);
% ax.YTick = [ ];

% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]

% print(FigName,'-dmeta');
% print(FigName,'-dpdf');
% Command =['pdfcrop ' FigName '.pdf ' FigName '.pdf'];
% system(Command);




%%
% -------------------------------------------------------------------------
% BER vs SNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ber_vs_snr_qam64_gray';
hfig = figure('Name',fig_name);
% yyaxis left
semilogy(esn0_db,ber,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
semilogy(esn0_db,pse/log2(m),'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');



hleg = legend('error counting','theory');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s / N_0$ (dB)');
% yyaxis left
y1 = ylabel('BER / $P_b$');
% yyaxis right
% y2 = ylabel(' ');

x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;


ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([4.5 19]);
ax.XTick = [5:2:19];

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-7 1]);
% ax.YTick = [ ];

grid on;

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'r';
% ylim([0 1.1]);
% ax.YTick = [ ];

% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]

% print(FigName,'-dmeta');
% print(FigName,'-dpdf');
% Command =['pdfcrop ' FigName '.pdf ' FigName '.pdf'];
% system(Command);


%%
% -------------------------------------------------------------------------
% BER vs EVM
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ber_vs_evmmax_qam64';
fhig = figure('Name',fig_name);
% yyaxis left
plot(100*evm_max_da,log10(ber),'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
% hold on;
% semilogy(EsN0_dB,Ps/log2(M),'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');



hleg = legend('64-QAM');
set(hleg,'Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$\textrm{EVM}_\textrm{max}$ (\%)');
% yyaxis left
y1 = ylabel('$10*\log_{10}\left(\textrm{BER}\right)$');
% yyaxis right
% y2 = ylabel(' ');

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

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([-10 0]);
ax.YTick = [-10:2:0];

grid on;

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'r';
% ylim([0 1.1]);
% ax.YTick = [ ];

% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]

% print(FigName,'-dmeta');
% print(FigName,'-dpdf');
% Command =['pdfcrop ' FigName '.pdf ' FigName '.pdf'];
% system(Command);

%%
% -------------------------------------------------------------------------
% EVM vs SNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_evm_vs_snr_qam64';
hfig = figure('Name',fig_name);
% yyaxis left
plot(esn0_db,100*evm_max_da,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
% hold on;
% semilogy(EsN0_dB,Ps/log2(M),'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');



hleg = legend('64-QAM');
set(hleg,'Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('$E_s/N_0$ (dB)');
% yyaxis left
y1 = ylabel('$\textrm{EVM}_\textrm{max}$ (\%)');
% yyaxis right
% y2 = ylabel(' ');

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

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([0 45]);
ax.YTick = [0:10:50];

grid on;

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'r';
% ylim([0 1.1]);
% ax.YTick = [ ];

% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]

% print(FigName,'-dmeta');
% print(FigName,'-dpdf');
% Command =['pdfcrop ' FigName '.pdf ' FigName '.pdf'];
% system(Command);


%%
% -------------------------------------------------------------------------
% EVM vs OSNR
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_evm_vs_osnr_qam64';
hfig = figure('Name',fig_name);
% yyaxis left
plot(osnr_db_20gbd,100*evm_max_da,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
% hold on;
% semilogy(EsN0_dB,Ps/log2(M),'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');



hleg = legend('64-QAM');
set(hleg,'Location','SouthEast','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);

x1 = xlabel('OSNR (dB)');
% yyaxis left
y1 = ylabel('$\textrm{EVM}_\textrm{max}$ (\%)');
% yyaxis right
% y2 = ylabel(' ');

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

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([0 45]);
ax.YTick = [0:10:50];

grid on;

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'r';
% ylim([0 1.1]);
% ax.YTick = [ ];

% fig.PaperUnits = 'centimeters';
% fig.PaperType = 'a4';
% fig.PaperOrientation = 'landscape';
% fig.PaperPosition = [5 5 17 8.5]

% print(FigName,'-dmeta');
% print(FigName,'-dpdf');
% Command =['pdfcrop ' FigName '.pdf ' FigName '.pdf'];
% system(Command);



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

