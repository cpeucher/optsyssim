% -------------------------------------------------------------------------
% Calculation of symbol-error ratio (SER) versus signal-to-noise ratio
% (SNR) and bit-error ratio versus SNR for quadrature-phase shift keying 
% (QPSK) modulation over the additive white gaussian noise channel (AWGN).
%
% Both natural and gray mapping are considered.
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


% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.font_name = 'arial';
fig.font_size = 15;
fig.interpreter = 'latex';

mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));




% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^25;
% Number of symbols. Should be power of 2.



% -------------------------------------------------------------------------
% Define constellation
% -------------------------------------------------------------------------
m = 4;
% Constellation order. M points with log2(M) even (square constellation)

[constellation_nat,~,~] = define_constellation('qpsk_natural',m);
% Define the constellation for natural mapping.

[constellation_gray,~,~] = define_constellation('qpsk_gray',m);
% Define the constellation for Gray mapping.



% -------------------------------------------------------------------------
% Binary data stream
% -------------------------------------------------------------------------
data = generate_binary(nsymbols,m);
% Generate binary data.


[words_dec,words_bin] = conv_bin2dec(data,log2(m));
% Convert binary data into 2-bit words.



% -------------------------------------------------------------------------
% Mapping
% -------------------------------------------------------------------------
sig_nat = constellation_nat(words_dec + 1);
sig_gray = constellation_gray(words_dec + 1);
% Mapping.

sig_ref_nat = sig_nat;
sig_ref_gray = sig_gray;
% Save reference signals.


% -------------------------------------------------------------------------
% Eb/N0 range
% -------------------------------------------------------------------------
esn0_db = [5:0.5:15];

% -------------------------------------------------------------------------
% Theoretical bit-error probabilities
% -------------------------------------------------------------------------
pb_qpsk_gray = 0.5*erfc(sqrt(10.^(esn0_db/10)/2));
ps_qpsk = 2*pb_qpsk_gray - pb_qpsk_gray.^2;
pb_qpsk_nat = 1.5*pb_qpsk_gray - pb_qpsk_gray.^2;


%--------------------------------------------------------------------------
% BER vs SNR
%--------------------------------------------------------------------------
for isnr = 1:length(esn0_db)
    % Loop over SNR.
    
    
    sig_rx_gray = add_awgn(sig_gray,esn0_db(isnr));
    sig_rx_nat = add_awgn(sig_nat,esn0_db(isnr));
    % Add white Gaussian noise.  
    
    
    vones = ones(1,nsymbols);
    
    [symbs_cx,bits_rx,words_dec_rx] = decision_qpsk_hard(sig_rx_gray,'qpsk_gray');
    % Decision.
    ber_gray(isnr) = calc_ber(bits_rx,data);
    % BER.
    ser_bits_gray(isnr) = sum(vones(abs(words_dec - words_dec_rx) ~= 0))/nsymbols;
    % SER after demapping.
    ser_symbs_gray(isnr) = sum(vones(abs(symbs_cx - sig_ref_gray) ~= 0))/nsymbols;
    % SER without demapping.
    
    [symbs_cx,bits_rx,words_dec_rx] = decision_qpsk_hard(sig_rx_nat,'qpsk_natural');
    % Decision.
    ber_nat(isnr) = calc_ber(bits_rx,data);
    % BER.
    ser_bits_nat(isnr) = sum(vones(abs(words_dec - words_dec_rx) ~= 0))/nsymbols;
    % SER after demapping.
    ser_symbs_nat(isnr) = sum(vones(abs(symbs_cx - sig_ref_nat) ~= 0))/nsymbols;
    % SER without demapping.    
    
end




%%
% -------------------------------------------------------------------------
% Plot results: symbol error probability / ratio
% -------------------------------------------------------------------------

fig_name = 'fig_dsp_ser_vs_snr_qpsk';
hfig1 = figure('Name',fig_name);
semilogy(esn0_db,ser_bits_gray,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
semilogy(esn0_db,ser_bits_nat,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
semilogy(esn0_db,ps_qpsk,'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hleg = legend('Gray','Natural','Theory');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);
x1 = xlabel('$E_s / N_0$ (dB)');
y1 = ylabel('SER / $P_\mathrm{se}$');
x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;
ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;
ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([4.5 15.5]);
ax.XTick = [5:2:15];

ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-7 1]);
% ax.YTick = [ ];

grid on

% print(fig_name,'-dmeta');
print(fig_name,'-dpdf');
Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
system(Command);


%%
% -------------------------------------------------------------------------
% Plot results: bit error probability / ratio
% -------------------------------------------------------------------------

fig_name = 'fig_dsp_ber_vs_snr_qpsk';
hfig2 = figure('Name',fig_name);
semilogy(esn0_db,ber_gray,'LineWidth',1.5,'Color','b','LineStyle','none','Marker','o','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on;
semilogy(esn0_db,ber_nat,'LineWidth',1.5,'Color','r','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
semilogy(esn0_db,pb_qpsk_gray,'LineWidth',1.5,'Color','b','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
semilogy(esn0_db,pb_qpsk_nat,'LineWidth',1.5,'Color','r','LineStyle',':','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hleg = legend('Gray','Natural','Theory (Gray)', 'Theory (Natural)');
set(hleg,'Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',fig.font_size,'NumColumns',1);
x1 = xlabel('$E_s / N_0$ (dB)');
y1 = ylabel('BER / $P_\mathrm{be}$');
x1.Interpreter = fig.interpreter;
y1.Interpreter = fig.interpreter;
y2.Interpreter = fig.interpreter;
ax = gca;
ax.LineWidth = 1;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;
ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
xlim([4.5 15.5]);
ax.XTick = [5:2:15];
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
% ax.YAxis(1).Color = 'k';
ylim([1e-7 1]);
% ax.YTick = [ ];

grid on;

% print(fig_name,'-dmeta');
print(fig_name,'-dpdf');
Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
system(Command);




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