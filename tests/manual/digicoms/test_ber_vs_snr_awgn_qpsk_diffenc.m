% -------------------------------------------------------------------------
% Calculation of BER vs Eb/N0 for differentially-encoded QPSK modulation
% with 'hard' and 'soft' differential decoding.
% 
% We reproduce Fig. 1 in
% M. Kuschnerov et al., “Low complexity soft differential decoding of QPSK
% for forward error correction in coherent optic receivers,” in European 
% Conference on Optical Communication (ECOC), Torino, Italy, Sep. 2010,
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
do_print = 0;

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
fig.interpreter = 'latex';
fig.font_size = 18;

mred = [178 0 24]/255;

% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));



% -------------------------------------------------------------------------
% Global parameters
% -------------------------------------------------------------------------
nsymbols = 2^25;
nsymbols = 2^20;
% Number of symbols. Should be power of 2.

% -------------------------------------------------------------------------
% Define constellations
% -------------------------------------------------------------------------
m = 4;
% Constellation order. M points with log2(M) even (square constellation)

[constellation_nat,~,~] = define_constellation('qpsk_natural',m);
% Define the constellation for natural mapping
[constellation_gray,~,~] = define_constellation('qpsk_gray',m);
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
symbs_nat = mapping(words_dec,constellation_nat);
symbs_gray = mapping(words_dec,constellation_gray);
% Mapping.

symbs_diff = diffenc_qpsk(words_dec);
% With differential encoding.


% -------------------------------------------------------------------------
% Save reference signals
% -------------------------------------------------------------------------
symbs_ref_nat = symbs_nat;
symbs_ref_gray = symbs_gray;

% -------------------------------------------------------------------------
% Eb/N0 range
% -------------------------------------------------------------------------
% esn0_db = 15;
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
    
    
    symbs_rx_gray = add_awgn(symbs_gray,esn0_db(isnr));
    symbs_rx_nat = add_awgn(symbs_nat,esn0_db(isnr));
    symbs_rx_diff = add_awgn(symbs_diff,esn0_db(isnr));
    % Add white Gaussian noise.  
    
    
    vones = ones(1,nsymbols);
    
    [symbs_cx,bits_rx,words_dec_rx] = decision_qpsk_hard(symbs_rx_gray,'qpsk_gray');
    % Decision.
    ber_gray(isnr) = calc_ber(bits_rx,data);
    % BER.
    ser_bits_gray(isnr) = sum(vones(abs(words_dec - words_dec_rx) ~= 0))/nsymbols;
    % SER after demapping.
    ser_symbs_gray(isnr) = sum(vones(abs(symbs_cx - symbs_ref_gray) ~= 0))/nsymbols;
    % SER without demapping.
    
    [symbs_cx,bits_rx,words_dec_rx] = decision_qpsk_hard(symbs_rx_nat,'qpsk_natural');
    % Decision.
    ber_nat(isnr) = calc_ber(bits_rx,data);
    % BER.
    ser_bits_nat(isnr) = sum(vones(abs(words_dec - words_dec_rx) ~= 0))/nsymbols;
    % SER after demapping.
    ser_symbs_nat(isnr) = sum(vones(abs(symbs_cx - symbs_ref_nat) ~= 0))/nsymbols;
    % SER without demapping.   
    
    
    [data_rx_hd,~] = diffdec_qpsk(symbs_rx_diff,'hard');
    % Hard differential decoding.    
    ber_diff_hd(isnr) = calc_ber(data_rx_hd,data(3:end));
    % BER calculation for hard differential decoding.
    
    [data_rx_sd,~] = diffdec_qpsk(symbs_rx_diff,'soft');
    % Soft differential decoding.    
    ber_diff_sd(isnr) = calc_ber(data_rx_sd,data(3:end));
    % BER calculation for soft differential decoding.

    

    

    
end




%%
% -------------------------------------------------------------------------
% Plot results: bit error probability / ratio
% -------------------------------------------------------------------------
fig_name = 'fig_dsp_ber_vs_snr_qpsk_diffenc';
figure('Name',fig_name)
semilogy(esn0_db - 10*log10(2),ber_gray,'LineWidth',1.5,'Color','k','LineStyle','-','Marker','s','MarkerSize',12,'MarkerFaceColor','k','HandleVisibility','on');
hold on;
% semilogy(esn0_db - 10*log10(2),ber_nat,'LineWidth',1.5,'Color','c','LineStyle','none','Marker','s','MarkerSize',12,'MarkerFaceColor','c','HandleVisibility','on');
% semilogy(esn0_db - 10*log10(2),pb_qpsk_gray,'LineWidth',1.5,'Color','k','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
% semilogy(esn0_db- 10*log10(2) ,pb_qpsk_nat,'LineWidth',1.5,'Color','c','LineStyle',':','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
semilogy(esn0_db - 10*log10(2),ber_diff_hd,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','o','MarkerSize',12,'MarkerFaceColor','b','HandleVisibility','on');
semilogy(esn0_db - 10*log10(2),ber_diff_sd,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','^','MarkerSize',12,'MarkerFaceColor','r','HandleVisibility','on');


xlabel('$E_b / N_0$ (dB)','Interpreter',fig.interpreter)
ylabel('BER','Interpreter',fig.interpreter)
% legend('QPSK - Gray','QPSK - Natural','QPSK - Gray - Theory', 'QPSK - Natural - Theory','DQPSK - hard','DQPSK - soft','Location','SouthWest','Box','on','Interpreter',fig_interpreter,'FontSize',15)
legend('QPSK - Gray','DQPSK - hard','DQPSK - soft','Location','SouthWest','Box','on','Interpreter',fig.interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig_interpreter)
xlim([3 9]);
ylim([1e-4 1e-1]);
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.XTick = [5:2:15];
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
grid on
% xlim([0 18])
% ylim([-0.5 0.05])
if do_print
    print(fig_name,'-dmeta');    
    print(fig_name,'-dpdf');
    Command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(Command);    
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

