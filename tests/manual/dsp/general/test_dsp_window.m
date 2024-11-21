% -------------------------------------------------------------------------
% Test of dsp_window function
% 
% We generate different standard windows (Hann, Kaiser, Hamming,
% Bartlett...) and compare them as well as their spectra.
%
% 2021-12-11
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

% -------------------------------------------------------------------------
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
fig_interpreter = 'latex';




% -------------------------------------------------------------------------
% Check defintions of the windows (periodic / symmetric)
% -------------------------------------------------------------------------
fft_length = 2^16;
% Number of points.

nwindow = 52;
% Windows size.



% We generate a Hann window using the Oppenheim & Schafer definition
% A.V. Oppenheim, R.W. Schafer, and J.R. Buck, Discrete-Time Signal 
% Processing, 2nd ed (Prentice Hall, Upper Saddle River, N.J, 1999).
w = 0.5 - 0.5*cos(2*pi*(0:nwindow)/nwindow);
% Generates a window of size M+1 with w[1] = 0 and w[M+1] = 0.

wt = w(1:nwindow);
% Truncate window to M points.

% We now use our old DSP_Window function.
params_window.type = 'hann';
params_window.length = nwindow;
params_window.symmetry = 'periodic';%'Symmetric';
params_window.beta = 2.4*pi;% For Kaiser window.
wf = dsp_window(params_window);

% Calculate spectrum with N-point FFT
W = fft(w,fft_length);
Ws = fft(wt,fft_length);
Wf = fft(wf,fft_length);



omega = (0:1/fft_length:(fft_length-1)/fft_length);
% Normalised angular frequency


figure('Name',['Hann window of size M = ' num2str(nwindow)])
subplot(1,2,1)
stem(w,'ro');
hold on
stem(wt,'bs')
stem(wf,'^g')
xlabel('sample')
ylabel('amplitude');
xlim([0.5 nwindow+1+0.5]);
legend('M+1 points (O&S)','M points (Matlab periodic)','from DSP function')
subplot(1,2,2)
plot(omega(1:fft_length/2),10*log10(abs(W(1:fft_length/2)).^2/max(abs(W(1:fft_length/2)).^2)),'r');
hold on;
plot(omega(1:fft_length/2),10*log10(abs(Ws(1:fft_length/2)).^2/max(abs(W(1:fft_length/2)).^2)),'b--');
xlabel('$\hat{\omega}/2\pi$','Interpreter','latex');
ylabel('magnitude (dB)');
xlim([0 1/2]);
ylim([-80 0]);



%%
% -------------------------------------------------------------------------
% Compare different windows
% -------------------------------------------------------------------------
fft_length = 2^16;
% FFT length for window spectrum calculation

params_window.length = 51;
params_window.symmetry = 'symmetric';%'Periodic';
% General parameters for all windows.
% Here we choose 'Symmetric' to ensure all windows can be directly
% compared.


indices = [1:params_window.length];
% Vector of indices for representation.


% Rectangular window
params_window.type = 'rectangular';
w_rectangular = dsp_window(params_window);
[neb_rectangular,mag_spec_rectangular,norm_ang_freq_rectangular] = dsp_window_properties(w_rectangular,fft_length);

fprintf('\n\n%s\t%s%3.2f%s\n',[params_window.type ':'],'NEB= ',neb_rectangular,' frequency bins');


% Hann window
params_window.type = 'hann';
w_hann = dsp_window(params_window);
[neb_hann,mag_spec_hann,norm_ang_freq_hann] = dsp_window_properties(w_hann,fft_length);

fprintf('%s\t\t%s%3.2f%s\n',[params_window.type ':'],'NEB= ',neb_hann,' frequency bins');



% Kaiser window
params_window.type = 'kaiser';
params_window.Beta = 2.4*pi;% For Kaiser window.
w_kaiser = dsp_window(params_window);
[neb_kaiser,mag_spec_kaiser,norm_ang_freq_kaiser] = dsp_window_properties(w_kaiser,fft_length);

fprintf('%s\t\t%s%3.2f%s\n',[params_window.type ':'],'NEB= ',neb_kaiser,' frequency bins');

% Hamming window
params_window.type = 'hamming';
w_hamming = dsp_window(params_window);
[neb_hamming,mag_spec_hamming,norm_ang_freq_hamming] = dsp_window_properties(w_hamming,fft_length);

fprintf('%s\t%s%3.2f%s\n',[params_window.type ':'],'NEB= ',neb_hamming,' frequency bins');

% Bartlett window
params_window.type = 'bartlett';
w_bartlett = dsp_window(params_window);
[neb_bartlett,mag_spec_bartlett,norm_ang_freq_bartlett] = dsp_window_properties(w_bartlett,fft_length);

fprintf('%s\t%s%3.2f%s\n',[params_window.type ':'],'NEB= ',neb_bartlett,' frequency bins');







% -------------------------------------------------------------------------
% Plot window 
% -------------------------------------------------------------------------
FIGPARAM.FontName = 'arial';
FIGPARAM.FontSize = 15;
FIGPARAM.Interpreter = 'latex';


FigName = ['Windows of size M = ' num2str(nwindow)];
fig = figure('Name',FigName);
% yyaxis left
plot(indices,w_rectangular,'LineWidth',1.5,'Color','k','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
plot(indices,w_bartlett,'LineWidth',1.5,'Color','g','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(indices,w_hann,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(indices,w_hamming,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(indices,w_kaiser,'LineWidth',1.5,'Color','c','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');


hleg = legend('rectangular','Bartlett','Hann','Hamming','Kaiser');
set(hleg,'Location','South','Box','off','Interpreter',FIGPARAM.Interpreter,'FontSize',FIGPARAM.FontSize,'NumColumns',1);

x1 = xlabel('sample');
% yyaxis left
y1 = ylabel('amplitude');
% yyaxis right
% y2 = ylabel(' ');

x1.Interpreter = FIGPARAM.Interpreter;
y1.Interpreter = FIGPARAM.Interpreter;
y2.Interpreter = FIGPARAM.Interpreter;


ax = gca;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = FIGPARAM.Interpreter;

ax.XAxis.FontName = FIGPARAM.FontName;
ax.XAxis.FontSize = FIGPARAM.FontSize;
xlim([0 nwindow]);
% ax.XTick = [ ];

% yyaxis left
ax.YAxis(1).FontName = FIGPARAM.FontName;
ax.YAxis(1).FontSize = FIGPARAM.FontSize;
% ax.YAxis(1).Color = 'k';
ylim([0 1.1]);
% ax.YTick = [ ];

% yyaxis right
% ax.YAxis(2).FontName = FIGPARAM.FontName;
% ax.YAxis(2).FontSize = FIGPARAM.FontSize;
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
% Plot spectrum 
% -------------------------------------------------------------------------
FigName = ['Spectra of windows of size M = ' num2str(nwindow)];
fig = figure('Name',FigName);
% yyaxis left
plot(norm_ang_freq_rectangular/2/pi,mag_spec_rectangular,'LineWidth',1.5,'Color','k','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
plot(norm_ang_freq_bartlett/2/pi,mag_spec_bartlett,'LineWidth',1.5,'Color','g','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(norm_ang_freq_hann/2/pi,mag_spec_hann,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(norm_ang_freq_hamming/2/pi,mag_spec_hamming,'LineWidth',1.5,'Color','r','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
plot(norm_ang_freq_kaiser/2/pi,mag_spec_kaiser,'LineWidth',1.5,'Color','c','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');


hleg = legend('rectangular','Bartlett','Hann','Hamming','Kaiser');
set(hleg,'Location','NorthEast','Box','off','Interpreter',FIGPARAM.Interpreter,'FontSize',FIGPARAM.FontSize,'NumColumns',1);

x1 = xlabel('$\hat{\omega}/2\pi$ (radians/sample)');
% yyaxis left
y1 = ylabel('magnitude (dB)');
% yyaxis right
% y2 = ylabel(' ');

x1.Interpreter = FIGPARAM.Interpreter;
y1.Interpreter = FIGPARAM.Interpreter;
y2.Interpreter = FIGPARAM.Interpreter;


ax = gca;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = FIGPARAM.Interpreter;

ax.XAxis.FontName = FIGPARAM.FontName;
ax.XAxis.FontSize = FIGPARAM.FontSize;
xlim([0 0.5]);
% ax.XTick = [ ];

% yyaxis left
ax.YAxis(1).FontName = FIGPARAM.FontName;
ax.YAxis(1).FontSize = FIGPARAM.FontSize;
% ax.YAxis(1).Color = 'k';
ylim([-80 0]);
% ax.YTick = [ ];

% yyaxis right
% ax.YAxis(2).FontName = FIGPARAM.FontName;
% ax.YAxis(2).FontSize = FIGPARAM.FontSize;
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
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

