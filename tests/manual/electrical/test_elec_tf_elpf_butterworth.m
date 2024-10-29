% -------------------------------------------------------------------------
% Test of low-pass Butterworth filter implementation
%
% We check the coefficients of the Butterworth polynomials used in the
% denominator of the transfer function and calculate selected transfer
% functions for orders n = 1 to 5.
% We also check the calculation of the group delay for a selected order.
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
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
fig_interpreter = 'latex';


% -------------------------------------------------------------------------
% First we check the calculation of the Butterwoth polynomials as employed
% in the elec_tf_elpf.m function
% -------------------------------------------------------------------------
for iorder = 1:10
    
    p = zeros(1,iorder + 1);
    % The coefficients of the Butterwoth polynomial will be stored in p.
    
    p(1) = 1;
    % Term of degree 0
    
    for ii = 2:iorder+1
        
        p(ii) = cos((ii - 2)*pi/2/iorder)/sin((ii - 1)*pi/2/iorder)*p(ii - 1);    
        % Term of degree ii - 1
        
    end
    
    butterworth_polynomials{iorder,:} = fliplr(p);
    % Matlab polynomial arrays go from coefficient of degree n to
    % coefficient of degree 0.
    
end

% Conclusion:
% OK, we retrieve the polynomials found in tables in the literature.

% -------------------------------------------------------------------------
% Next we calculate the transfer function (magnitude and phase) for various
% filter orders.
% -------------------------------------------------------------------------
params_elpf.type = 'butterworth';
params_elpf.f3dB = 1;
% Cutt-off frequency
order_max = 5;
% Maximum filter order that will be considered
npoints = 1000;
% Number of points at which the transfer function will be evaluated

freq = logspace(-3,3,npoints);
% Frequency axis

tf = zeros(order_max,npoints);
% Pre-initialize data

for iorder = 1:5
    % Calculate the transfer functions for the various orders
    
    params_elpf.order = iorder;
    
    tf(iorder,:) = elec_tf_elpf(params_elpf,freq);
    
end


fig_name = [file_name_core_figure '_amplitude_response'];
figure('Name',fig_name)
semilogx(freq,10*log10(abs(tf(1,:)).^2),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
hold on
semilogx(freq,10*log10(abs(tf(2,:)).^2),'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,10*log10(abs(tf(3,:)).^2),'Color','g','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,10*log10(abs(tf(4,:)).^2),'Color','c','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,10*log10(abs(tf(5,:)).^2),'Color','k','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
xlabel('frequency $f$ (Hz)','Interpreter',fig_interpreter)
ylabel('$\left|H\left(f\right)\right|^2$ (dB)','Interpreter',fig_interpreter)
legend('n = 1','n = 2','n = 3','n = 4','n = 5','Location','SouthWest','Box','on','Interpreter',fig_interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig_interpreter)
ax = gca;
ax.TickLabelInterpreter = fig_interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([min(freq) max(freq)])
ylim([-65 5])
grid on
if do_print print(fig_name,'-dmeta'); end


fig_name = [file_name_core_figure '_phase_response'];
figure('Name',fig_name)
semilogx(freq,unwrap(angle(tf(1,:))),'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
hold on
semilogx(freq,unwrap(angle(tf(2,:))),'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,unwrap(angle(tf(3,:))),'Color','g','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,unwrap(angle(tf(4,:))),'Color','c','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
semilogx(freq,unwrap(angle(tf(5,:))),'Color','k','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
xlabel('frequency $f$ (Hz)','Interpreter',fig_interpreter)
ylabel('$\phi\left(f\right)$ (rad)','Interpreter',fig_interpreter)
legend('n = 1','n = 2','n = 3','n = 4','n = 5','Location','SouthWest','Box','on','Interpreter',fig_interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig_interpreter)
ax = gca;
ax.TickLabelInterpreter = fig_interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
xlim([min(freq) max(freq)])
ylim([-10 1])
grid on
if do_print print(fig_name,'-dmeta'); end



%%
% -------------------------------------------------------------------------
% Calculation of the group delay for specified order
% -------------------------------------------------------------------------
freq2 = linspace(0,5,npoints);
% A linear frequency grid will make group delay calculations better...
params_elpf.order = 3;
% Filter order
    
tf2 = elec_tf_elpf(params_elpf,freq2);
% Transfer function calculation
phase2 = unwrap(angle(tf2));
% The phase is defined here as the argument of the complex transfer
% function

[freq2_delay,group_delay2] = num_diff(1,freq2,phase2);
group_delay2 = -group_delay2/2/pi;
% With our phase definition, the group delay is minus the derivative of the
% phase with respect to angular frequency

phase_delay2 =  -phase2./freq2/2/pi;
% We can also calculate the phase delay. The minus sign makes it compatible
% with the above definition of the phase


fig_name = [file_name_core_figure '_group_delay'];
figure('Name',fig_name)
yyaxis left
plot(freq2,abs(tf2).^2,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
ylabel('magnitude $\left|H\left(f\right)\right|^2$ (dB)','Interpreter',fig_interpreter)
yyaxis right
plot(freq2_delay,group_delay2,'Color','r','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12)
ylabel('group delay $\tau_g$ (s)','Interpreter',fig_interpreter)
xlabel('frequency $f$ (Hz)','Interpreter',fig_interpreter)
legend('magnitude','delay','Location','NorthEast','Box','on','Interpreter',fig_interpreter,'FontSize',15)
% title('','FontWeight','Normal','Interpreter',fig_interpreter)
ax = gca;
ax.TickLabelInterpreter = fig_interpreter;
ax.FontSize = 15;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color  = 'r';
xlim([min(freq2) max(freq2)])
%ylim([-10 1])
grid on
if do_print print(fig_name,'-dmeta'); end

[freq_tf,magnitude,phase,phase_delay,freq_delay,group_delay] = elec_tf_analysis(freq2,params_elpf,1000,1);
% Test elec_tf_analysis function

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');