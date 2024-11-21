% -------------------------------------------------------------------------
% Some tests of dsp_fir_linear function...
%
% 2021-06-24
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
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;


% -------------------------------------------------------------------------
% Test of moving average 
% -------------------------------------------------------------------------
n = [1:1:100];
% Indices.

x = cos(n/7);
% Sequence.

L = 5;
% Number of taps.

y = dsp_fir_linear(x,ones(1,L),zeros(1,L))/L;
% FIR filter aplied to input sequence. Moving average.

yc = circshift(y,-0.5*(L-1),2);
% Non-causal compensation of the group delay.

xc = dsp_delay(x,0.5*(L-1));
% Delayed original sequence.

figure('Name','Moving average using FIR filter - test on sinusoid')
plot(n,x,'b')
hold on
plot(n,y,'r')
plot(n,yc,'k--')
plot(n,xc,'g')
legend('original sequence','filtered sequence','delay-compensated (non causal) filtered sequence','delayed original sequence','Location','SouthEast')

% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');

