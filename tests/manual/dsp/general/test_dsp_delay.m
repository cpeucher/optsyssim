% -------------------------------------------------------------------------
% Test of dsp_delay. m function
%
% 2021-06-23
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
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;


% -------------------------------------------------------------------------
% Input sequence is line vector
% -------------------------------------------------------------------------
a = linspace(1,10,10)
% Initial sequence.

b = dsp_delay(a,5)
% Apply 5 sample delay with default initial conditions. 

c = dsp_delay(a,5,[1 2 3 4 5])
% Apply 5 sample delay with specified initial conditions.

% We also compare with the outcome of the dsp_fir_linear function
z = [0 -1 -2 -3 -4 -5];
% Initial FIR buffer.
b = [0 0 0 0 0 1];
% Tap coefficients

d1 = dsp_fir_linear(a,b,z)
% Apply FIR filter.
d2 = dsp_delay(a,5,fliplr(z(1:end-1)))
% Apply delay.


% -------------------------------------------------------------------------
% Parallel processing: input is matrix
% -------------------------------------------------------------------------
a = linspace(1,50,50);
a = reshape(a,5,10)
% Create a matrix of 5 sequences to delay.

b = dsp_delay(a,3)
% Delay each sequence by 3 samples with default initial conditions.

c = dsp_delay(a,3,[9 10 11])
% Delay each sequence by 3 samples, with identical specified initial
% conditions,

d = dsp_delay(a,3,[9 10 11;12 13 14;15 16 17;18 19 20;21 22 334])
% Delay each sequence by 3 samples, with different specified initial
% conditions for each sequence.

e = dsp_delay(a,3,[NaN NaN NaN])
% That also works if the initial state is set to NaN...



% -------------------------------------------------------------------------
% Save data
% -------------------------------------------------------------------------
% file_name_data = [file_name_core_data '_' time_stamp]
% save(file_name_data,'');