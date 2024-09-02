% -------------------------------------------------------------------------
% Test of jones_prod_mv function
%
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
% 2024-07-31
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
% Switches
% -------------------------------------------------------------------------
do_debug = 1;
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;
fig.interpreter = 'latex';


% -------------------------------------------------------------------------
% Test 1: time-dependent matrix and vector
% -------------------------------------------------------------------------
A = [1 1 2; 
     3 2 1];
% Time-dependent vector.

J(1,1,:) = [2 1 5];
J(1,2,:) = [5 5 3];
J(2,1,:) = [3 5 1];
J(2,2,:) = [5 3 5];
% Time-dependent matrix.

actB = jones_prod_mv(J,A);
% Product of vector A by matrix J.

expB = [17 11 13;
        18 11 7];
% Expected result    

err = sum((actB - expB).^2,'all')
% Verify that the expected and actual results are identical.


% -------------------------------------------------------------------------
% Test 2: time-dependent vector, but constant matrix
% -------------------------------------------------------------------------
J1 = J(:,:,1);
% In this case the J matrix is a 2x2 matrix. 
% We sample it from the J matrix.

actB1 = jones_prod_mv(J1,A);
% Actual result

expB1 = [17 12 9;
         18 13 11];
% Expected result

err1 = sum((actB1 - expB1).^2,'all')
% Difference between expected and actual results.


% -------------------------------------------------------------------------
% Test 3: we swap the polarization components
% -------------------------------------------------------------------------
Js = [0 1; 1 0];
% Static Jones matrix

actBs = jones_prod_mv(Js,A);
% Actual result

expBs = [3 2 1; 1 1 2];
% Expected result

errs = sum((actBs - expBs).^2,'all')
% Difference between expected and actual results.




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

