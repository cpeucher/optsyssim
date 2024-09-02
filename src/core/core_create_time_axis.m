function [time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate)
% Create time and frequency axes
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function creates time and frequency axes made of regularly spaced 
% samples. It also returns the frequency spacing between samples in the 
% time and frequency domains.
% This function is adapted to the simulation of digital communication
% systems where the generated signals will consist of an integer number of
% symbols.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% nsamples_per_symbol = 32;
% nsymbols = 128;
% symbol_rate = 25e9;
% [time_array,dt,frequency_array,df] = core_create_time_axis(nsamples_per_symbol,nsymbols,symbol_rate);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% nsamples_per_symbol   number of samples per symbol [integer scalar]
%                           A power of 2 is higly recommended.
%
% nsymbols              number of symbols [integer scalar]
%                           A power of 2 is highly recommended.
%
% symbol_rate           symbol rate in baud (aka symbols per second)
%                           [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% time_array            time sample values, in s [real vector]
%
% dt                    time interval between 2 consecutive samples, in s
%                           [real scalar]
%
% frequency_array       relative frequency sample values, in Hz
%                           [real vector]
%
% df                    frequency interval between 2 consecutive samples, 
%                           in Hz [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = nsamples_per_symbol*nsymbols;
% Total number of samples. Should be a power of 2 if both the number of
% samples per symbol and the number of symbols are power of 2, as
% recommanded.
dt = 1/(nsamples_per_symbol*symbol_rate);
% Time interval between two consecutive samples.
df = symbol_rate/nsymbols;
% Frequency interval between two consecutive samples.
time_array = (0:(nsamples - 1))*dt;
% Build vector containing all the time samples.
frequency_array = (-nsamples/2:nsamples/2 - 1)*df;
% Build vector containing all the frequency samples.

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------