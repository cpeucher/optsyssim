function samps = add_cfo(samps,cfo_normalised)
% Add carrier frequency offset (CFO) to the signal samples
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function adds carrier frequency offset (CFO) to the signal samples.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% cfo_normalised = cfo_absolute/symbol_rate;
% symbs = add_cfo(symbs,cfo_normalised);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% samps             signal samples [complex vector]
%
% cfo_normalised    normalised carrier frequency offset [real scalar] 
%                       or [real vector]
%
%                       cfo_normalised is defined as the product of the
%                       carrier frequency offset (in Hz) and the
%                       sampling interval of the signal (in s)
%                       cfo_normalised = Df * Tsa
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% samps             signal samples with added CFO [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% samps = samps.*exp(1i*2*pi*cfo_normalised*[0:1:length(samps) - 1]);

samps = samps.*exp(1i*2*pi*cfo_normalised.*[1:1:length(samps)]);

end
