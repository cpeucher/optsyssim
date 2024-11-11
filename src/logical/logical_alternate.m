function seq = logical_alternate(type,sequence_length)
% Generation of alternate binary sequence
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function generates a binary pattern made of alternate 1's and 0's.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% seq = logical_alternate(alternate0,nsymbols);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% type              type of the generated pattern [string]
%
%                       If type = 'alternate0' the pattern starts with a 0.
%                       If type = 'alternate1' the pattern starts with a 1.
%
% sequence_length   number of bits in the generated pattern [integer]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% seq               generated pattern [binary vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

seq = [1:sequence_length];

if strcmp(type,'alternate0')
    seq = ~mod(seq,2);   
elseif strcmp(type,'alternate1')
    seq = mod(seq,2);    
else
    disp('logical_alternate: sequence type not implemented.');
end

end