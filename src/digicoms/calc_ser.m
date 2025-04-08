function ser = calc_ser(symbs_cx,symbs_ref)
% Calculation of symbol-error-ratio (SER)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns the symbol-error-ratio of received symbols by 
% comparison with a vector of reference symbols.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% ser = calc_ser(symbs_cx,symbs_ref); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs_cx          complex symbols after decision [complex vector]
%
% symbs_ref         reference transmitted symbols [complex vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% ser               symbol-error-ratio [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if length(symbs_cx) ~= length(symbs_ref)    
    error('calc_ser: the received and reference complex symbols vectors do not have the same lengths.');    
end

nsymb = length(symbs_ref);
% Number of symbols

vones = ones(1,nsymb);
ser = sum(vones(abs(symbs_cx - symbs_ref) ~= 0))/nsymb;
% Calculate symbol error ratio


end