function [symbs_cx,bits,words_dec] = decision_qpsk_hard(symbs_rx,mapping)
% Hard decision for QPSK
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function performs a hard decision on a QPSK constellation 
% that has already been normalised to the standard constellation, i.e. with
% points centered around (1+1i, -1+1i, -1-1i, 1-1i).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% es_rx = sum(abs(symbs_rx).^2)/length(symbs_rx);
% symbs_rx = symbs_rx*sqrt(2/es_rx);
% [symbs_cx,bits,words_dec] = decision_qpsk_hard(symbs_rx,'qpsk_gray');
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs_rx          received complex symbols [complex vector]
%
%                       They have already been normalised so that they are
%                       centered around the standard constellation points,
%                       i.e. (1+1i, -1+1i, -1-1i, 1-1i).
%
%                       If AWGN is the main degradation, this can be done
%                       by ensuring that the energy per symbol of the 
%                       received constellation assuming all symbols are
%                       equally represented) is equal to 2.
%
% mapping           mapping of the constellation [string]
%
%                       'qpsk_natural'
%                       'qpsk_gray'
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% symbs_cx          result of decision based on symbol location 
%                       [complex vector]
%
%                       Vector with elements in (1+1i, -1+1i, -1-1i, 1-1i).
%                       after quadrant decision.
%                       Here the actual mapping does not matter.
%
% bits              recovered bit stream after demapping [binary vector]
%
%                       Demapping is performed.
%
% words_dec         reconstituted symbols after demapping [integer vector]
%
%                       The symbols are expressed in decimal
%                       representation, i.e. 0, 1, 2 3.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

symbs_cx = sign(real(symbs_rx)) + 1i*sign(imag(symbs_rx));
% Symbol decision in the complex plane. The resulting vector contains
% points in the standard constellation that correspond to the result of
% the hard decision.


switch mapping
    
    case 'qpsk_natural'        
      
        msb = imag(symbs_rx) > 0;        
        lsb = msb .* (real(symbs_rx)<0) + not(msb) .* (real(symbs_rx)>0);        
        
    case 'qpsk_gray'
        
        msb = imag(symbs_rx) > 0;
        lsb = real(symbs_rx) > 0;
        
    otherwise
        
        error('decision_qpsk_hard: mapping not implemented.')        
        
end

        words_dec = 2*msb + lsb;
        % Decimal representation of the words after demapping
        bits =  logical(reshape(cat(1,msb(:)',lsb(:)'),1,2*length(msb)));
        % Recovered bit stream after demapping

end