function ber = dsp_ber(bits,bits_ref)
% Bit-error-ratio calculation by comparison of two binary vectors
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns the bit-error-ratio of received binary data by 
% comparison with a vector of reference binary data.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% ber = calc_ber(bits,bits_ref);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% bits              received bit stream [logical vector]
%                       
% bits_ref          reference bit stream [logical vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% ber               calculated bit-error-ratio [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if length(bits) ~= length(bits_ref)
    error('dsp_ber: the received and reference binary bit vectors do not have the same lengths.');
end

ber = sum(xor(bits,bits_ref))/length(bits);
% That's it

end