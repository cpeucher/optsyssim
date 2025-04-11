function [bits,symbs_diff] = diffdec_qpsk(symbs,decision)
% Differential decoding for QPSK
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements differential decoding of differentially-encoded
% QPSK symbols.
% This function is directly compatible with the differential encoding
% function:
% symbs = diffenc_qpsk(words_dec)
% Two approaches are implemented:
% 1. Initial quadrant decision of the received symbols followed by
% determination of the differential phase, or rather, equivalently of the
% 2-bit words corresponding to the differentially-encoded symbols. In this
% case, the determination of the 2-bit word is based on symbols that belong
% to the set {1+j, -1+j, -1-j, 1-j}, or equivalently the differential
% phases belong to the set {0, pi/2, pi, 3*pi/2}. No other values are
% possible since a hard quadrant decision is performed as the first step.
% This approach is referred to as 'hard' decision.
% 2. No hard decision is performed as an initial step. The differential
% phases / corresponding 2-word bits are determined directly from the
% received samples. In this case the differential phase is not limited to
% values in the set {0, pi/2, pi, 3*pi/2} when the input symbols are noisy,
% or, in other terms, we do not work on the constellation {1+j, -1+j, -1-j,
% 1-j} to determine the differential phase / 2-bit words, but on a noisy
% constellation consisting in a cloud of symbols around those points. 
% This approach is referred to as 'soft' decision, even though a hard
% decision is ultimately made to determine the 2-bit symbols. Nevertheless,
% this approach could allow the use of soft-decision codes, hence its name.
% The two proposed approaches correspond to the ones whose performance are
% represented in Fig. 1 in the paper:
% M. Kuschnerov et al., “Low complexity soft differential decoding of QPSK
% for forward error correction in coherent optic receivers,” in European 
% Conference on Optical Communication (ECOC), Torino, Italy, Sep. 2010,
% paper Th.9.A.6. doi: 10.1109/ECOC.2010.5621498.
% Approach 1 "hard" decision => binary DQPSK
% Approach 2 "soft" decision => soft DQPSK
% See the scripts:
% test_dsp_qpsk_diffenc.m for an illustration of the difference
% between the two approaches
% test_ber_vs_snr_awgn_qpsk_diffenc.m for a calculation of the BER vs Es/N0
% curves and comparison with non differentially-encoded QPSK (i.e. we
% reproduce Fig. 1 in the above reference by Kuschnerov et al.).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% nsymbols = 2^10;
% m = 4;
% esn0_db = 20;
% bits = generate_binary(nsymbols,m);
% [words_dec, words_bin] = conv_bin2dec(bits,log2(m));
% symbs = diffenc_qpsk(words_dec);
% symbs_rx = add_awgn(symbs,esn0_db);
% [bits_rx, symbs_diff] = diffdec_qpsk(symbs_rx,'hard');
% ber = calc_ber(bits_rx,bits(3:end));
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             input differentially encoded symbols [complex vector]
%
%                       It is not necessary for these symbols to be
%                       normalized so that their energy corresponds to that
%                       of the standard constellation for QPSK (i.e. 2).
%                       However, this is necessary if ones wants the 
%                       differential constellation to fall on the standard
%                       constellation (looks nicer...)
%
% decision          decision method [string]
%
%                       decision = 'hard'       
%                       decision = 'soft'
%
%                       See above for a description of the two decision
%                       approaches.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% bits               recovered bits [integer vector]
%
% symbs_diff         recovered differential constellation [complex vector]
%
%                       For representation or further processing, for
%                       instant using a soft-decision code.
%                       The differential constellation is
%                       c[k] = r[k] r*[k-1] exp(j pi/4) sqrt(2) where r[k]
%                       are the received symbols in case of 'soft' 
%                       decision, or
%                       c[k] = d[k] d*[k-1] exp(j pi/4) sqrt(2) where 
%                       d[k] are the hard-decided symbols in case of 'hard' 
%                       decision. 
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(decision,'hard')   
    
    symbs = sign(real(symbs)) + 1i*sign(imag(symbs));
    % Hard decision on the received symbols
    % Symbols after hard-decision d[k]
    
end

symbs_diff = symbs.*dsp_delay(conj(symbs),1)*exp(1i*pi/4)/sqrt(2);
% Differential decoding
% Those are the symbols:
% c[k] = r[k] r*[k-1] exp(j pi/4) sqrt(2) where r[k] are the received
% symbols in case of 'soft' decision, or
% c[k] = d[k] d*[k-1] exp(j pi/4) sqrt(2) where d[k] are the hard-decided
% symbols in case of 'hard' decision. 
% We referred to c[k] as the differential constellation, since 
% c[k] = sqrt(2) exp(j pi/4) exp(j Delta  Phi[k]) where Delta  Phi[k] is
% the differential phase.
% The product of the symbols by the 1-symbol delayed symbols leads to phase
% shifts of 0, pi/2, pi, 3*pi/2 (hard decision) or around (soft decision,
% in the presence of noise) and to a differentially-decoded constellation
% of sqrt(2)*{1,j,-1,-j} (hard decision) or clouds of points around (soft
% decision, in the presence of noise). The constellation is then rotated by
% pi/4 so that it is {1+j, -1+j, -1-j, 1-j} (hard) or around (soft).
% The sqrt(2) factor ensures we are on a standard constellation (if the
% input symbols have been properly normalised so that their energy is that
% of the standard constellation). This is not important for the
% determination of the binary data (msb,lsb), but looks nice if one wants
% to represent the differential constellation c[k].

symbs_diff = symbs_diff(2:end);
% Remove the first symbol, which is 0 due to the initial state of the
% delay

msb = imag(symbs_diff) < 0;
lsb = real(symbs_diff) < 0;
% Hard decision for the msb and lsb. This corresponds to the differential
% encoding implemented in symbs = diffenc_qpsk(words_dec).

bits = reshape([msb;lsb],[1,2*length(msb)]);
% Interleave the msb and lsb in a single bit stream

end