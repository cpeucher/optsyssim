function [bits,symbs_diff] = diffdec_qam(symbs,m)
% Differential decoding of square m-QAM constellation.
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements differential decoding of square m-QAM
% constellations according to the method described in 
% J.-K. Hwang, et al., Angle differential-QAM scheme for resolving phase 
% ambiguity in continuous transmission system, International Journal of 
% Communication Systems, vol. 21, no. 6, pp. 631-641, Jun. 2008, 
% doi: 10.1002/dac.914.
% This is an extension of the method implemented for QPSK modulation in our
% diffdec_qpsk.m function.
% Decoding is performed iteratively on the 0.5*log2(m) auxiliary
% subconstellations that are used to build square m-QAM constellations in
% the diffenc_qam.m function.
% For the time being only hard decision is implemented.
% See diffdec_qpsk.m for more details.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% nsymbols = 2^20;
% m = 64;
% esn0_db = 20;
% bits = generate_binary(nsymbols,m);
% symbs = diffenc_qam(bits,m);
% symbs_rx = add_awgn(symbs,esn0_db);
% [bits_rx,symbs_diff] = diffdec_qam(symbs_rx,m);
% ber = calc_ber(bits_rx,bits(log2(m)+1:end)) 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             input differentially encoded symbols [complex vector]
%
%                       These symbols are normalized so that their energy 
%                       corresponds to that of the standard constellation.
%
% m                 order of the modulation format [integer scalar]
%
%                       To generate m-QAM constellation.
%                       Only square constellations are generated with this
%                       function, i.e. m = 4, 16, 64...
%
% decision          decision method [string]
%
%                       decision = 'hard'       
%                       decision = 'soft'
%                       [not currently implemented]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% bits               recovered bits [integer vector]
%
% symbs_diff         recovered differential constellation [complex vector]
%
%                       For representation or further processing.
%                       Seer dsp_qpsk_diffdec.m for details.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


nsymbols = length(symbs);
% Number of symbols

nstages = 0.5*log2(m);
% Number of required decoding stages
% nstages = 1 for QPSK
% nstages = 2 for 16-QAM
% nstages = 3 for 64-QAM
% etc...

R = sqrt(2)*2.^(0:nstages - 1);
% Radii of the auxiliary sub-constellations
% R(nstages)        corresponds to the radius of the auxiliary 
%                   constellation defining the quadrant
%                   i.e. R(nstages) = sqrt(2)   for QPSK
%                        R(nstages) = 2*sqrt(2) for 16-QAM
%                        R(nstages) = 4*sqrt(2) for 64-QAM
% R(1)              corresponds to the radius of the smallest auxiliary
%                   constellation
%                        R(1) = sqrt(2)   for all formats

cp = zeros(nstages+1,nsymbols);
% Preallocate matrix that will be used to save the recovered auxiliary
% constellations
% We let an initial stage be present (stage 0) for ease of implementation
% of the algorithm
dcp = zeros(nstages,nsymbols);
% Corresponding differential constellations
% No need for stage 0

bits = zeros(2*nstages,nsymbols);
% Preallocate vector where decoded bits will be stored



for istage = 1:nstages
    % Loop over number of differential decoding stages.
    
    symbs = symbs - cp(istage,:);
    % Subtract the centers of the previoulsy calculated auxiliary
    % constellation from the symbols. This is done iteratively for each
    % stage.
    % For iteration 2, it should be ensured that the input symbs are
    % normalised to the relevant standard constellation.
    
    cp(istage+1,:) = R(nstages - istage + 1)/sqrt(2)*(sign(real(symbs)) + 1i*sign(imag(symbs)));
    % Hard decision on auxiliary constellation
    
    dcp(istage,:) = cp(istage + 1,:).*dsp_delay(conj(cp(istage + 1,:)),1)*exp(1i*pi/4)/R(nstages - istage + 1);
    % Differential constellation
    
    bits(1 + 2*(istage - 1),:) = imag(dcp(istage,:)) < 0;
    bits(2 + 2*(istage - 1),:) = real(dcp(istage,:)) < 0;
    % Decision on the lsb and msb for the current differential
    % constellation
    
    
    
end
% End of loop over differential decoding stages


bits = reshape(bits,[1 nsymbols*log2(m)]);
% We combine the msb and lsb recovered after each decoding stage

bits = bits(log2(m)+1:end);
% The log2(m) bits corresponding to the first symbol are not meaningful

symbs_diff = sum(dcp,1);
% We reconstruct the differential constellation from the differential
% auxiliary constellations

symbs_diff = symbs_diff(2:end);
% We remove the first symbol, which is 0 due to the use of the dsp_delay
% function in the calculation of each auxiliary differential constellation    

end