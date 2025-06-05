function [cfo_estimate,estimator_delay] = cfo_diffphase(samps,mpower,sample_delay,block_length,varargin)
% CFO estimation for M-PSK by differential phase algorithm (Leven)
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the carrier frequency offset (CFO) estimation 
% method based on the differential phase algorithm described in:
% A. Leven, N. Kaneda, U.-V. Koc, and Y.-K. Chen, "Frequency estimation in 
% intradyne reception," IEEE Photonics Technology Letters, vol. 19, no. 6,
% pp. 366-368, Mar. 2007.
% http://dx.doi.org/10.1109/LPT.2007.891893
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tk = [0:1:length(symbs) - 1]/symbol_rate;
% mpower = 4;
% sample_delay = 1;
% block_length = 1000;
% [cfo_estimate, cfo_estimation_delay] = cfo_diffphase(symbs,mpower,sample_delay,block_length);
% symbs_comp = dsp_delay(symbs,cfo_estimation_delay).*exp(-1j*2*pi*cfo_estimate.*dsp_delay(tk,cfo_estimation_delay)*symbol_rate);
% symbs_comp = symbs_comp([2*cfo_estimation_delay + 1:end]);
% 
% nblocks = floor(length(symbs)/block_length);
% cfo_estimate = cfo_diffphase(symbs,mpower,sample_delay,block_length,'block');
% symbs_comp = symbs(1:nblocks*block_length).*exp(-1j*2*pi*cfo_estimate.*tk(1:nblocks*block_length)*symbol_rate);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% samps             signal samples x[k] [complex vector]
%
% mpower            power M to which x[k]x*[k - d] will be elevated to
%                       suppress the modulation in an M-PSK signal
%                       [integer scalar]
%
%                       Use mpower = 4 for QPSK.
%
% sample_delay      delay d between the samples x[k] and x[k - d] used to
%                       perform the estimation [integer scalar]
%
%                       In general d = 1 is used.
%
% block_length      block length L over which averaging is performed
%                       [integer scalar]
%
%                        By default and when varargin = 'fir', an L-tap FIR
%                        filter is used to perform the averaging 
%
%                        Since the delay of the estimator is (L + d - 1)/2,
%                        one should ensure it is an integer number in
%                        order to ease the compensation (otherwise some
%                        interpolation would be required).                     
%
%                        L + d should therefore be an odd number.
%
%                        When varargin = 'block', averaging is performed
%                        block-by-block over adjacent blocks of length
%                        block_length.
%                        In this case the estimator delay is irrelevant.
%
% varargin            averaging method [string / optional parameter]
%
%                        Use moving average over the entire vector of input
%                        samples with a FIR filter (physical 
%                        implementation) or block-by-block averaging.
%
%                        varargin = 'fir'
%                           moving average with an FIR filter of length
%                           block_length
%
%                        varargin = 'block' 
%                           block-by-block averaging over block_length 
%                           samples
%                           Note that if the number of input samples is not
%                           an integer multiple of block_length, the vector
%                           of input samples will be truncated and the
%                           estimator will be returned over the largest
%                           integer multiple of block_length.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% cfo_estimate      normalised estimate of the CFO [real vector]
%
%                       The length of cfo_estimate is the same as the
%                       length of the input vector samps when the averaging
%                       method is 'fir'. It is truncated to an integer
%                       number of blocks of size block_length when the
%                       averaging method is 'block'.
%
%                       The CFO is normalised by the sampling
%                       frequency. 
%
%                       If fs is the sampling frequency and
%                       Ts = 1/fs is the time interval between
%                       consecutive samples, then the actual CFO (in
%                       Hz) is
%                       CFO [Hz] = cfo_estimate*fs.
%
%                       Due to the delay of the estimator, the first 
%                       samples are not necessarily meaningful.
%
% estimator_delay   delay of the estimator, in number of samples
%                       [real scalar]
%
%                       Unfortunately, if L + d is even, the estimator
%                       delay is not an integer number of samples.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if (nargin == 4)

    method = 'fir';

elseif (nargin == 5) 

    method = varargin{1};

else

    error('cfo_diffphase: wrong number of input arguments');

end


z = (samps.*conj(dsp_delay(samps,sample_delay))).^mpower;
% Multiplication of samples and delayed-sample elevated to the m^th power


switch method
    % Averaging method

    case 'fir'
        % Runing averaging using an FIR filter

        cfo_estimate = angle(dsp_fir_linear(z,ones(1,block_length),zeros(1,block_length))/block_length)/2/pi/mpower/sample_delay;
        % Normalised CFO estimate f/fs
        % The actual CFO (in Hz) is this normalised estimate divided by the
        % sampling period (i.e. typically multiplied by the symbol rate for symbol
        % rate processing).

        estimator_delay = (block_length + sample_delay - 1)/2;
        % Estimator delay, in samples

    case 'block'
        % Block averaging

        nblocks = floor(length(samps)/block_length);
        % Number of complete blocks of length block_length in the input data

        cfo_estimate_block = zeros(1,nblocks);
        % Pre-initialise

        for iblock = 1:nblocks

            cfo_estimate_block(iblock) = angle(mean(z((iblock - 1)*block_length + 1:iblock*block_length)))/2/pi/mpower/sample_delay;

        end
        % CFO estimation by blocks

        cfo_estimate = cfo_estimate_block(ones(1,block_length),:);
        cfo_estimate = cfo_estimate(:).';
        % Create a vector of CFO estimates that has the same size as the
        % number of (truncated) input samples

        estimator_delay = 0;
        % Not relevant

    otherwise

        error('cfo_diffphase: averaging option not implemented');

end


end