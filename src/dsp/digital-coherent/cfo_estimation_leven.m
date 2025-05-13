function [cfo_estimate,estimator_delay] = cfo_estimation_leven(samps,mpower,sample_delay,block_length)
% CFO estimation for M-PSK modulation according to Leven
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements the carrier frequency offset (CFO) estimation 
% method based on the phase differential algorithm described in:
% A. Leven, N. Kaneda, U.-V. Koc, and Y.-K. Chen, "Frequency estimation in 
% intradyne reception," IEEE Photonics Technology Letters, vol. 19, no. 6,
% pp. 366-368, Mar. 2007.
% http://dx.doi.org/10.1109/LPT.2007.891893
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% tk = [0:1:length(symbs) - 1]/symbol_rate;
% cfo_estimation_mpower = 4;
% cfo_estimation_sample_delay = 1;
% cfo_estimation_block_length = 1000;
% [cfo_estimate, cfo_estimation_delay] = cfo_estimation_leven(symbs_rx,cfo_estimation_mpower,cfo_estimation_sample_delay,cfo_estimation_block_length);
% symbs_comp = dsp_delay(symbs_rx,cfo_estimation_delay).*exp(-1j*2*pi*cfo_estimate.*dsp_delay(tk,cfo_estimation_delay)*symbol_rate);
% symbs_comp = symbs_comp([2*cfo_estimation_delay + 1:end]);
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
%                        An L-tap FIR filter is used to perform the
%                        averaging.
%
%                        Since the delay of the estimator is (L + d - 1)/2,
%                        one should ensure it is an integer number in
%                        order to ease the compensation (otherwise some
%                        interpolation would be required).
%
%                        L + d should therefore be an odd number.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% cfo_estimate      normalised estimate of the CFO [real vector]
%
%                       The length of cfo_estimate is the same as the
%                       length of the input vector samps.
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
%                       Infortunately, if L + d is even, the estimator
%                       delay is not an integer number of samples.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

z = (samps.*conj(dsp_delay(samps,sample_delay))).^mpower;

cfo_estimate = angle(dsp_fir_linear(z,ones(1,block_length),zeros(1,block_length))/block_length)/2/pi/mpower/sample_delay;
% Normalised CFO estimate f/fs
% The actual CFO (in Hz) is this normalised estimate divided by the
% sampling period (i.e. typically multiplied by the symbol rate for symbol
% rate processing).  

estimator_delay = (block_length + sample_delay - 1)/2;
% Estimator delay, in samples
        
end