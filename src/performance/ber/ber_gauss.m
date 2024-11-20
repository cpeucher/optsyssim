function  [out_ber,out_threshold,out_sample] = ber_gauss(sig,pattern,params)
% BER estimation using Gaussian approximation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function estimates the BER of a binary modulated electrical signal
% by fitting standard distributions to the marks and spaces. BER as a
% function of threshold and /or sampling instant, as well as optimum BER
% can be obtained.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_bergauss.ignore_bits_start = 0;
% params_bergauss.distribution = 'gauss';
% params_bergauss.threshold_mode = 'optimum';%'optimum_search';%'fixed';
% params_bergauss.threshold = 0.5;
% params_bergauss.sample_mode = 'optimum';%'fixed';
% params_bergauss.sample_index = [ ];
% [BER,threshold,sample] = ber_gauss(sig,bit_pattern,params_bergauss); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
%                       The signal does not need to be retimed before
%                       calling this function.
%
% pattern           bit sequence against which the received signal should 
%                       be compared [binary vector]
%
% params            BER calculation parameters [structure]
%
%                       params.ignore_bits_start
%                           number of bits to ignore at the start of the 
%                           sequence for BER evaluation [integer scalar]
%
%                       params.distribution
%                           type of distribution fitted to the received 
%                           data for BER evaluation [string]
%
%                           params.distribution = 'gauss' 
%                               a single Gaussian distribution is used for 
%                               the marks, and another Gaussian 
%                               distribution is used for the spaces.
%
%                               Known to overestimate the BER in case ISI 
%                               is present.
%
%                       params.threshold_mode
%                           method used to determine the decision threshold
%                           [string]
%
%                           params.threshold_mode = 'optimum'
%                               the optimum threshold is calculated 
%                               analytically based on the estimated means 
%                               and standard deviations of the marks and 
%                               the spaces.
%
%                           params.threshold_mode = 'optimum_search'
%                               the optimum threshold is searched 
%                               numerically by sweeping over threshold 
%                               values and identify the ones that returns
%                               the lowest estimate of the BER.
%
%                           params.threshold_mode = 'fixed'
%                               fixed threshold values, specified in the 
%                               parameter params.threshold are used.
%
%                       params.threshold
%                           vector of threshold values in case fixed 
%                           thresholds are used []
%
%                       params.sample_mode
%                           method used to determine the sampling time
%                           [string]
%
%                           params.sample_mode = 'optimum'
%                               optimum sampling time, in which case the 
%                               BER is calculated as a function of sampling 
%                               time and the minimum value is returned.
%
%                           params.sample_mode = 'fixed'
%                               fixed sampling time, as specified in the 
%                               params.sample_index parameter.
%
%                       params.sample_index
%                           values of the sample index at which
%                           the BER will be evaluated. 
%                           Single value or vector of values in
%                           [1 nsamples_per_symbol].
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% out_ber           estimated bit-error-ratio [real scalar]
%
% out_threshold     corresponding decision threshold [real scalar]
%
% out_sample        corresponding sample index [integer scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = length(sig);
% Number of samples in the input signal
nbits = length(pattern);
% Number of bits in the input sequence
nsamples_per_bit = nsamples/nbits;
% Number of samples per bit

ref_sig = elec_coder_nrz(pattern,nsamples_per_bit);
% Create electrical reference signal corresponding to the expected
% sequence

[~,indexmax] = max(xcorr(sig,ref_sig));
% Calculate the delay shift that maximises cross-correlation between the
% signal and the reference signal
sig = circshift(sig,[0 -indexmax]);
% Resynchronise the signal. 

pattern = pattern(1 + params.ignore_bits_start:nbits);
% Remove the bits that are ignored at the begining of the sequence.
sig = sig(1 + params.ignore_bits_start*nsamples_per_bit:length(sig));
% Remove the part of the input signal corresponding to the ignored bits
nbits = nbits-params.ignore_bits_start;
% Update the number of bits


if strcmp(params.distribution,'gauss')
    % Single Gaussian distributions for the marks and for the spaces
    
    if (strcmp(params.sample_mode,'optimum') && isempty(params.sample_index))
        sample_index = 1:1:nsamples_per_bit;
    else
        sample_index = params.sample_index;
    end
    % In case the optimum sampling time is looked for and a range of sample
    % indices is not specified, one should loop over all samples in a bit
    % duration
    % Else the specified sample_index array is considered.
    
    if (strcmp(params.threshold_mode,'optimum_search') && isempty(params.threshold))
        %[maxspaces,minmarks]=BER_ThresholdRange(sig,prbs,sampleindex,BitsToIgnoreStart,Nbits)
        threshold = linspace(min(sig),max(sig),100);
    else
        threshold = params.threshold;
    end
    % In case the optimum threshold is searched (i.e. not calculated) and a
    % range of thresholds is not calculated
    % Else the specified Threshold array is considered.
    
    
    for isample = 1:length(sample_index)
        % Loop over sampling time.
        
        marks = [];
        spaces = [];
        % Initialise arrays to collect the mark and spaces samples
        
        sig_samples = sig(sample_index(isample)+((1:nbits)-1)*nsamples_per_bit);
        % Collect samples of the signal, corresponding to a given
        % SampleIndex.
        
        ref_sig_samples = ref_sig(sample_index(isample)+((1:nbits)-1)*nsamples_per_bit);
        % Collect samples of the reference signal, corresponding to the
        % same SampleIndex.
        
        marks = sig_samples(ref_sig_samples == 1);
        spaces = sig_samples(ref_sig_samples == 0);
        % Distribute the signal samples in 2 categories: marks and
        % spaces.
        
        I1 = mean(marks);
        % Mean of the marks

        I0 = mean(spaces);
        % Mean of the spaces

        sigma1 = std(marks,1);
        % Standard deviation of the marks

        sigma0 = std(spaces,1);
        % Standard deviation of the spaces
        % For a normal distribution, the sample mean (given by the
        % Matlab mean(X) function), and the sample variance (given by
        % the square of the Matlab std(X,flag=1) function) are maximum
        % likelihood estimates.
        
        % We now have the means and standard deviations of the marks and
        % space at the specified sampling time. We can either look for the
        % optimum threshold and corresponding BER, or sweep over specified
        % threshold values
        
        
        if strcmp(params.threshold_mode,'optimum')
            % In case the optimum threshold is used, no need to loop over
            % thresholds
            
            id(isample) = (sigma0*I1 + sigma1*I0) / (sigma0 + sigma1);
            % Determination of optimum threshold
            Q = (I1 - I0)/(sigma1 + sigma0);
            % Q-factor at the optimum threshold
            ber(isample) = 0.5*erfc(Q/sqrt(2));
            % BER at the optimum threshold
            
            
        elseif (strcmp(params.threshold_mode,'fixed') || strcmp(params.threshold_mode,'optimum_search'))
            % In case of fixed specified threshold values, we loop over
            % those values.          
            
            for ithreshold = 1:length(threshold)
                % Loop over all thresholds.
                ber(ithreshold,isample) = 0.25*(erfc((I1-threshold(ithreshold))/sqrt(2)/sigma1)+erfc((threshold(ithreshold)-I0)/sqrt(2)/sigma0));
                % Estimate BER
            end
            % End of loop over thresholds
            
        else
            error('BER_Gauss: threshold mode not implemented.');
        end
        % End of selection of threshold mode
        
    end
    % End of loop over sampling times
    
    
    %----------------------------------------------------------------------
    % Post-processing of the BER vs threshold and sampling time data,
    % depending on the options.
    %----------------------------------------------------------------------    
    if strcmp(params.sample_mode,'optimum')
        if strcmp(params.threshold_mode,'optimum')
            % If params.threshold_mode='optimum', then BER is an array 
            % containing the optimum BER at each of the sampling times. 
            % The optimum threshold for each sampling time is saved in the 
            % Id array.
            [out_ber,out_sample] = min(ber);
            out_threshold = id(out_sample);
        elseif strcmp(params.threshold_mode,'optimum_search')
            % If params.threshold_mode='optimum_search', then BER is a 2 dimensional
            % array containing the BER as a function of threshold and
            % sampling time: BER(nthresh,nsample). The function should then
            % look for the optimum BER and return the corresponding
            % sampling time and threshold
            [ber,isample] = min(ber,[],2);
            % Returns the BER at the optimum sampling time for each
            % threshold value. isample corresponds to the index of the
            % sample where the minimum BER is reached (for each threshold).
            [out_ber,ithresh] = min(ber);
            % Returns the minimum BER. ithresh corresponds to the index of
            % the optimum threshold value
            out_threshold = threshold(ithresh);
            out_sample = isample(ithresh);
        elseif strcmp(params.threshold_mode,'fixed')
            % If params.threshold_mode='fixed', then BER is a 2 dimensional array
            % containing the BER as a function of threshold and sampling
            % time: BER(nthresh,nsample). The function should then return
            % the BER at the specified threshold values and at the optimum
            % sampling time. The OutSample array will contain the optimum
            % sample index for each threshold.
            [out_ber,out_sample] = min(ber,[],2);
            out_threshold = threshold;
        end
    elseif strcmp(params.sample_mode,'fixed')
        if strcmp(params.threshold_mode,'optimum')
            % If params.threshold_mode='optimum', then BER is an array containing
            % the optimum BER at each of the sampling times. The optimum
            % threshold for each sampling time is saved in the Id
            % array. The function should return the optimum threshold for
            % each sampling time and the corresponding BER.
            out_ber = ber;
            out_threshold = id;
            out_sample = sample_index;
        elseif strcmp(params.threshold_mode,'optimum_search')
            % If params.threshold_mode='optimum_search', then BER is a 2 dimensional
            % array containing the BER as a function of thereshold and
            % sampling time: BER(nthresh,nsample). The funciton should then
            % look for the optimum BER and return the corresponding
            % sampling time and threshold.
            [out_ber,ithresh] = min(ber,[],1);
            % Retuns the optimum BER for each sampling time. ithresh is the
            % index of the sample where the optimum is reached.
            out_threshold = threshold(ithresh);
            out_sample = sample_index;
        elseif strcmp(params.threshold_mode,'fixed')
            % If params.threshold_mode='fixed', then BER is a 2 dimensional array
            % containing the BER as a function of threshold and sampling
            % time: BER(nthresh,nsample). The function should then return
            % the BER at the specified threshold values and at the
            % specified sampling times.
            % The OutSample array will simply reproduce the input Sample
            % array. The OutThreshold array will reproduce the input
            % Threshold array.
            out_ber = ber;
            out_sample = sample_index;
            out_threshold = threshold;
        end
    end
    % Decide what to output depending on the options.
    
    
else
    
    error('ber_gauss: distribution not implemented (yet).');
    % So far only Gaussian distribution is implemented.
    
end
% End of selection of distribution.

end