function eyeop = eopw(sig,params)
% Calculation of eye opening penalty
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the eye opening and eye opening penalty for
% electrical binary signals (i.e. after photodetection). It can search
% either the maximum eye opening at a single sampling point, or return the
% height of the rectangle of maximum area with a specified time width that
% fits within the eye.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_eop.norm_power = 2*char_opt_average_power(sig);
% sig = abs(sig.x).^2;
% params_eop.method = 'eow';%'cr';
% params_eop.bit_pattern = bit_pattern;
% params_eop.eo_b2b = -1;%eo_norm;
% params_eop.display_eye = 1;
% params_eop.eye_width_samples = 12;
% params_eop.bits_to_ignore_start = 0;
% params_eop.eye_display_name = 'Eye diagram';
% eop = eopw(sig,params_eop); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input electrical signal [real vector]
%
% params            eye opening calculation parameters [structure]
% 
%                       params.method
%                           method used for the eye opening 
%                           calculation [string]
%
%                           params.method = 'cr'
%
%                               The method first written in C by Christian 
%                               Rasmussen, in his PhD thesis "Transmission 
%                               analysis in WDM networks", DTU, 1999, then
%                               adapted to Matlab for use with VPI 
%                               co-simulation by Beata Zsigri is adopted.
%
%                               For this method, the eye opening is defined
%                               as the maximum height of the eye over a
%                               single sample.
% 
%                               This legacy method is faster if one is
%                               interested in the eye opening defined at a
%                               single sample.
%                       
%                               This is just a call to eops function.
%                          
%
%                           params.method ='eow'
%
%                               Eye opening for a specified opening width.
%                               The eye opening is defined at the maximum
%                               height that fits the eye with a given width, 
%                               specified in terms of number of samples. 
%
%                               If the width is equal to 1 sample, then the 
%                               definition coincides with that of the 'cr' 
%                               method, even though the eye opening is 
%                               calculated in a different way.
%
%                       params.bit_pattern
%                           binary sequence carried by the signal 
%                               [binary vector]
%                           
%                       params.norm_power
%                           optical power used for normalisation, in W
%                               [real scalar]
%
%                           This normalisation value is set as an input 
%                           variable since this may allow the definition 
%                           of suitable eye openings for detection schemes 
%                           using a delay interferometer and balanced 
%                           detection, for instance.
%
%                           For an OOK signal, a standard normalisation is 
%                           2*Pav where Pav is the optical signal average
%                           power. Then a signal with infinite extinction 
%                           ratio and no distortion will have a normalised
%                           eye opening of 1 (or 0 dB), since 
%
%                           EOnorm = (P1 - P0)/2Pav
%                                  = (P1 - P0)/(P1 + P0)
%                                  = (ER - 1)/(ER + 1)
%
%                       params.eo_b2b
%                           back-to-back normalised eye opening, in linear
%                           units [real scalar]
%
%                           params.eo_b2b = -1
%                               should be used if the back-to-back eye 
%                               opening is unknown and should be determined
%                               by calling this function.
%
%                           params.eo_b2b= eo_norm; 
%                               should be used if the eye opening eo_norm 
%                               is known from a previous call to the 
%                               function, for instance in order to 
%                               determine the eye opening penalty. 
%
%                               eo_norm corresponds to the value 
%                               eyeop.eye_norm 
%                               returned by the function.
%
%
%                       params.display_eye [0/1]
%
%                           Optionally display the eye diagram of the 
%                           analysed signal together with the eye opening 
%                           limits.
%
%                           This option is only implementd for 
%                           params.method = 'eow'.
%                       
%                       params.eye_width_samples
%                           width of the searched eye opening, expressed in 
%                           terms of number of samples [integer scalar]
% 
%                           Relevant when 
%                               params.method='eow'
%
%                           Using params.method = 'cr' or 
%                           params.method ='eow'
%                           with params.eye_width_samples = 1 should result
%                           result in the same output.
%
%                       params.bits_to_ignore_start
%                           number of bits to ignore at the start of the 
%                           signal [integer scalar] 
%
%                           Only used for the 'eow' method.
%
%                           In case the 'cr' method is used, then the bits
%                           to ignore need to be truncated in both the 
%                           input signal and the expected bit pattern 
%                           before calling the function. 
%
%                       params.eye_display_name
%                           name of the eye diagram display [string]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% eyeop             eye opening and related quantities [structure]
%
%                       eyeop.eop
%                           eye opening penalty, in dB [real scalar]
%
%                           Returns 0 in case the back-to-back eye opening
%                           is calculated.
%
%                       eyeop.eo_norm
%                           eye opening normalised to the normalisation 
%                           power params.norm_power, in linear scale
%                           [real scalar]
%
%                       eyeop.opt_sampling
%                           index of the optimum sample where the eye 
%                           opening is maximal [integer scalar]
%
%                           The index is in the range [1 nsamples_per_bit]
%
%                           In case params.method = 'eow', this is the 
%                           index of the first sample in the eye opening 
%                           rectangle.
%
%                           The width of rectangle is 
%                           (eye_width_samples - 1)*dt.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                    time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt

nsamples = length(sig);
% Number of samples in the input signal

nbits = length(params.bit_pattern);
% Number of bits in the input binary pattern

nsamples_per_bit = nsamples/nbits;
% Number of samples per bit



if strcmp(params.method,'eow')
    % Check the existence of required input parameters when the calculation
    % method is set to 'eow' 
    
    if ~isfield(params,'eye_width_samples')
        % Check whether the width of the searched eye opening is defined       
        error('eopw: the width of the searched eye opening should be specified in eow mode.');
        % We could arbitrarily set it to 0 so that the results are
        % equivalent to the ones returned using the 'cr' method, but we
        % prefer to make it straight and return an error.        
    else
        % The eye width is indeed specified as an input parameter.        
        if params.eye_width_samples > nsamples_per_bit
            % Return an error if the eye width is larger than the number of
            % samples per bit.
            error('eopw: eye opening width should be smaller than eye width.')
        end        
    end

    
    if ~isfield(params,'bits_to_ignore_start')
        % Check whether the number of bits to ignore at the start of the
        % sequence is defined.
        
        params.bits_to_ignore_start = 0;
        % If not, arbitrarily set the number of bits to ignore at the start
        % of the sequence to 0.        
    end
        
end
   
    
if strcmp(params.method,'cr')
    % The classic method by Christian Rasmussen. 
    % The eye opening corresponds to the maximum height of the eye taken at
    % one single sample value.
    
    [eyeop.eop,eyeop.eo_norm,eyeop.opt_sampling] = eops(sig,params.bit_pattern,params.norm_power,params.eo_b2b);
    % We simply call the orignal function adapted from Christian Rasmussen.
    
    
elseif strcmp(params.method,'eow')
    % The eye opening corresponds to the maximum height of a rectangle of
    % width params.eye_width_samples fitting into the eye diagramme.
    
    [sig,~] = rx_resynchronise(sig,params.bit_pattern,'elec');
    % Coarse resynchronisation
    
    bit_pattern_truncated = params.bit_pattern(params.bits_to_ignore_start + 1:nbits);
    sig = sig(params.bits_to_ignore_start*nsamples_per_bit + 1:nsamples_per_bit*nbits);
    nbits = nbits - params.bits_to_ignore_start;
    % Truncate pattern and signal to ignore the specified number of bits at
    % the start
    
    [sig,~] = rx_resynchronise(sig,bit_pattern_truncated,'elec');
    % Fine synchronisation.
    
    % Note: the synchronisation above is a little bit dirty since one
    % performs synchronisation before truncation. The untruncated waveform
    % may affect the outcome of the cross correlation used for
    % synchronisation. However experience shows it works fine provided the
    % number of bits is sufficienly large. We will very likely always work
    % with more than 128 bits, where it should be ok.
    
    for istart = 1:nsamples_per_bit - params.eye_width_samples + 1
        % Start looping over the position of the eye width
        
        for ibit = 1:nbits
            % Start looping over the bits

            samples(:,ibit) = sig(istart + (ibit - 1)*nsamples_per_bit:istart + params.eye_width_samples -1 + (ibit - 1)*nsamples_per_bit);
            % Collect samples belonging to the specified eye opening width
            % for each bit
            % samples is a matrix of nbits columns, each column containing
            % the samples collected for a given bit.
        end
        % End of loop over the bits
        
        
        minvalues = min(samples,[],1);
        % Minimum values of the samples collected for each bit
        % minvalues is a line vector of nbits columns. 
        % Each element contains the minimum value of the samples gathered 
        % for that particular bit.

        maxvalues = max(samples,[],1);
        % Maximum values of the samples collected for each bit. 
        % maxvalues is a line vector of nbits columns. 
        % Each element contains the maximum value of the samples gathered 
        % for that particular bit.
        
        % The eye opening is defined as the minimum of the ones minus the
        % maximum of the zeros (over the eye opening width, i.e. not
        % necessarily occuring at the same sample, unless the eye opening
        % width is equal to 1 sample).
        
        
        min_array = bit_pattern_truncated.*minvalues;
        % min_array contains the minimum value taken by a given bit if this
        % bit is a 1, or 0 if this bit is a 0.
        min_ones = min_array(find(min_array));
        % Removes the 0's from the min_array.
        % min_ones contains all the minimum values taken by all the 1 bits.
        
        max_array = double(not(bit_pattern_truncated)).*maxvalues;
        % max_array contains the maximum value taken by a given bit if this 
        % bit is a 0, or 0 if this bit is a 1.
        max_zeros = max_array(find(max_array));
        % Removes the 0's from the max_array. 
        % max_zeros contains all the maximum values taken by all the 0 bits.
        
        eo(istart) = min(min_ones) - max(max_zeros);
        baseline(istart) = max(max_zeros);
        topline(istart) = min(min_ones);
        % This is the eye opening for this particular position in the eye.       
        
    end
    % End of loop over the position of the eye width
    
    [max_eye_opening,eyeop.opt_sampling] = max(eo);
    % Maximum eye opening over all positions in the eye
    
    eyeop.eo_norm = max_eye_opening/params.norm_power;
    % Normalised to the signal power
    
    
    if max_eye_opening >= 0
        if params.eo_b2b ~= -1
            % The back-to-back eye opening is already known and the function
            % returns the eye opening penalty.
            eyeop.eop = 10*log10(params.eo_b2b/eyeop.eo_norm);
        else
            % The normalised back-to-back eye opening is returned by the
            % function. The eye opening penalty is therefore zero dB.
            eyeop.eop = 0;
        end
    else
        % No eye opening is found. The EOP is set to + infinity
        eyeop.eop = +Inf;
        eyeop.eo_norm = NaN;
        eyeop.opt_sampling = NaN;
    end
    
    
    %----------------------------------------------------------------------
    % display eye diagram
    %----------------------------------------------------------------------    
    if params.display_eye == 1         
        
        params_eye.neyes = 1;
        params_eye.nsamples_per_symbol = nsamples_per_bit;
        params_eye.save.txt = 0;
        params_eye.save.emf = 0;
        params_eye.save.jpg = 0;
        params_eye.save.vertical_scale_type = 'auto';%'fixed';
        params_eye.save.vertical_scale_max = 1.0e-3;
        params_eye.save.display_baseline = 1;
        params_eye.colour_grade = 0;
        params_eye.name = params.eye_display_name;
        meas_eye(sig,params_eye);
        
        hold on;
        
        if params.eye_width_samples == 1
            % If the eye opening is defined at only a single sample we
            % display a line.
            line(ones(1,10)*(eyeop.opt_sampling - 1)*dt/1.0e-12,linspace(baseline(eyeop.opt_sampling),topline(eyeop.opt_sampling),10),'LineWidth',2,'LineStyle','-','Color','r')        
        else
            % Else we display the limits of the eye opening vector.
            rectangle('Position',[(eyeop.opt_sampling-1)*dt/1.0e-12,baseline(eyeop.opt_sampling),(params.eye_width_samples-1)*dt/1.0e-12,max_eye_opening],'LineWidth',2,'LineStyle','-','EdgeColor','r');
        end       
        
    end
    % End of optional display eye diagram
    
    
else
    error('eopw: eye opening calculation method not implemented');
end

end