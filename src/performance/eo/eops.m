function [eop,eo_norm,opt_sampling] = eops(sig,bit_pattern,norm_power,eo_b2b)
% Calculation of eye opening penalty at optimum sampling time
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates eye opening and eye opening penalty for
% electrical binary signals (i.e. after photodetection). 
% The code was first written in C by Christian Rasmussen in his PhD work
% "Transmission analysis in WDM networks", Technical University of Denmark,
% 1999, then adapted to Matlab for use with VPI co-simulation by Beata 
% Zsigri.
% In this method, the eye opening is defined as the maximum height of the
% eye over a single sample. The function therefore returns the eye opening
% at the optimum sampling time.
% This function is kept for compatibility. The eopw function allows more 
% options for the eye opening calculations and calls the present function
% for method = 'cr'.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [eop,eo_norm,opt_sampling] = eops(abs(sig.x).^2,bit_pattern,params_eop.norm_power,-1)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               electrical signal to be characterised [real vector] 
%
% bit_pattern       binary sequence carried by the signal [binary vector]
%
% norm_power        optical power used for normalisation [real scalar]
%
%                       This normalisation value is set as an input 
%                       variable since this may allow the definition of 
%                       suitable eye openings for detection schemes using 
%                       a delay interferometer and balanced detection, 
%                       for instance.
%
% eo_b2b            back-to-back normalised eye opening [real scalar]
%
%                       Use -1 if the back-to-back eye opening is unknown 
%                       and should be determined by calling this function.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% eop               eye opening penalty, in dB [real scalar]
%
%                       Returns 0 in case the back-to-back eye opening is 
%                       calculated.
%
% eo_norm           eye opening normalised to the normalisation power 
%                   norm_power, in linear scale [real scalar]
%
% opt_sampling      index of the optimum sample where the eye opening is
%                       maximal [integer scalar]
%
%                       The index is in the range [1 nsamples_per_bit]
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

nbits = length(bit_pattern);
% Number of bits in the input binary pattern

nsamples_per_bit = nsamples/nbits;
% Number of samples per bit

nzeros = nbits - sum(bit_pattern);
% Number of 'spaces' in the pattern    

max_eye_opening = 0;
opt_sampling = -1; 
% Initialisation of the eye opening and optimum sampling time values

for sampling_time = 1:1:nsamples_per_bit
    
   sampled_sig = sig(sampling_time:nsamples_per_bit:end);
   % Collect samples of the signal at one particular sampling_time
   % (expressed in terms of samples) within the bit slot
   
   ra = sort(sampled_sig);
   % Sort the samples in ascending order.  
   
   level = 0.5*(ra(nzeros) + ra(nzeros + 1));
   % Threshold level, taken at the middle of the highest 'space' and lowest
   % 'mark'
   
   retrieved = (sampled_sig > level);
   % Retrieved bit pattern
   % Here it is somehow assumed that the samples corresponding to 
   % indices 1...nzeros effectively correspond to 'spaces'
   % while the samples corresponding to indices nzeros + 1...nbits are
   % 'marks'. This would not be the case if the eye is closed. However this
   % will be checked later since in this case one would not retrieve the
   % proper bit sequence.
   
   % We check synchronisation between the retrieved data and the original
   % bit pattern:
   sync = 0;
   for shift = 0:1:nbits - 1
       sync = (sum(xor(circshift(bit_pattern,[0 shift]),retrieved)) == 0);
       if sync
           break
       end       
   end        
     
   if sync
       % In case the initial bit pattern has indeed been retrieved, we
       % check if the eye opening at this particular sampling time is 
       % larger than the one previously stored. If yes, we update the eye 
       % opening and the optimum sampling time. Otherwise we do nothing and 
       % let the algorithm run again for the next sampling time (if any).
      if (ra(nzeros + 1) - ra(nzeros) > max_eye_opening)
         max_eye_opening = ra(nzeros + 1) - ra(nzeros);
         opt_sampling = sampling_time;
      end
   end
end
% End of loop over sampling time

eo_norm = max_eye_opening/norm_power;
% Normalise the eye opening with the normalisation power


if max_eye_opening ~= 0
   if eo_b2b ~= -1
      % The back-to-back eye opening is already known and the function
      % returns the eye opening penalty.
      eop = 10*log10(eo_b2b/eo_norm);
   else
      % The normalised back-to-back eye opening is returned by the 
      % function. The eye opening penalty is therefore zero dB.  
      eop = 0;
   end
else
   % No eye opening is found. The EOP is set to + infinity.
   eop = +Inf;
   eo_norm = +Inf;
   opt_sampling = NaN;
end