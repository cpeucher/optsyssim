function [out,freq] = psd_welch(data,block_size,overlap_samples,window_type,fft_length,nsided,out_type,varargin)
% Power spectral density estimation using the Welch periodogram method
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function estimates the power spectral density (PSD) of a random
% sequence by calculating its periodogram according to Welch's method.
% See for instance:
% P. Stoica and R.L. Moses, Spectral Analysis of Signals (Pearson/Prentice 
% Hall, Upper Saddle River, N.J, 2005), chapter 1 and 2 for more details on
% the method and the definitions used in this function.
% The block size, number of overlapping samples, FFT length and type of
% windowing function can all be set. 
% Both single-sided and double-sided PSDs can be returned for normalised
% and physical frequency axes, depending on whether the sampling frequency
% fs is provided as an input parameter.
% 
% Let phi(^w) be the "power spectral density" as defined in Stoica,
% eqs. (1.3.7) and (1.3.10) as a function of normalised angular frequency
% -pi < ^w < pi. 
% The average power of the input sequence {y[n]} is 
% r(0) = E[|y[n]|^2] = (1/2/pi) int_{-\pi}^{pi} phi(^w) d^w      (A)
% The periodogram as calculated in this function is an estimate of phi(^w).
%
% However, the vector out returned by the function is NOT the periodogram.
%
% 1) In case no sampling frequency is provided as an optional input, the
% function returns phi(^w) /2/pi. This is the quantity to integrate in
% order to obtain the power in a given bandwidth according to equation (A)
% above.
% phi(^w) d^w /2/pi is the infinitesimal power in a bandwidth d^w around
% ^w.
% Observe that even though the definition of the PSD provides phi(^w),
% one needs to integrate phi(^w) /2/pi to obtain the power in a given
% normalised angular frequency bandwidth.
% The value returned by the function is actually phi(^w) /2/pi, so that it
% indeed corresponds to a spectral density in units of power /
% radian-per-sample.
% For a zero-mean process, var{y} = \int_{^w1}^{^w2} out d^w.
% If {y[n]} is a white noise sequence, the density out returned by the
% function can be predicted to be out_predicted that can be calculated
% according to:
% var{y} = out_predicted*pi if single-sided, or
% var{y} = our_predicted*2*pi if double-sided.
%
% 2) In case fs is provided as an input, the function returns the power
% spectral density in units of power/Hz. The power in a given bandwidth can
% directly be obtained by integration of this density.
% If {y[n]} is a white noise sequence resulting from the sampling of a real
% process y(t) with sampling frequency fs, then the density returned by the
% function can be predicted according to:
% var{y} = out_predicted*fs/2 if single-sided, or
% var{y} = out_predicted*fs if double-sided.
%
% Some further conversions can be made, depending on the considered
% frequency axis.
% In what follows [out, freq] are the vectors returned by a call to the
% function psd_welch.
% A) if no sampling frequency is provided as optional input parameter, and
% if we wish the PSD as a function of
% - normalised angular frequency axis (-pi < ^w < pi), in
%   radians-per-sample
%   out is the PSD in power/radian-per-sample, i.e. phi(^w)/2/pi.
%   PSD = out               => PSD in power/radian-per-sample
%   freq                    => normalised angular frequency in
%                              radian-per-sample
% - normalised frequency axis (-1/2 < ^f < 1/2), in cycles-per-sample
%   The PSD in power/cycle-per-sample can be obtained by multipling out
%   obtained in a normalised call (no fs provided) by 2*pi and scaling the
%   frequency axis.
%   PSD = out*2*pi          => PSD in power/cycle-per-sample
%   freq = freq/2/pi        => normalised frequency in cycle-per-sample
% B) if fs is provided as an optional input parameter, and if we wish the
% PSD as a function of
% - physical frequency axis (-fs/2 < f < fs/2), in Hz
%   PSD = out               => PSD in power/Hz
%   freq                    => frequency in Hz
% - physical angular frequency axis (-pi*fs < w < pi*fs), in
%   radians-per-second
%   PSD = out/2/pi          => PSD in power/radian-per-second
%   freq = freq*2*pi        => angular frequency in radian-per-second
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% block_size = 100;
% overlap_samples = 50;
% fft_length = 640;
% window_type = 'hamming';
% [psd,norm_freq] = psd_welch(data,block_size,overlap_samples,window_type,fft_length,'one_sided','psd');
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% data              input sequence [real vector]
%
% block_size        size of the block over which periodograms will be
%                       calculated [integer scalar]
%
% overlap_samples   number of overlapping samples between consecutive
%                       blocks [integer scalar]
%
% window_type       type of window, e.g. 'hann','hamming','bartlett'...
%                       [string]
%
% fft_length        length of FFT used for the calculation of each
%                       periodogram
%
% nsided            decides whether one-sided or two-sided psd is 
%                       calculated [string]
%
%                       nsided = 'one_sided';    (default value)
%                       nsided = 'two_sided';
%
% out_type          determines wich quantity is returned by the 
%                       function [string]
%
%                       out_type = 'psd';    (default value)  
%                           The PSD is returned.
%                       out_type = 'power';  
%                           The power spectrum is returned.
%
% varargin          sampling frequency, in Hz [real scalar]
%
%                       To be input in case the function 
%                       should return the PSD / power spectrum as a 
%                       function of the actual frequency, and not the
%                       normalised angular frequency.
%                       Providing a sample frequency value as optional
%                       parameter will automatically result in PSD /
%                       power spectrum as a function of actual
%                       physical frequency being returned.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% out               power spectral density [real vector] in
%                       - power/radian-per-sample if no sampling frequency
%                         is provided as an optional input parameter
%                       - power/Hz if a sampling frequency is provided as
%                         an optional input parameter 
%
% freq              frequency-related quantity [real vector]
%                       - normalised angular frequency, in
%                         radian-per-sample if no sampling frequency
%                         is provided as an optional input parameter
%                       - frequency in Hz if a sampling frequency is 
%                         provided as an optional input parameter 
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% CHECK INPUTS
% -------------------------------------------------------------------------
if nargin < 1
    error('psd_welch: not enough input arguments.');
end

if nargin < 2
    block_size = 256;
end

if nargin < 3
    overlap_samples = floor(block_size/2);
end

if nargin < 4
    window_type = 'hann';
end

if nargin < 5
    fft_length = 2^nextpow2(block_size);
end

if nargin < 6
    nsided = 'one_sided';
end

if nargin < 7
    out_type = 'psd';
end

if nargin < 8
    normalised = 1;
    % If the sampling frequency is not provided, the normalised PSD is
    % returned.
elseif nargin == 8
    % If the sampling frequency is provided, we assume the user is
    % interested in the non-normalised PSD.
    normalised = 0;
    fs = varargin{1}; 
    % Sampling frequency, in Hz.
else
    error('psd_welch: too many input arguements.');
end



if isempty(block_size)    
    block_size = 256;    
end

if isempty(overlap_samples)
    overlap_samples = floor(block_size/2);
end

if isempty(window_type)
    window_type = 'hann';
end

if isempty(fft_length)
    fft_length = 2^nextpow2(block_size);
end

if isempty(nsided)
    nsided = 'onesided';
end

if isempty(out_type)
    out_type = 'psd';
end

if normalised == 0 && isempty(fs)
    error('psd_welch: sample rate is required for non-normalised output.');
end

data = data(:).';
% Ensure data is a line vector

% -------------------------------------------------------------------------
% FUNCTION CORE
% -------------------------------------------------------------------------
nsamples = length(data);
% Number of samples in the input data

x = data(((1:block_size).'+(0:block_size - overlap_samples:nsamples - block_size)).');
% Split the input data in (overlaping) blocks
% The blocks are saved as the lines of a matrix of dimension nblocks x
% block_size
dimx = size(x);
nblocks = dimx(1);
% Number of blocks

window.type = window_type;
window.length = block_size;
window.symmetry = 'periodic';
w = dsp_window(window);
% Create window


s1 = sum(abs(w));
% Normalisation factor. Linked to the so-called "coherent gain" of the
% window: CG = S1 / N where N is the window length.
s2 = sum(abs(w).^2);
% Normalisation factor. Linked to the so-called "noise gain" of the window:
% NG = S2 / N where N is the window length.

yj = abs(fft(x.*w,fft_length,2)).^2 / s2;
% Calculate the periodogram of each windowed bock (windowed periodogram).
% The periodogram is in fft order.

pp = mean(yj,1);
% Average over all periodograms


norm_freq = (0:1:fft_length-1)/fft_length;
% Normalised frequency axis, in increasing frequency order, in cycles per
% sample
% The actual frequency is fs * norm_freq.

% if mod(fft_length,2)
%     % odd value of fft_length
%     % mod(fft_length,2) = 1 
%     
%     normFreq((fft_length+1)/2+1:end) = normFreq((fft_length+1)/2+1:end) - 1;  
%     
% else
%     % even value of fft_length
%     % mod(fft_length,2) = 0
%     
%     normFreq(fft_length/2+1:end) = normFreq(fft_length/2+1:end) - 1;   
%     
% end
% Code above can be replaced by single line below.

norm_freq((fft_length+mod(fft_length,2))/2+1:end) = norm_freq((fft_length+mod(fft_length,2))/2+1:end) - 1;   
% normFreq now contains the normalised frequencies in fft order, i.e.
% matching the order of vector pp.
% norm_freq(1) = 0
% norm_freq(2) = 1/fft_length
% norm_freq(fft_length/2) = 0.5 - 1/fft_length                 fft_length even
% norm_freq(fft_length/2+1) = -0.5                             fft_length even
% norm_freq(fft_length) = -1/fft_length                        fft_length even
% norm_freq((fft_length+1)/2) =  0.5 - 1/fft_length/2          fft_length odd
% norm_freq((fft_length+1)/2+1) = -0.5 + 1/fft_length/2        fft_length odd
% norm_freq(fft_length) = -1/fft_length                        fft_length odd


if strcmp(nsided,'one_sided')
    % Single-sided PSD
    out = pp(1:(fft_length+mod(fft_length,2))/2);
    % Part of the 'spectrum' corresponding to positive frequencies.
    % In case fft_length is odd, the Nyquist frequency does not appear in
    % the frequency vector.
    % In case fft_length is even, the Nyquist frequency (-fs/2 or -1/2 if
    % normalised frequencies are considered) appear in the frequency array
    % at index fft_length/2 +1.
    % In both cases, positive frequencies correspond to the indices:
    % 1:(fft_length+mod(fft_length,2))/2
    % However, we want to add the Nyquist frequency when fft_length is even,
    % to replicate the behaviour of the Matlab Signal Processing Toolbox 
    % pwelch function.
    % Therefore we need to consider the range of indices:
    % 1:(fft_length+mod(fft_length,2))/2+1-mod(fft_length,2)
    % which is:
    % 1:(fft_length-mod(fft_length,2))/2+1
    out = pp(1:(fft_length-mod(fft_length,2))/2+1);
    % This includes the value of the periodogram at the Nyquist frequency
    % -fs/2 when fft_length is even, but is strictly limited to positive
    % frequencies when fft_length is even (in which case the Nyquist
    % frequency does not appear in the frequency vector).
    norm_freq = norm_freq(1:(fft_length-mod(fft_length,2))/2+1);
    % These are the positive frequencies if fft_length is odd, or
    % the positive frequencies + -1/2 if fft_length is even.
    if mod(fft_length,2)
        % Odd value of fft_length
        out(2:end) = 2*out(2:end);
        % Factor 2 for single-sided PSD for all points apart the one
        % corresponding to dc. The Nyquist frequency is not present in the
        % frequency vector for odd number of points.        
    else
        % Even value of fft_length
        norm_freq(end) = norm_freq(end) + 1;
        % Replace the Nyquist frequency value -1/2 by 1/2
        out(2:end-1) = 2*out(2:end-1);
        % Factor 2 for single-sided PSD for all points apart those
        % corresponding to dc and the Nyquist frequency, which appear only
        % once in the spectrum.        
    end

elseif strcmp(nsided,'two_sided')
    % Double-sided PSD
    % The norm_freq and out vectors are in increasing frequency order.
    out = fftshift(pp);
    norm_freq = fftshift(norm_freq);
    
else
    error('psd_welch: the spectrum should either be "one_sided" or "two_sided". No other choice is possible!');    
end
% End of switch over single/double-sided PSD / spectrum


switch out_type        
    
    case 'psd'
        % The output is the power spectal density.                 
        
        if normalised
            out = out/2/pi;
            % psd in power/radian-per-sample
            freq = 2*pi*norm_freq;
            % Normalised angular frequency radian-per-sample
        else
            out = out/fs;
            % psd in power/Hz
            freq = norm_freq*fs;
            % Frequency in Hz
        end             
  
        
    case 'power'
        % The output is the power.
        
        if normalised
            out = out * s2 /s1^2;            
            freq = 2*pi*norm_freq;            
            
        else            
            out = out * s2 /s1^2;            
            freq = norm_freq*fs;
        
        end

end
% End of switch over output type

end