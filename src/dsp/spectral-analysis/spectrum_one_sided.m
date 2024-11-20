function [f,X1] = spectrum_one_sided(x,fs,varargin)
% Quick calculation of one-sided spectrum of a signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns the one-sided spectrum of the input signal as well
% as the corresponding positive frequency axis.
% It is designed to work regardless of whether the number of samples in the
% signal is odd or even.
% The power spectrum is abs(X1).^2.
% If one wishes to use an fft size different from the length of the input
% signal, the fft length can be specified as an optional input parameter.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [f,X1] = spectrum_one_sided(x,fs);
% [f,X1] = spectrum_one_sided(x,fs,2^10);
% [f,X1] = spectrum_one_sided(x,fs,2^10,1);
% [f,X1] = spectrum_one_sided(x,fs,[],1); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                samples of the signal [real vector]
%
% fs               sampling frequency of the input signal, in Hz
%                       [real scalar]
%
%                       To work with normalised quantities, input fs = 1.
%
% varargin         optional input arguments
%
%                       varargin{1}
%                           size of the fft [integer scalar]
%
%                           If no optional argument is specified, the size
%                           of the fft is taken equal to the length of the 
%                           signal.
%
%                       varargin{2}
%                           include the Nyquist frequency in the spectrum
%                           [0/1]
%
%                           Applies only when the fft length is even.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% f                 positive frequency vector [real vector]
%                       - in Hz if fs is provided in Hz
%                       - in normalised frequencies units (cycles per
%                         sample) if fs = 1 is provided as an input
%
% X1                one-sided spectrum
%                       One should calculate abs(X).^2 to get the power
%                       spectrum.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = length(x);
% Number of samples in the input signal

% 1st optional argument: FFT length (nfft)
% 2nd optional argument: include Nyquist (include_nyquist)

if nargin == 4
    
    if isempty(cell2mat(varargin(1)))        
        nfft = nsamples;        
    else
        nfft = varargin{1};
    end
    
    include_nyquist = varargin{2} ;
    
elseif   nargin == 3
    
    nfft = varargin{1};
    include_nyquist = 0;
    
elseif nargin == 2
    
    nfft = nsamples;
    include_nyquist = 0;
    
else
    error('spectrum_one_sided: wrong number of input arguments.');
end

if nfft < nsamples
    nnorm = nfft;
    % In case the FFT length is smaller than the samples vector length, the 
    % samples are truncated to the FFT length. We should normalize the
    % outcome of the fft function by the FFT length, not the original
    % samples vector length.
else
    nnorm = nsamples;
    % Otherwise, when the FFT length is larger than the number of samples,
    % the proper normalisation is the samples vector length.
end


X = fft(x,nfft)/nnorm;
% Double-sided spectrum, in fft order

X1 = X(1:ceil(nfft/2));
% Only retain the components corresponding to the positive frequencies.
% Valid for odd and even numbers of samples / fft lengths.
% In our convention, we do not include the Nyquist frequency fs/2 as a 
% positive frequency for an even number of samples (for an odd number of 
% samples, the Nyquist frequency does not belong to the sampled frequencies
% anyway).

X1(2:end) = sqrt(2)*X1(2:end);
% Ensure that the power is doubled for all frequency components, apart from
% DC, which only appears onece in the double-sided spectrum. As stated
% above, we do not need to take care of the Nyquist frequency since it is
% actually kept as -fs/2 and does not appear among the frequencies
% considered for the single-sided spectrum.

f = (0:ceil(nfft/2) - 1)*fs/nfft;
% Corresponding frequency array


if ~rem(nfft,2) && include_nyquist
    % We add Nyquist frequency if asked. 
    % Only relevant if the FFT length is even
    f = [f,fs/2];
    % Add the Nyquist frequency to the frequency array
    X1 = [X1,X(nfft/2+1)];
    % Add the value of the spectrum at the Nyquist frequency
    % We do not double the power of the spectral component at the Nyquist
    % frequency.
end

end