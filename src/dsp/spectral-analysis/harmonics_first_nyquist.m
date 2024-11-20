function [freq,cf,cfs] = harmonics_first_nyquist(f,fs,nmax)
% Calculate the frequencies of aliases of harmonics in 1st Nyquist zone
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the frequencies of the aliases of the harmonics
% in the first Nyquist zone of a signal of frequency f sampled at fs.
% Since the use of the function is to calculate the total harmonic
% distortion of a sampled signal, it also returns the frequencies of the
% harmonics (non aliased) that fall into the first Nyquist zone. However,
% it does not return the frequencies of the aliased signal if f > fs/2. 
% The function also returns the integer cs and cfs so that the returned
% frequencies are
% freq = cf*f + cfs*fs
% with cf*cfs < 0 (i.e. the signs of cf and cfs are opposite).
% It is thus possible to distinguish between the aliased and non-aliased
% harmonics by looking into the indices cf and cfs.
% If cfs = 0, the corresponding returned frequency is that of a non-aliased
% harmonic of the signal within the 1st Nyquist zone.
% As stated above, the function only returns the frequencies of harmonics
% and aliased harmonics, not of aliases of the signal when f > fs/2.
% Therefore cf > 1.
% If the frequencies of the aliases of the signal are of interest, the
% function can be modified by starting the loop over harmonic numbers from
% the value 1 (right now it is  for n = 2:nmax ... end)
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% f                 frequency of the signal [real scalar]
%
%                       The signal frequency is expected to be f < fs/2
%                   
% fs                sampling frequency of the signal, in the same unit as
%                       the one in which f is specified [real scalar]
%
% nmax              maximum order of the harmonic to consider [integer]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% freq              frequencies of the harmonics of the signal as well as 
%                       the aliased of the harmonics of the signal in the 
%                       first Nyquist zone, i.e. between 0 and fs/2
%                       [real vector]
%
% cf                corresponding value of cf so that freq = cf*f + cfs*fs
%                       [integer vector]
%
% cfs               corresponding value of cfs so that freq = cf*f + cfs*fs
%                       [integer vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

freq = [];
cf = [];
cfs = [];
% Initialise ouputs to empty vectors
% We will return
% freq = cf* f + cfs*fs

for n = 2:nmax
    % Loop over harmonics orders    
    
    nnyquist = ceil(2*n*f/fs);
    % Determine the Nyquist zone to which the harmonic belongs
    
    if mod(nnyquist,2)
        % Odd Nyquist zone
        % f = n*f - k*fs
        
        k = fix(nnyquist/2);
           
        freq = [freq,n*f - k*fs];
        cf = [cf,n];
        cfs = [cfs,-k];
        
    else
        % Even Nyquist zone
        % f = k*fs - n*f
        
        k = nnyquist/2;
        
        freq = [freq,k*fs - n*f];
        cf = [cf,-n];
        cfs = [cfs,k];
    end
end
% End of loop over harmonics orders        

end