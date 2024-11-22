function y = dsp_overlap_save(x,h,L)
% Digital filtering in the frequency domain using overlap-and-save method
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements digital filtering in the frequency domain using
% the overlap-save method, which is particularly well suited for
% processing long streams of data in real time (not here...)
% See e.g.
% J. G. Proakis and D. G. Manolakis, Digital signal processing, Fourth 
% edition, Pearson new international edition. Harlow, Essex: Pearson, 2014.
% chapter 7, pp. 498-...
% Here we have a simulation / off-line implementation where the entire data
% to process is available when calling the function, i.e. the segmentation
% is performed within the function.
% The present implementation is not vectorised, but loops over the input
% blocks. A vectorized implementation may require large memory in order to
% build a matrix of blocks and filter them in the frequency domain by
% fft/ifft.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% y = dsp_overlap_save(x,h,L);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input signal [vector]  
%
% h                 filter tap coefficients [vector]
%
% L                 block length in the overlap-save scheme [integer]
%                       For 50% overlap use L = M - 1, where M is the
%                       filter length M = length(h).
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% y                 output signal [vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

x = x(:).';
% Ensure x is a line vector

h = h(:).';
% Ensure h is a line vector

M = length(h);
% FIR filter length

nblocks = floor(length(x)/L);
% Number of blocks to process

N = L + M - 1;
% Size of the FFT/IFFT

xapp = zeros(1,M - 1);
% Prepare block of M - 1 zeros to append at the start of the data block in
% the first iteration

H = fft([h,zeros(1,N - M)]);
% Zero padding of the filter impulse response so that it has N elements,
% and calculation of its DFT
% This is done once for all.

y = [];
% Create empty vector to store the output


for iblock = 1:nblocks
    
    xx = x((iblock - 1)*L + 1:iblock*L);
    % Extract a block of length L from the input data
    
    xb = [xapp,xx];
    % Data block to process
    % The length of xb is N = L + M - 1.
    
    xapp = xx(L - M + 2:L);
    % Prepare the overlap to append in the next iteration
    % These are the last M - 1 elements from the current block.
    
    yb = ifft(fft(xb).*H);
    % Filter output calculated in the frequency domain for the block of
    % length L + M - 1
    
    y = [y,yb(M:end)];
    % We discard the first M -1 samples and update the output by
    % concatenation with the previously calculated samples.
    
end
% End of loop over blocks



end