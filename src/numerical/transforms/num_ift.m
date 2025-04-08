function x = num_ift(X,dt,t0)
% Numerical estimation of the inverse continuous Fourier transform based on IDFT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function estimates the inverse Fourier transform x(t) of a function 
% X(w) based on its samples X[k] = X(2*pi*(k -1)*df), k = 1...N
% The samples of the inverse Fourier transform are 
% x[l] = x((l - 1)*dt), l = 1...N
% The input spectrum samples are provided in increasing frequency order.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% f = num_ift(F,dt);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% X                 samples of the Fourier transform of x(t), in increasing
%                       frequency order [complex vector]
%                       X[k] = X(2*pi*(k -1)*df), k = 1...N
%                       with df*dt = 1/N
%
% dt                samples separation [real scalar]
%
% t0                centre of the integration window [real scalar]
%                       See description of num_ft function.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% x                 samples of the time-domain function x(t)
%                       [complex vector]
%                       x[l] = x((l - 1)*dt), l = 1...N
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

switch nargin
    
    case 2
        x = ifft(fftshift(X))/dt;
        % We keep this special case only for backward compatibility, even
        % though it is incorporated in the more general case below.
        % The function was first implemented by assuming the [0 T] interval. 
        
    case 3
        w = 2*pi*(0:length(X) - 1)/length(X)/dt;
        % Frequencies, in FFT order.        
        x = ifft(fftshift(X).*(-1).^(0:length(X) - 1).*exp(1j*w*t0))/dt;
        % OBS: fftshift same as ifftshift for even number of samples.
        
    otherwise
        
        error('num_ift: number of input parameters should be 2 or 3.');

end