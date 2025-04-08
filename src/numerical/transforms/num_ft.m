function X = num_ft(x,dt,t0)
% Numerical estimation of the continuous-time Fourier transform based on DFT
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function estimates the Fourier transform X(w) of a function x(t)
% based on its samples x[l] = x((l - 1)*dt), l = 1...N
% The samples of the spectrum X[k] = X(2*pi*(k -1)*df), k = 1...N are
% returned.
% The output spectrum samples are in increasing frequency order.
% We integrate
% X(w) = \int_{t0 - T/2}^{t0 + T/2} x(t) e^{-jwt} dt
% where T is the duration of the time window of the simulation.
% By default it is assumed that t0 = T/2.
% Otherwise, one needs to specify t0, e.g. t0 = 0.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% F = num_ft(f,dt);
% % For time interval [0 T]
% F = num_ft(f,dt,0);
% % For time interval [-T/2 T/2]
% F = num_ft(f,dt,(T1 + T2)/2);
% % For time interval [T1 T2]
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 samples of the original function  [complex vector]
%                       x[l] = x((l - 1)*dt), l = 1...N
%
% dt                samples separation [real scalar]
%
% t0                centre of the integration window [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% X                 samples of the Fourier transform of x(t), in increasing
%                       frequency order [complex vector]
%
%                       X[k] = X(2*pi*(k -1)*df), k = 1...N
%                       with df*dt = 1/N
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
        X = fftshift(fft(x))*dt;
        % We keep this special case only for backward compatibility, even
        % though it is incorporated in the more general case below.
        % The function was first implemented by assuming the [0 T] interval. 
        
    case 3
        % Calculation of the Fourier integral over the [t0 - T/2, t0 + T/2]
        % interval.
        w = 2*pi*(0:length(x) - 1)/length(x)/dt;
        X = fftshift((-1).^(0:length(x) - 1).*exp(-1j*w*t0).*fft(x))*dt;
        
    otherwise
        
        error('num_ft: number of input parameters should be 2 or 3.');
        
end