function Sav = calc_stokes_integrate(S,nav)
% Function template
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements ...
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       This is it.
%
% params            model parameters [structure]
%
%                       params.field1       
%                           this is it [real scalar]
% 
%                       params.field2 
%                           this is it [real scalar]
%
%                           params.field2 = 'value'
%                               in this case the model does this
%   
%                   possible formats
%                       [real vector], [real scalar], [integer vector]
%                       [integer scalar], [structure], [complex vector]
%                       [optical signal structure], [binary vector] 
%                       [string]
%                       [0/1]
%                       [complex 3D array]
%                       [complex 2D array]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%                       This is it.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% CONSTANT              essential physical constants [structure]
%
% time_array            time samples, in s [real vector]
%
% dt                    time samples separation, in s [real scalar]
%
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% df                    frequency samples separation, in Hz [real scalar]
%
% space_grid            2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[~,nsamples] = size(S)




averaging_window = zeros(1,nsamples);
averaging_window(1:nav) = 1/nav;
% define window corresponding to the integration time; to be used for
% moving window average.

% Y = fft(X,n,dim) returns the Fourier transform along the dimension dim. For example, if X is a matrix, then fft(X,n,2) returns the n-point Fourier transform of each row.



Sav(1,:) = real(ifft(fft(S(1,;)).*fft(averaging_window)));
Sav(2,:) = real(ifft(fft(S(2,:)).*fft(averaging_window)));
Sav(3,:) = real(ifft(fft(S(3,:)).*fft(averaging_window)));
Sav(4,:) = real(ifft(fft(S(4,:)).*fft(averaging_window)));
% moving window averging over the specified integration time 



end