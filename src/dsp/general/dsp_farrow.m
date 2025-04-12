function y = dsp_farrow(x,c,binit,mu)
% Implementation of Farrow structure, e.g. for interpolation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a Farrow structure made of N FIR filter
% consisting each of M taps.
% The structure is as follows:
%
%                       in ----------------------------------------
%                           |      |                |      |      |
%                           C_N   C_{N -1}          C3     C2     C1
%                           |      |                |      |      |
%                           |      |                |      |      |                           
%                           ---X------X--        ------X------X------> out
%                              |      |                |      |
%                       mu ------------------------------------
% 
% The C_j, j = 1... N represent the M-tap FIR filters.
% N is the number of filters.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% Linear interpolator (2 basepoints - 1 sample delay):
% nfilters = 2;
% ntaps = 2;
% c0 = [0, 1];
% c1 = [1, -1];
% c = [c0;c1].';
%
% Parabolic interpolator (3 basepoints - 1 sample delay):
% nfilters = 3;
% ntaps = 3;
% c0 = [0,   1,  0];
% c1 = [0.5, 0,  -0.5];
% c2 = [0.5, -1, 0.5];
% c = [c0;c1;c2].';
%
% Cubic interpolator (4 basepoints - 2 sample delay): 
% nfilters = 4;
% ntaps = 4;
% c0 = [0,    0,     1,    0];
% c1 = [-1/6, 1,     -1/2, -1/3];
% c2 = [0,   1/2,    -1,   1/2];
% c3 = [1/6, -1/2,   1/2,  -1/6]; 
% c = [c0;c1;c2;c3].';
%
% Parabolic interpolator (4 basepoints - 2 sample delay)
% ntaps = 4;
% nfilters = 3;
% alpha = 0.5;
% c = [0,  -alpha,     alpha;
%     0,   alpha+1,  -alpha;
%     1,   alpha-1   -alpha;
%     0,  -alpha,     alpha];
%
% y = dsp_farrow(x,c,zeros(ntaps,nfilters),mu);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% x                 input samples [vector]
%
% c                 coefficients of the FIR filters in the structure 
%                       [matrix] 
%
%                       c is a matrix with ntaps lines and nfilters
%                       columns, where ntaps is the number of taps in each
%                       filter and nfilters is the number of FIR filters.
%                       The tap coefficient are ordered as required by the
%                       dsp_fir_linear function.
%                       c(:,j) contains the coefficients of filter C_j for 
%                       j= 1...N.
%
% binit             initial values of the buffer for the FIR filters 
%                       [matrix]
%
%                       The size of binit is the same as the size of c.
%                       binit(:,j) contains the initial values of the
%                       buffer in filter C_j, for j = 1...N
%
% mu                fractional delay [real scalar or vector] 
%
%                       If a vector, its length should be the same as the
%                       length of the input samples vector x.
%                       The fractional delay takes values within 
%                       0 <= mu < 1.
%                       It defines the interpolants in terms of fraction of
%                       the interval between x[nl] and x[nl+1], normalized
%                       by the input samples sampling interval Tin.
%                       The position of the interpolant with basepoint
%                       index nl and fractional delay mu[nl] is
%                       t[nl] = (nl + mu[nl])*Tin;
%                       The vector mu therefore contains as many elements
%                       as the vector x. However, depending on the degree
%                       size of the FIR filters in the structure, the
%                       interpolation may not be correct for the first
%                       values of x[nl], mu[nl]. 
%                       If mu is a scalar, the same value of fractional
%                       interval is applied to all samples. 
%                       mu = 0.5 and mu = 0.5*ones(1,length(x)) are
%                       therefore equivalent.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% y                 output samples [vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if length(mu) > 1 && length(mu)~=length(x)
    error('dsp_farrow: the length of mu should be one or length(x)');
end

[ntaps,nfilters] = size(c);
% Determine the number of filters

if length(mu) > 1
    mu = [zeros(1,floor(ntaps/2)),mu];
    mu = mu(1:length(x));
end
% Ensure correspondance between the fractional delay and the sample index
% due to initialisation of FIR filters


y = dsp_fir_linear(x,c(:,nfilters),binit(:,nfilters));
% Initialise the Farrow structure output with the output of filter
% #nfilters
% Would work fine in case only one filter is present



for ifilters = 1:nfilters-1
    % Horner calculation of the output.    
    y = dsp_fir_linear(x,c(:,nfilters - ifilters),binit(:,nfilters - ifilters )) + mu.*y;    
end
% Loop from filter nfilters - 1+1 to filter 1.
   

end