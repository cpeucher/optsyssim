%--------------------------------------------------------------------------
% Matrix of weights for Simpson's rule double integral calculations
% From:
% Introduction to Numerical Methods and Matlab Programming for Engineers
% Todd Young & Martin Mohlenkamp
% Copyright © 2008, 2009, 2011 Todd Young & Martin Mohlenkamp
% Available: http://www.math.ohiou.edu/courses/math344/ [Retrieved
% 2011-12-12]
% Available: http://www.ohiouniversityfaculty.com/youngt/IntNumMeth/
% [Retrieved 2021-11-24]
%--------------------------------------------------------------------------
function W = mydblsimpweights(m,n)
% Produces the m by n matrix of weights for Simpson's rule
% for double integrals
% Inputs: m -- number of intervals in the row direction.
%              must be even.
%         n -- number of intervals in the column direction.
%              must be even.
% Output: W -- a (m+1)x(n+1) matrix of the weights
if rem(m,2)~=0 | rem(n,2)~=0
    error('m and n must be even')
end
u = mysimpweights(m);
v = mysimpweights(n);
W = u*v';