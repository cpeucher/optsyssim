%--------------------------------------------------------------------------
% Weights for Simpson's rule integral calculations
% From:
% Introduction to Numerical Methods and Matlab Programming for Engineers
% Todd Young & Martin Mohlenkamp
% Copyright © 2008, 2009, 2011 Todd Young & Martin Mohlenkamp
% Available: http://www.math.ohiou.edu/courses/math344/ [Retrieved
% 2011-12-12]
% Available: http://www.ohiouniversityfaculty.com/youngt/IntNumMeth/
% [Retrieved 2021-11-24]
%--------------------------------------------------------------------------
function w = mysimpweights(n)
if rem(n,2) ~= 0
    error('n must be even')
end
w = ones(n+1,1);
for i = 2:n
    if rem(i,2)==0
        w(i)=4;
    else
        w(i)=2;
    end
end


