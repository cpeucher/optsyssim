function x = prepad (x, L)
%POSTPAD   Pads or truncates a vector x to a specified length L.
%
%   This is a simpler version of the file distributed with octave.
%   The function does not check input parameters.

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHOR : Peter Soendergaard.
%   TESTING: OK
%   REFERENCE: NA

xl=size(x,1);

if xl==L
  % Do nothing, bails out as quickly as possible.
  return;
end;

xw=size(x,2);

if ndims(x)>2
  error('Postpad of multidim not done yet.');
end;

if xl<L
  x=[zeros(L-xl,xw); x];
else
  %x=x(1:L,:);
  x=x(xl-L+1:xl,:);
end;
