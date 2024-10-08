function sig = pol_rotator(sigin,theta)
% Polarization rotator
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This functions rotates the input optical field by a given angle in the
% counter clockwise direction, i.e. the rotation matrix:
% R(theta)=[cos(theta) -sin(theta)  
%           sin(theta) cos(theta)]
% is applied to the input signal.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = pol_rotator(sig,pi/2);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% theta             rotation angle, in radians [real vector]
%
%                       For a static rotator, theta is a scalar
%
%                       For a dynamic rotator, theta is an 1xN vector,
%                       where N is the number of samples
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

sig.x = cos(theta).*sigin.x - sin(theta).*sigin.y;
sig.y = sin(theta).*sigin.x + cos(theta).*sigin.y;

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------