function R = jones_rotation(theta)
% Rotation matrix
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This functions constructs the 2x2 rotation matrix corresponding to a 
% given rotation angle defined in the counter-clockwise direction.
% The rotation angle can be time-dependent.
% The rotation matrix is
% R(theta) = [cos(theta)  sin(theta);  
%             -sin(theta) cos(theta)];
% 
% This function is to be used to convert Jones vectors and Jones matrices
% between two orthonormal basis.
% If (x,y) is the canonical orthonormal basis on which a Jones vector J and
% a Jones matrix T are expresses.
% If (x',y') is an orthonormal basis that is obtained from (x,y) by a
% rotation of angle theta. If J' and T' are the corresponding Jones vector
% and Jones matrix, respectively, expressed on the basis (x',y').
% Then:
% J' = R(theta) J
% J = R(-theta) J'
% T' = R(theta) T R(-theta)
% T = R(-theta) T' R(theta)
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% R = jones_rotation(theta); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% theta             rotation angle, in radians [real vector]
% 
%                       For a static rotation, theta is a scalar.
%
%                       For a dynamic rotation, theta is an 1xN vector,
%                           where N is the number of samples.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% R                 rotation matrix [2x2xN real matrix]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% REMARKS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% TO DO:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% CREDITS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% AUTHOR:
% -------------------------------------------------------------------------
% Christophe Peucheret (christophe.peucheret@univ-rennes.fr)
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

theta = theta(:).';
% Ensure angle is a line vector.

R(1,1,:) = cos(theta);
R(1,2,:) = sin(theta);
R(2,1,:) = -sin(theta);
R(2,2,:) = cos(theta);

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------