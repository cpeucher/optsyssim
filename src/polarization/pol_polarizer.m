function sig = pol_polarizer(sig,theta,il,per)
% polarizer
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a linear polarizer, defined by its orientation, to
% an optical signal. 
% Non-ideal behaviour of the polarizer, such as insertion loss and finite
% polarisation extinction ratio can also be accounted for.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% sig = pol_polarizer(sig,angle,il,per); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal, expressed on the canonical basis
%                       [optical signal structure]
%
% theta             angle of the polarizer with respect to the (0,x) 
%                       axis of the canonical basis (x,y), in rad
%                       [real vector]
%
%                       For a static polarizer, theta is a scalar
%
%                       For a dynamic polarizer, theta is an 1xN vector,
%                           where N is the number of samples
%
% il                insertion loss of the polarizer, in dB [real scalar]
%
% per               polarisation extinction ratio, in dB [real scalar]
%
%                       Use Inf for an ideal polarizer. 
%                       This quantity has to be positive.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal, expressed on the canonical basis
%                       [optical signal structure]
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

in = [sig.x; sig.y];
% Build Jones vector of the input signal.

% We calculate:
% eout = jones_rotation(theta)*jones_polarizer_x(il,per)*jones_rotation(-theta)*ein;

tmp1 = jones_prod_mv(jones_rotation(theta),in);
tmp2 = jones_prod_mv(jones_polarizer_x(il,per),tmp1);
out = jones_prod_mv(jones_rotation(-theta),tmp2);

sig.x = out(1,:);
sig.y = out(2,:);
% Convert from Jones vector to the usual optical field structure defined in
% the basis (x,y).

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------