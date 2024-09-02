function sig = pol_retarder(sig,azimuth,retardation)
% Optical wave retarder
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a wave retarder whose slow axis
% makes an angle azimuth with the -x axis of the canonical basis (x,y)
% and whose phase shift between the slow and fast components (also known as 
% retardation) is equal to retardation.
% In particular:
% for retardation = 2*pi   ->   wave plate
% for retardation = pi     ->   half wave plate
% for retardation = pi/2   ->   quarter wave plate
% The azimuth angle can change with time.
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% azimuth = 0;
% retardation = pi;% pi for half-wave plate, pi/2 for quarter-wave plate
% sig = pol_retarder(sig,azimuth,retardation);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal, whose components are expressed 
%                       in the canonical basis (x,y)
%                       [optical signal structure]
%
% azimuth           angle between the slow axis of the retarder and the 
%                       -x axis of the canonical basis [real vector]
%
%                       For a static retarder, azimuth can be a scalar
%
% retardation       phase retardation between the slow and fast axis 
%                       [real scalar]
%
%                       In a real implementation of a wave plate: 
%                       retardation=2*pi/lambda (n_s -n_f)*d where n_s and 
%                       n_f are the refractive indices associated with the 
%                       propagation of the slow and fast components,
%                       respectively, and d is the thickness of the plate.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal, whose components are expressed 
%                       in the canonical basis (x,y)
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

tmp1 = jones_prod_mv(jones_rotation(azimuth),in);
tmp2 = jones_prod_mv(jones_retarder(retardation),tmp1);
out = jones_prod_mv(jones_rotation(-azimuth),tmp2);
% Apply change of base and apply wave retarder.
% We calculate:
% eout = jones_change_coordinate(-azimuth)*jones_retarder(retardation)*jones_change_coordinate(azimuth)*ein;

sig.x = out(1,:);
sig.y = out(2,:);
% Convert from Jones vector to the usual optical field structure defined in
% the basis (x,y).

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------