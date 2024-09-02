function [sig_i, sig_j] = pol_pbs(sig,theta,il,per)
% Polarization beam splitter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function simulates a polarization beam splitter. The components of
% the signal along the -i and -j polarization axes, as defined by the angle
% of the polarization beam splitter with respect to the canonical -x axis
% are separated.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% [sig_x, sig_y] = pol_pbs(sig,theta,il,per); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
% theta             angle of the polarization beam splitter defined in 
%                       the counter-clockwise direction from the axis -x of 
%                       the canonical basis (x,y), in radians [real scalar]
%                       theta=(x^i)=(y^j)
%
% il                insertion loss of the PBS, in dB [real scalar]
%                       The loss is defined as the ratio of the sum of the
%                       power at the 2 outputs of the PBS to the power at
%                       the input. 
%                       This is therefore the "common loss" to the 2 output 
%                       ports.
%
% per               polarization extinction ratio, in dB [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig_i             output optical signal resolved along the axis -i of the
%                       PBS, which makes an angle theta with the axis -x 
%                       of the canonical basis (x,y). The signal is
%                       expressed in the canonical basis (x,y)
%                       [optical signal structure]
%
% sig_j             optical signal resolved along the axis -j of the 
%                       PBS, which makes an angle theta with the axis -y 
%                       of the canonical basis (x,y). The signal is 
%                       expressed in the canonical basis (x,y)
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

power_loss_linear = 10^(-il/10);
% Power loss in linear scale.

if per == Inf
    t_max = 1;
else
    epsilon = 10^(per/10);
    % polarization extinction ratio, in linear units.
    t_max = epsilon/(1 + epsilon);
end
% Calculate the maximum transmission (in power) that corresponds to the
% desired polarization.

in = [sig.x; sig.y];
% Build Jones vector of the input signal.

tmp1 = jones_prod_mv(jones_rotation(-theta),in);
tmp2 = jones_prod_mv(jones_polarizer_x(0,per),tmp1);
out = jones_prod_mv(jones_rotation(theta),tmp2);
% Jones matrix calculation of the output along -i, in the canonical
% basis (x,y).
% Loss is not included. But imperfect (non infinite) extinction ratio of
% the polarizer is.
% We calculate:
% eout = jones_rotation(theta)*jones_polarizer_x(0,per)*jones_rotation(-theta)*ein;

sig_i.x = sqrt(power_loss_linear*t_max)*out(1,:);
sig_i.y = sqrt(power_loss_linear*t_max)*out(2,:);
% Convert from Jones vector to the usual optical field strcuture defined in
% the basis (x,y).


tmp1 = jones_prod_mv(jones_rotation(-theta - pi/2), in);
tmp2 = jones_prod_mv(jones_polarizer_x(0,per),tmp1);
out = jones_prod_mv(jones_rotation(theta + pi/2),tmp2);
% Jones matrix calculation of the output along -j, in the canonical
% basis (x,y).
% Loss is not included. But imperfect (non infinite) extinction ratio of
% the polarizer is.
% We calculate:
% eout = jones_rotation(theta+pi/2)*jones_polarizer_x(0,per)*jones_rotation(-theta-pi/2)*ein;

sig_j.x = sqrt(power_loss_linear*t_max)*out(1,:);
sig_j.y = sqrt(power_loss_linear*t_max)*out(2,:);
% Convert from Jones vector to the usual optical field structure defined in
% the basis (x,y)


end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------