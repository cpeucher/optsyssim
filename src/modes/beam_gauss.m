function u = beam_gauss(pos,k,w0,z)
% Scalar field of a fundamental Gaussian beam
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the field of a paraxial Gaussian
% beam. In case of propagation, the envelope needs to be multiplied by a
% term in exp(-1i*k*z), in the rare cases when this matters.
% The scalar field u is normalised so that \iint |u|^2 dx dy = 1. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% position = [0 0];
% lambda = 1550e-9;
% k = 2*pi/lambda;
% w0 = 10e-6;
% z = 0;
% power = 1e-3;
% u = sqrt(power)*beam_gauss(position,k,w0,z);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% pos               position of the beam [2-element real vector]
%
%                       pos = [x0,y0]
%
% k                 wave vector in the medium [real scalar]
%
%                       k = 2*pi*n / lambda,
%
%                       where n is the refractive index and lambda the
%                       wavelength
%
% w0                beam waist, in m [real scalar]
%
% z                 axial distance, in m [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% u                 field of Gaussian beam, without the propagating term 
%                       in e^(-1i*k*z) [complex matrix]
%
%                       It is normalised so that
%                       \iint |u|^2 dx dy = 1
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% space_grid            2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global space_grid


zr = k*w0*w0/2;
% Rayleigh range (also known as confocal parameter)

w = w0 * sqrt(1 + (z/zr)^2);
% Spot size

R = z + zr*zr /z;
% Wavefront radius of curvature

r2 = (space_grid.X - pos(1)).^2 + (space_grid.Y - pos(2)).^2;

u = sqrt(2/pi)* exp(1i*atan(z/zr) - r2/w^2 -1i*k*r2/2/R)/w;
% Complex amplitude


end