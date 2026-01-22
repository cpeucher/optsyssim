function field = fibre_si_lp_field(params,l,b,V,sym)
% Calculate scalar field of the LPlm mode of a step index fibre
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the scalar electric field of the guided LPlm 
% mode in a step-index fibre with perfect circular geometry (valid under
% the weakly-guiding approximation).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% xrange = [-1, 1]*50e-6;  
% yrange = [-1, 1]*50e-6; 
% nxpoints = 2001;
% nypoints = 2001;
% global space_grid
% space_grid = create_space_grid(xrange,yrange,nxpoints,nypoints);
% [space_grid.THETA,space_grid.RHO] = cart2pol(space_grid.X,space_grid.Y); 
% 
% params_fibre.n1 =  1.45;
% params_fibre.n2 =  1.4361;
% params_fibre.a = 25e-6;
% mode_type = 'LP';
% mode_l = 1;             % Index l of mode LPlm
% mode_m = 1;             % Index m of mode LPlm. Check that mode_m <= nmodes
% lambda = 850e-9; 
% V = 2*pi*params_fibre.a*sqrt(params_fibre.n1^2 - params_fibre.n2^2)/lambda;
% [b,nmodes] = fibre_si_b('LP',mode_l,V,params_fibre);
% % Normalised propagation constant
% field = fibre_si_lp_field(params_fibre,mode_l,b(mode_m),V,'even'); 
% % Mode field
% 
% visparams.limit_radius = Inf;
% visparams.show_core_limit = 0;
% visparams.save = 0;
% visparams.colormap = 'jet';%'hot';
% visparams.name = ['Mode field distribution for ' mode_type num2str(mode_l) num2str(mode_m)];
% fibre_si_plot_mode(field,params_fibre,visparams);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            optical fibre parameters [structure]
%
%                       params.n2: refractive index of the cladding 
%                           [real scalar]
%
%                       params.n1: refractive index of the core
%                           [real scalar]
%
%                       params.a: core radius, in m [real scalar]
% 
% l                 azimutal mode order [positive integer]
%
% b                 normalised propagation constant, for the mode LPlm as 
%                       obtained by calling 
%                       [b,nmodes] = fibre_si_b(mode,l,V,params) 
%                       [integer scalar]
%
% V                 normalised frequency [real scalar]
%
% sym               mode symmetry [string]
%
%                       sym = 'even';     for an even mode
%                       sym = 'odd';      for an odd mode
%
%                       Use 'even' for modes of the type LP0m.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% field             modal field [real matrix]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% space_grid        2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global space_grid



core = space_grid.X.^2 + space_grid.Y.^2 <= params.a^2;
% Signature of the core of the fibre
% Matrix whose elements are 1 inside the core, 0 in the cladding
clad = 1 - core;
% Signature of the cladding
% Matrix whose elements are 0 inside the core, 1 in the cladding

u = V*sqrt(1 - b);
v = V*sqrt(b);
% Transverse propagation constants u and v
% We have u^2 + v^2 = V^2.

fcore = besselj(l,u*space_grid.RHO/params.a);
% Radial component of the field in the core
fclad = besselj(l,u)*besselk(l,v*space_grid.RHO/params.a)/besselk(l,v);
% Radial component of the field in the cladding

fclad(abs(space_grid.RHO) < eps) = 0;
% Since we use the expression of the field in the cladding over the entire
% space grid, then multiply by the signature of the cladding clad, we will
% have some Inf * 0 = NaN indetermination if R = 0 is a point of the grid,
% due to the divergence of the besselk function.
% We need to anticipate that by letting fclad(R = 0)= 0 or something else
% than Inf.

field = fcore.*core + fclad .*clad;
% Total field (core + cladding).


switch sym
    % Distinguish the cases of odd and even modes.
    
    case 'even'
        
        field = field.*cos(l*space_grid.THETA);
        % Multiply by the azimutal component of the field for even mode       
        
    case 'odd'
        
        field = field.*sin(l*space_grid.THETA);
        % Multiply by the azimutal component of the field for odd mode  
        
    otherwise
        
        error('fibre_si_lp_field: mode should be even or odd. Use even for LP0m.')
        
end