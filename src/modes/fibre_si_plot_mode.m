function fibre_si_plot_mode(field,params,visparams)
% Plot field or intensity distribution (irradiance pattern) of fibre mode
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function plots the field or intensity distribution of a fibre mode.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% visparams.limit_radius = Inf;% 2*params_fibre.a;
% visparams.show_core_limit = 0;
% visparams.save = 0;
% visparams.colormap = 'jet';%'hot';
% visparams.name = 'Mode field distribution';
% fibre_si_plot_mode(field,params_fibre,visparams);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% field             quantity to be plotted [real matrix]
%
%                       This could be the mode field or mode power
%                       distibution.
%
%                       The field is defined over the space_grid Cartesian 
%                       grid.
% 
% params            optical fibre parameters [structure]
%
%                       params.a: core radius, in m [real scalar]
%
% visparams         plotting options [structure]
%
%                       visparams.limit_radius
%                            radius over which the field will be
%                            represented [real scalar]
%
%                            This allows to limit the representation to a
%                            circular area, for instance corresponding to
%                            the cladding dimension, instead of
%                            representing the mode on the entire
%                            rectangular (square) space_grid.
%
%                           visparams.limit_radius = params.a;
%                               will only represent the field in the core.
%
%                           visparams.limit_radius = Inf; 
%                               will represent the field over the entire 
%                               space_grid
%
%                           visparams.show_core_limit
%                               displays the limit of the core as a dashed 
%                               white line [0/1]
%
%                           visparams.save
%                               saves the output in .emf and .jpg formats
%                               [0/1]
%
%                           visparams.colormap
%                               color map used to display the field
%                               [string]
%  
%                               e.g. visparam.colormap = 'hot';'jet';% etc
% 
%                           visparams.name
%                               name of the figure [string]%                               
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% space_grid        2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global space_grid

vis_range = double(space_grid.X.^2 + space_grid.Y.^2 < visparams.limit_radius^2);
vis_range(vis_range == 0) = NaN;
% Set radius of visualisation. Only the space within this radius will be
% represented.

figure('Name',visparams.name)
pcolor(space_grid.X,space_grid.Y,field.*vis_range);
hold on

if visparams.show_core_limit
    plot_core_limit_x = linspace(-params.a,params.a,1001);
    plot(plot_core_limit_x,sqrt(params.a^2-plot_core_limit_x.^2),'LineWidth',5,'Color','w','LineStyle',':');
    plot(plot_core_limit_x,-sqrt(params.a^2-plot_core_limit_x.^2),'LineWidth',5,'Color','w','LineStyle',':');
end

axis square
shading flat
colormap(visparams.colormap)
axis off

if visparams.save
    print(gcf,'-dmeta','-r72',[visparams.name '.emf']);
    print(gcf,'-djpeg','-r72',[visparams.name '.jpg']);
end

end