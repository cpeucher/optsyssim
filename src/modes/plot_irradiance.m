function hfig = plot_irradiance(F,visparams)
% Plot intensity distribution
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function represents the distribution of a real quantity F in a
% Cartesian coordinate system.
% Optional slices at specific x and/or y values can also be plotted.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% visparams.display = '2d';% '2d_slicex','2d_slicey','2d_slicexy'
% visparams.name = 'Intensity distribution for xxx beam';
% visparams.xslice = 0;
% visparams.yslice = 10e-6;
% hfig = plot_irradiance(I,visparams); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% F                 quantity to represent [real matrix]
%
%                       F is a matrix defined on the grid specified by the
%                       global variable space_grid
%
%                       F should only contain real numbers.
%
% visparams         visualisation parameters [structure]
%
%                       visparams.display   
%                           display mode [string]
%
%                           visparams.display = '2d'
%                           visparams.display = '2d_slicex'
%                           visparams.display = '2d_slicey'
%                           visparams.display = '2d_slicexy'
%
%                       visparams.name      
%                           name of the figure [string]
%
%                       visparams.xslice    
%                           value of x for slice representation 
%                           (in '2d_slicex' and '2d_slicexy' modes)
%                           [real scalar]
%
%                       visparams.yslice
%                           value of y for slice representation
%                           (in '2d_slicey' and '2d_slicexy' modes
%                           [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% hfig                  figure handle [string]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% space_grid            2D space grid in the transverse plane [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global space_grid


switch visparams.display
    
    case '2d'
        % Single 2D plot
        
        hfig = figure('Name',visparams.name);
        pcolor(space_grid.x,space_grid.y,F);
        xlabel('x');
        ylabel('y');
        % colormap(jet);
        axis square
        shading flat
        colormap hot
        % axis off
        
        
    case '2d_slicex'
        % 2D plot + horizontal slice
        
        Yq = space_grid.Y;
        Xq = visparams.xslice*ones(length(space_grid.y),length(space_grid.x));
        Fx = interp2(space_grid.X,space_grid.Y,F,Xq,Yq);
        Fx = Fx(:,1);
        % Interpolation of the data for the specified x value
        
        
        hfig = figure('Name',visparams.name);
        subplot(1,2,1)
        pcolor(space_grid.x,space_grid.y,F);
        vline(visparams.xslice,'w');
        xlabel('x');
        ylabel('y');
        % colormap(jet);
        axis square
        shading flat
        colormap hot
        % axis off
        subplot(1,2,2)
        plot(Fx,space_grid.y)
        axis square
        ylim([space_grid.y(1) space_grid.y(end)])
        ylabel('y');
               
        
    case '2d_slicey'
        % 2D plot + vertical slice       
        
        Xq = space_grid.X;
        Yq = visparams.yslice*ones(length(space_grid.y),length(space_grid.x));
        Fy = interp2(space_grid.X,space_grid.Y,F,Xq,Yq);
        Fy = Fy(1,:);
        % Interpolation of the data for the specified y value
        
        
        hfig = figure('Name',visparams.name);
        subplot(2,1,1)
        pcolor(space_grid.x,space_grid.y,F);
        hline(visparams.yslice,'w');
        xlabel('x');
        ylabel('y');
        %colormap(jet);
        axis square
        shading flat
        colormap hot
        %axis off
        subplot(2,1,2)
        plot(space_grid.x,Fy)
        axis square
        xlim([space_grid.x(1) space_grid.x(end)]);
        xlabel('x');       
        
        
    case '2d_slicexy'
        % 2D plot + horizontal and vertical slices
        
        Yq = space_grid.Y;
        Xq = visparams.xslice*ones(length(space_grid.y),length(space_grid.x));
        Fx = interp2(space_grid.X,space_grid.Y,F,Xq,Yq);
        Fx = Fx(:,1);
        % Interpolation of the data for the specified x value
        
        Xq = space_grid.X;
        Yq = visparams.yslice*ones(length(space_grid.y),length(space_grid.x));
        Fy = interp2(space_grid.X,space_grid.Y,F,Xq,Yq);
        Fy = Fy(1,:);
        % Interpolation of the data for the specified y value
        
        
        hfig = figure('Name',visparams.name);
        subplot(2,2,1)
        pcolor(space_grid.x,space_grid.y,F);
        hline(visparams.yslice,'w');
        vline(visparams.xslice,'w');
        xlabel('x');
        ylabel('y');
        %colormap(jet);
        axis square
        shading flat
        colormap hot
        %axis off
        subplot(2,2,2)
        plot(Fx,space_grid.y)
        axis square
        ylim([space_grid.y(1) space_grid.y(end)])
        ylabel('y');
        subplot(2,2,3)
        plot(space_grid.x,Fy)
        axis square
        xlim([space_grid.x(1) space_grid.x(end)]);
        xlabel('x');
        
        
    otherwise
        error('plot_irradiance: display format not defined.');
end
% End switch over display format


end