function fighandle = plot_poincare_sphere(FigName)
% Generate Poincaré sphere to plot Stokes vector data
%
% This code is adapted from:
% Poincare Sphere Plot of Polarimetry Stokes Vectors
% Version 1.0.0.0 (255 KB) by Jeffrey Roth
% https://se.mathworks.com/matlabcentral/fileexchange/10979-poincare-sphere-plot-of-polarimetry-stokes-vectors
% [accessed 2024-08-02]
%
% We use the Poincare sphere representation.
% 
% See Copyright notice in license.txt file.
%
% The function plots the axes and the sphere itself.
% The Stokes data needs to be plotted with full control from outside the
% function by using:
% plot3(S(2,:)./S(1,:),S(3,:)./S(1,:),S(4,:)./S(1,:),'LineStyle','-','LineWidth',3,'Color','r','Marker','none','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r')
% for trajectories, or 
% plot3(S(2)/S(1),S(3)/S(1),S(4)/S(1),'Marker','o','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
% for points(states of polarization)

fighandle = figure('Name',FigName);
[X,Y,Z] = sphere(20);
X = X;
Y = Y;
Z = Z;
Hs = mesh(X,Y,Z,'facecolor','w','edgecolor',[0.5 0.5 0.5]);  % set grid facecolor to white
caxis([1.0 1.01]);  % set grid to appear like all one color
alpha(0.70);  % set opacity of sphere to 70%
axis equal;  % make the three axes equal so the ellipsoid looks like a sphere
set(gcf,'Renderer','opengl');
hold on;
% Draw x- and y- and z-axes
Hx = plot3([-1.5 1.5], [0 0], [0 0],'k-');
set(Hx,'linewidth',2,'linestyle','-','color','k');
ht_x = text(1.75,0,0,'+S_1','fontweight','bold','fontsize',12,'fontname','arial');
Hy = plot3([0 0], [-1.5 1.5], [0 0],'k-');
set(Hy,'linewidth',2,'linestyle','-','color','k');
ht_y = text(0.1,1.6,0,'+S_2','fontweight','bold','fontsize',12,'fontname','arial');
Hz = plot3([0 0], [0 0], [-1.5 1.5],'k-');
set(Hz,'linewidth',2,'linestyle','-','color','k');
ht_z = text(-0.05,0,1.35,'+S_3','fontweight','bold','fontsize',12,'fontname','arial');

ht_lcp = text(-0.05,0.0,1.1,'RCP','fontweight','bold','fontsize',12,'fontname','arial','color','k');
% Draw a bold circle about the equator (2*epsilon = 0)
x_e = (-1:.01:1);
for i = 1:length(x_e)
    z_e(i) = 0;
    y_e_p(i) = +sqrt(1 - x_e(i)^2);
    y_e_n(i) = -sqrt(1 - x_e(i)^2);
end
He = plot3(x_e,y_e_p,z_e,'k-',x_e,y_e_n,z_e,'k-');
set(He,'linewidth',2,'color','k');
% Draw a bold circle about the prime meridian (2*theta = 0, 180)
y_pm = (-1:.01:1);
for i = 1:length(x_e)
    x_pm(i) = 0;
    z_pm_p(i) = +sqrt(1 - y_pm(i)^2);
    z_pm_n(i) = -sqrt(1 - y_pm(i)^2);
end
Hpm = plot3(x_pm,y_pm,z_pm_p,'k-',x_pm,y_pm,z_pm_n,'k-');
set(Hpm,'linewidth',2,'color','k');
set(gca,'fontweight','bold','fontsize',12,'fontname','arial');
view(135,20);