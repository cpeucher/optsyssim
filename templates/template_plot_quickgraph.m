% -------------------------------------------------------------------------
% 
%
%
% 2024-08-20
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Clean-up
% -------------------------------------------------------------------------
clear all
close all

% -------------------------------------------------------------------------
% Switches
% -------------------------------------------------------------------------
do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;


% -------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------
% data_path = '../data/';
% 
% data_fname = 'filename.mat';
% load([data_path data_fname]);
% % Load Matlab mat file
 
% data_fname = 'filename.dat';
% data_format = '%f %f %f';
% data_size =[3 Inf]; % Number of columns.
% data_fid = fopen([data_path data_fname],'r');
% fgets(data_fid); % Ignore first line.
% data = fscanf(data_fid,data_format,data_size);
% data = data';
% fclose(data_fid);
% % Load text file

% data_fid = fopen([data_path data_fname]);
% C = textscan(data_fid,'%u %s %f %f %f %f %f %f','delimiter',',','headerlines',1,'collectoutput',0);
% fclose(data_fid);




%%
% -------------------------------------------------------------------------
% Fonts etc.
% -------------------------------------------------------------------------
fig.font_name = 'times';
fig.font_size = 18;
fig.interpreter = 'latex';
mred = [178 0 24]/255;

% -------------------------------------------------------------------------
% Choice of colormap
% -------------------------------------------------------------------------
% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));

% -------------------------------------------------------------------------
% File name
% -------------------------------------------------------------------------
file_name_core_figure = strrep(mfilename,'plot','fig');
file_name_core_data = strrep(mfilename,'plot','data');
time_stamp = string(datetime('now','TimeZone','Z','Format','yyyyMMdd''T''HHmmss''Z'));


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Processing
% -------------------------------------------------------------------------

x = [0:0.01:2]*2*pi;
y1 = sin(x);
y2 = cos(x);



%%
% -------------------------------------------------------------------------
% Figure: 
% -------------------------------------------------------------------------
fig_name = [file_name_core_figure];
hfig = figure('Name',fig_name);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
% yyaxis left
h1 = plot(x,y1,'LineWidth',1.5,'Color','b','LineStyle','-','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');
hold on
h2 = plot(x,y2,'LineWidth',1.5,'Color','r','LineStyle','--','Marker','none','MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on');

% -------------------------------------------------------------------------
% Legend
% -------------------------------------------------------------------------
hleg = legend([h1,h2],{'',''});

hleg.Title.String = '';
hleg.Location = 'SouthEast';
hleg.Orientation = 'vertical';
hleg.NumColumns = 1;

hleg.Box = 'off';
hleg.EdgeColor = 'black';
hleg.TextColor = 'black';
hleg.Color = 'white';
hleg.LineWidth = 1;

hleg.Interpreter = fig.interpreter;
hleg.FontName = fig.font_name;
hleg.FontSize = fig.font_size;

hleg.Visible = 'on';


% -------------------------------------------------------------------------
% Axes
% -------------------------------------------------------------------------
ax = gca;

ax.LineWidth = 1.5;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.TickLength = [0.01 0.025]*1.5;
ax.TickLabelInterpreter = fig.interpreter;

ax.XAxis.FontName = fig.font_name;
ax.XAxis.FontSize = fig.font_size;
ax.XAxis.Color = 'k';
% xlim([]);
% ax.XTick = [ ];
% xtickformat(ax,'%3.3f')
% ax.XTickLabel = {'','',''}

% yyaxis left
ax.YAxis(1).FontName = fig.font_name;
ax.YAxis(1).FontSize = fig.font_size;
ax.YAxis(1).Color = 'k';
% ylim([ ]);
% ax.YTick = [ ];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

% yyaxis right
% ax.YAxis(2).FontName = fig.font_name;
% ax.YAxis(2).FontSize = fig.font_size;
% ax.YAxis(2).Color = 'k';
% ylim([ ]);
% ax.YTick = [ ];
% ytickformat(ax,'%3.3f')
% ax.YTickLabel = {'','',''}

% ax.Layer = 'top';
% Whenever needed, ensure that the axes are on top of e.g. color boxes.


% -------------------------------------------------------------------------
% Axes labels
% -------------------------------------------------------------------------
x1 = xlabel('');
x1.Interpreter = fig.interpreter;
x1.FontName = fig.font_name;
x1.FontSize = fig.font_size;
x1.FontWeight = 'normal';
x1.Color = 'k';

% yyaxis left
yl = ylabel('');
yl.Interpreter = fig.interpreter;
yl.FontName = fig.font_name;
yl.FontSize = fig.font_size;
yl.FontWeight = 'normal';
yl.Color = 'k';

% yyaxis right
% y2 = ylabel('');
% y2.Interpreter = fig.interpreter;
% y2.FontName = fig.font_name;
% y2.FontSize = fig.font_size;
% y2.FontWeight = 'normal';
% y2.Color = 'k';



% -------------------------------------------------------------------------
% Grid lines
% -------------------------------------------------------------------------
ax.XGrid = 'off';
ax.XMinorGrid = 'off';

ax.YGrid = 'off';
ax.YMinorGrid = 'off';

ax.GridColorMode = 'manual';
ax.GridColor = 'k';

ax.MinorGridColorMode = 'manual';
ax.MinorGridColor = 'k';

ax.GridAlphaMode = 'manual';
ax.GridAlpha = 0.25;

ax.MinorGridAlphaMode = 'manual';
ax.MinorGridAlpha = 0.25;

ax.GridLineStyle = '-';
ax.MinorGridLineStyle = '--';


% -------------------------------------------------------------------------
% Add text
% -------------------------------------------------------------------------
% tt1 = text(3,6,'text','Interpreter',fig.interpreter,'FontName',fig.font_name,'FontSize',fig.font_size,'HorizontalAlignment','left');

% tt = text(0.63,-0.5,'$x \simeq 0.6$');
% tt.Interpreter = fig.interpreter;
% tt.FontSize = fig.font_size;
% tt.FontName = fig.font_name;
% tt.FontWeight = 'normal';
% tt.Color = 'black';
% tt.HorizontalAlignment = 'left';


% -------------------------------------------------------------------------
% Horizontal or vertical lines
% -------------------------------------------------------------------------
% xline = 7.5;
% yy = get(gca,'ylim');
% plot([xline xline],yy,'k--','LineWidth',1);
% % Vertical line
% 
% yline = 0.4;
% xx = get(gca,'xlim');
% plot(xx,[yline yline],'k','Linewidth',1);
% % Horizontal line


% -------------------------------------------------------------------------
% Adjust aspect ratio of the figure, if needed.
% -------------------------------------------------------------------------
if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

% -------------------------------------------------------------------------
% print to file
% -------------------------------------------------------------------------
if do_add_figsize_to_filename
    [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
    fig_name = [fig_name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
end

if do_print
    print(fig_name,'-dmeta');  
    print(fig_name,'-djpeg');
    print(fig_name,'-dpdf');
    crop_command =['pdfcrop ' fig_name '.pdf ' fig_name '.pdf'];
    system(crop_command);
end