
% line_color = linspecer(6,'qualitative');
% Use plot(x,y,'Color',line_color(1,:));

do_print = 0;
do_add_figsize_to_filename = 1;
margin_figure = 0;

fig.interpreter = 'latex';

fig_name = [file_name_core_figure '_'];
hfig = figure('Name',fig_name);
plot(x,y,'Color','b','LineStyle','-','Marker','none','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor','w','HandleVisibility','on')
xlabel(' ','Interpreter',fig.interpreter)
ylabel(' ','Interpreter',fig.interpreter)
legend('','Location','NorthEast','Box','off','Interpreter',fig.interpreter,'FontSize',18)
% title('','FontWeight','Normal','Interpreter',fig.interpreter)
ax = gca;
ax.TickLabelInterpreter = fig.interpreter;
ax.FontSize = 18;
ax.LineWidth = 1.5;
ax.TickLength = [0.01 0.025]*1.5;
ax.XMinorTick = 'on';
% xlim([0 18])
% ylim([-0.5 0.05])

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

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




hfig.PaperUnits = 'centimeters';
hfig.PaperType = 'a4';
hfig.PaperOrientation = 'landscape';
hfig.PaperPosition = [5 5 17 10];