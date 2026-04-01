function print_figure(hfig,do_print,do_add_figsize_to_filename,margin_figure)
% Print figure according to our very own needs and specifications
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function prints a figure following our very own specifications for
% e.g. notes or lecture notes.
% The point is to remove clutter from scripts.
% Can be made better.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% do_print = 1;
% do_add_figsize_to_filename = 0;
% margin_figure = 0;
% print_figure(hfig,do_print,do_add_figsize_to_filename,margin_figure)
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% hfig              figure handle [string]
%
% do_print          printing switch [0/1]
%
% do_add_figsize_to_filename       as it says switch [0/1]
%
% margin_figure     switch for adjusting dimensions of the figure (higher
%                       than wide) so that it can fit in the margin of some
%                       document templates (e.g. kaobook)
%   
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if margin_figure
    hfig.Position(3:4) = [1 0.9]*hfig.Position(4);
end

if ~do_print
    return
else

    if do_add_figsize_to_filename
        [figure_ratio_numerator,figure_ratio_denominator] = rat(hfig.Position(3)/hfig.Position(4));
        file_name = [hfig.Name '_w' num2str(figure_ratio_numerator,'%u') 'h' num2str(figure_ratio_denominator,'%u')];
    else
        file_name = hfig.Name;
    end

    print(file_name,'-dmeta');
    print(file_name,'-djpeg');
    print(file_name,'-dpdf');
    crop_command =['pdfcrop ' file_name '.pdf ' file_name '.pdf'];
    system(crop_command);
end


end