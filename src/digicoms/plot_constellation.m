function plot_constellation(symbs,type,figure_name,varargin)
% Quickly plot digital constellation diagram
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function plots a received constellation diagram with various
% standard representations.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% constellation_type = 'plain';%'heat';'cluster';
% constellation_name = 'received constellation';
% plot_constellation(symbs,constellation_type,constellation_name,[-5:1:5]);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% symbs             complex symbols to represent [complex vector], or
%                   clusters of complex symbols to represent [cell array]
%
%                       In the latter case symbs{ii,:} contains the
%                       elements (complex symbols) of the ii-th cluster 
%                       to represent. 
%
% type              type of visualisation [string]
%
%                       type = 'plain'      plain blue dots
%                       type = 'heat'       density map
%                       type = 'cluster'    one colour per cluster
%
% figure_name       name of the figure [string]
%
% varargin          tick positions [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

switch nargin
    case 4
        figure_ticks = varargin{1};
    otherwise
        figure_ticks = [-15:2:15];
end

figure('Name',figure_name)

switch type
    
    case 'plain'
        
        scatter(real(symbs),imag(symbs),'b','filled');
        
        max_iq = max(max(abs(real(symbs))),max(abs(imag(symbs))));
        
    case 'heat'
        
        scatplot(real(symbs), imag(symbs));
        % Makes us of 
        % Alex Sanchez (2021). scatplot 
        % https://www.mathworks.com/matlabcentral/fileexchange/8577-scatplot
        % MATLAB Central File Exchange. Retrieved October 8, 2021. 
        
        max_iq = max(max(abs(real(symbs))),max(abs(imag(symbs))));
        
    case 'cluster'
        
        hold on;
        max_iq = 0;
        for ii =1:numel(symbs)
            scatter(real(symbs{ii,:}),imag(symbs{ii,:}),'filled');            
            max_iq = max(max_iq,max(max(abs(real(symbs{ii,:}))),max(abs(imag(symbs{ii,:})))));            
        end       
        
    otherwise
        
        error('plot_constellation: visualisation type not defined.');        
        
end


xlim(1.1*[-1 1]*max_iq);
ylim(1.1*[-1 1]*max_iq);
axis square
xticks(figure_ticks);
yticks(figure_ticks);

end