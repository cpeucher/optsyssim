function plot_constellation_mapping(constellation,params)
% Plot original digital constellation with mapping
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function plots the constellations of digital modulation formats as
% defined by the define_constellation.
% It is not to be used to plot noisy / distorted constellations, but just
% the theoretical constellations with the chose bit mapping.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params.name = ['Constellation for '];
% params.labels.dec = 1; 
% params.labels.bin = 1;
% params.axes = 0;
% params.save = 0;
% plot_constellation_mapping(constellation,params); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% constellation     constellation [complex vector] 
%                       As obtained by the dsp_constellation_define
%                       function.
%
% params            plot parameters  [structure]
%
%                      params.name               name of the figure
%                      params.labels.dec = 0,1   to plot decimal labels
%                      params.labels.bin = 0,1   to plot binary labels
%                      params.axes = 0,1         to plot axis
%                      params.save = 0,1         to save in emf format
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

symb = [0:length(constellation) - 1];
% Symbols indices for the given constellation
symb = dec2bin(symb);
% The symbols are now expressed in binary form


fig = figure('Name',params.name);
s = scatter(real(constellation),imag(constellation));
s.Marker = 'o';
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
s.SizeData = 100;
hold on
axis square
dx = 0.2;
dy = 0.2;
if params.labels.bin
    text(real(constellation) + dx, imag(constellation) + dy, cellstr(symb));
end
if params.labels.dec
    text(real(constellation) + dx, imag(constellation) - dy, num2str(bin2dec(cellstr(symb))),'Color','red');
end

xlim([min(real(constellation))-1 max(real(constellation))+1]);
ylim([min(imag(constellation))-1 max(imag(constellation))+1]);
xticks([min(real(constellation)):2:max(real(constellation))]);
yticks([min(imag(constellation)):2:max(imag(constellation))]);

if ~params.axes
    axis off
end

if params.save
    saveas(fig,['constellation_' params.name '.emf']);
end

end