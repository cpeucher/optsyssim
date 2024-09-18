function [sig_x_i,sig_x_q,sig_y_i,sig_y_q] = rx_coherent(sig,sig_lo,params)
% Front-end for polarisation- and phase-diversity coherent receiver
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a coherent receiver front-end made from two
% polarisation beam splitters, two 90-degree optical hybrids and 4 pairs of
% balanced photodiodes.
% The front-end is ideal in the sense that the optical hybrid circuit is 90
% degrees and all photodiodes are identical. No amplitude and phase
% imbalance is therefore accounted for in this model.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_coh.pbs_sig.angle = 0;
% params_coh.pbs_sig.loss = 0;
% params_coh.pbs_sig.extinction_ratio = Inf;
% params_coh.pbs_lo.angle = -pi/4;
% params_coh.pbs_lo.loss = 0;
% params_coh.pbs_lo.extinction_ratio = Inf;
% params_coh.pd.responsivity = 1;
% params_coh.pd.thermal_noise_density = 0*10e-12;
% params_coh.pd.shot_noise = 'off';
% params_coh.pd.dark_current = 0;
% params_coh.elpf.type = 'bessel';
% params_coh.elpf.order = 4;
% params_coh.elpf.f3dB = 0.7*symbol_rate;
% [sig_x_i,sig_x_q,sig_y_i,sig_y_q] = rx_coherent(sig,sig_lo,params_coh);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig              received optical signal [optical signal structure]
%                    
% sig_lo            local oscillator optical signal 
%                       [optical signal structure]
%
% params            receiver parameters [structure]
%
%                       params.pbs
%                           polarization beam splitter parameter
%                           [structure]
%
%                           See pol_pbs for more details.
%                           params.pbs.angle
%                           params.pbs.loss 
%                           params.pbs.extinction_ratio
%
%                       params.pd
%                           photodetectors parameters [structure]
%
%                           See rx_pin for more details.
%                           params.pd.responsivity
%                           params.pd.thermal_noise_density
%                           params.pd.shot_noise
%                           params.pd.dark_current
%
%                       params.elpf
%                           post-detection filter parameters [structure]
%
%                           See elec_elpf for more details.
%                           params.elpf.type
%                           params.elpf.order
%                           params.elpf.f3dB
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig_x_i           -x polarisation, in-phase current [real vector]
%
% sig_x_q           -x polarisation, quadrature current [real vector]
%
% sig_y_i           -y polarisation, in-phase current [real vector]
%
% sig_y_q           -y polarisation, quadrature current [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

[sig_x,sig_y] = pol_pbs(sig,params.pbs_sig.angle,params.pbs_sig.loss,params.pbs_sig.extinction_ratio);
% PBS applied to the signal
[sig_lo_x,sig_lo_y] = pol_pbs(sig_lo,params.pbs_lo.angle,params.pbs_lo.loss,params.pbs_lo.extinction_ratio);
% PBS applied to the local oscillator

fprintf('\n\n%s\n','Signal PBS');
fprintf('%s\t%3.2f%s\n','-x o/p: power //x:',mean(abs(sig_x.x).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-x o/p: power //y:',mean(abs(sig_x.y).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-y o/p: power //x:',mean(abs(sig_y.x).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-y o/p: power //y:',mean(abs(sig_y.y).^2)/1.0e-3,' mW');

fprintf('\n\n%s\n','LO PBS');
fprintf('%s\t%3.2f%s\n','-x o/p: power //x:',mean(abs(sig_lo_x.x).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-x o/p: power //y:',mean(abs(sig_lo_x.y).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-y o/p: power //x:',mean(abs(sig_lo_y.x).^2)/1.0e-3,' mW');
fprintf('%s\t%3.2f%s\n','-y o/p: power //y:',mean(abs(sig_lo_y.y).^2)/1.0e-3,' mW');


[sigoutx1,sigoutx2,sigoutx3,sigoutx4] = opt_hybrid90_2x4(sig_x,sig_lo_x);
% optical 90 degree 2x4 hybrid - x
[sigouty1,sigouty2,sigouty3,sigouty4] = opt_hybrid90_2x4(sig_y,sig_lo_y);
% optical 90 degree 2x4 hybrid - y

sig_x_i = rx_pin(sigoutx2,params) - rx_pin(sigoutx1,params);
sig_x_q = rx_pin(sigoutx4,params) - rx_pin(sigoutx3,params);
% Balanced photodiodes - x


sig_y_i = rx_pin(sigouty2,params) - rx_pin(sigouty1,params);
sig_y_q = rx_pin(sigouty4,params) - rx_pin(sigouty3,params);
% Balanced photodiodes - y

end