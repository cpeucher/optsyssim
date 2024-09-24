function [sig_through,sig_drop] = opt_mrr(sig_in,sig_add,params)
% Micro-ring resonator filter in add-drop configuration
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function applies a micro-ring resonator (MRR) in add-drop 
% configuration to optical signals. The response of the resonator is
% assumed to be identical for the 2 polarisations of the signal.
%
%  drop <-----------------------------------< add
%                       -      -
%                    -           -
%                  -               -
%                 -                 -
%                  -               -
%                    -           -
%                       -      -
%  input>-----------------------------------> through
%
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_mrr.field_round_trip_loss = 0.96;
% params_mrr.power_coupling_1 = 0.63;
% params_mrr.power_coupling_2 = params_mrr.power_coupling_1;
% params_mrr.centre_frequency = 0;
% params_mrr.fsr = 100e9;
% params_mrr.visualiser_status = 0;
% params_mrr.save_tf.status = 0;
% params_mrr.save_tf.file_name = 'mrr_tf.dat';
% sig_add = opt_no_sig;
% [sig_through,sig_drop] = opt_mrr(sig,sig_add,params_mrr); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig_in            optical signal at the "input" port of the MRR
%                       [optical signal structure]
%
% sig_add           optical signal at the "add" port of the MRR
%                       [optical signal structure]
%
% params            MRR parameters [structure]
%                       See the description of opt_tf_mrr for more details
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig_through       optical signal at the through port of the MRR
%                       [optical signal structure]
%
% sig_drop          optical signal at the drop port of the MRR
%                       [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array


tf = opt_tf_mrr(frequency_array,params,params);
% Calculate the transfer functions of the MRR

spectrum_in_x = fft(sig_in.x);
spectrum_in_y = fft(sig_in.y);
spectrum_add_x = fft(sig_add.x);
spectrum_add_y = fft(sig_add.y);
% Calculate spectra of the input signals

sig_through.x = ifft(spectrum_in_x.*fftshift(tf.s11)) + ifft(spectrum_add_x.*fftshift(tf.s12));
sig_through.y = ifft(spectrum_in_y.*fftshift(tf.s11)) + ifft(spectrum_add_y.*fftshift(tf.s12));
sig_drop.x = ifft(spectrum_in_x.*fftshift(tf.s21)) + ifft(spectrum_add_x.*fftshift(tf.s22));
sig_drop.y = ifft(spectrum_in_y.*fftshift(tf.s21)) + ifft(spectrum_add_y.*fftshift(tf.s22));
% Apply the MRR transfer funcitons to the input signals and revert to 
% time domain.

end