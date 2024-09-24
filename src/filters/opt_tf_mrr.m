function tf = opt_tf_mrr(freq,params,vis)
% Transfer function of add-drop micro-ring resonator filter
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the elements of the scattering matrix of a 
% single ring resonator add-drop filter.
% See e.g. C. K. Madsen and J. H. Zhao, "Optical filter design and analysis 
% - A signal processing approach," Wiley, New York, 1999
% and more specifically Section 3.3.2 A single stage AR filter 
% (pp. 131-133).
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_mrr.field_round_trip_loss = 0.96;
% params_mrr.power_coupling_1 = 0.63;
% params_mrr.power_coupling_2 = params_mrr.power_coupling_1;
% params_mrr.centre_frequency = 0;
% params_mrr.fsr = 100e9;
% vis_mrr.visualiser_status = 1;
% vis_mrr.save_tf.status = 0;
% vis_mrr.save_tf.file_name = 'mrr_tf.dat';
% tf = opt_tf_mrr(frequency_array,params_mrr,vis_mrr); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              frequencies at which the transfer function 
%                       will be calculated, in Hz [real vector]
%
%                       This should be specified as relative frequencies
%                       with respect to the centre frequency of the
%                       simulation bandwidth (global parameter
%                       reference_frequency).
%
% params            MRR parameters [structure]
%
%                       params.field_round_trip_loss
%                           field loss for one round trip in the resonator.
%                           [real scalar]
%
%                           The couplers are assumed to be lossless. 
%
%                           params.field_round_trip_loss = 1 means that the 
%                               resonator is loss less
%
%                           params.field_round_trip_loss = 0.1 means that
%                               the power loss through one round trip in
%                               the resonator is equal to 0.1*0.1 = 0.01
%
%                       params.power_coupling_1
%                           power coupling factor to the cross port of the 
%                           input coupler (from the "in" port to the ring 
%                           or from the ring to the "through" port)
%                           [real scalar]
%
%                       params.power_coupling_2
%                           power coupling factor to the cross port of the
%                           through coupler (from the ring to the "drop" 
%                           port or from the "add" port to the ring)
%                           [real scalar]
%
%                       params.centre_frequency
%                           centre frequency of one of the resonances,
%                           in Hz [real scalar]
%
%                           This value is specified in terms of relative
%                           frequency with respect to the centre frequency
%                           of the simulation bandwidth (global parameter 
%                           reference_frequency)
%                           If the resonance is tuned to the centre of the 
%                           simulation bandwidth, i.e. to 
%                           reference_frequency, then set
%                           params.centre_frequency = 0
%
%                       params.fsr
%                           free spectral range of the resonator, in Hz
%                           [real scalar]
% 
% vis               transfer function saving and visualisation parameters
%                       [structure]
%
%                       vis.visualiser_status
%                           specifies whether the transfer functions are 
%                           displayed or not [0/1]
%
%                       vis.save_tf
%                           transfer function file saving information
%                           [structure]
%
%                               vis.save_tf.status = 0/1
%                                   specifies whether the transfer
%                                   functions should be saved in a text 
%                                   file [0/1]
%
%                               vis.save_tf.file_name
%                                   name of the file where the TF will be
%                                   saved [string]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                transfer functions between the different ports of the 
%                       add-drop ring resonator and calculated metrics
%                       [structure]
%
%                       tf.s11
%                           transfer function from "in" to "through"
%                           [complex vector]
%
%                       tf.s12
%                           transfer function from "add" to "through"
%                           [complex vector]
%
%                       tf.s21
%                           transfer function from "in" to "drop"
%                           [complex vector]
%       
%                       tf.s22
%                           transfer function from "add" to "drop"
%                           [complex vector]
%
%                       tf.finesse_analytical_exact
%                           finesse of the MRR (analytical, exact)
%                           [real scalar]
%
%                       tf.finesse_analytical_approximated
%                           finesse of the MRR (analytical, approximated)
%                           [real scalar]
%
%                       tf.q_analytical_exact
%                           Q-factor of the MRR centre resonance
%                           (analytical, exact) [real scalar]
%
%                       tf.q_analytical_approximated
%                           Q-factor of the MRR centre resonance 
%                           (analytical, approximated) [real scalar]
%
%                       tf.extinction_ratio_s21
%                           extinction ratio of the drop transfer function,
%                           in dB [real scalar]
%
%                       tf.extinction_ratio_s11
%                           extinction ratio of the through transfer 
%                           function, in dB [real scalar]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global reference_frequency


phi = 2*pi*(freq - params.centre_frequency)/params.fsr;
% Phase shift

z = exp(-1i*phi);
% Phase shift term

c1 = sqrt(1 - params.power_coupling_1);
c2 = sqrt(1 - params.power_coupling_2);
% Field transfer of the input couplers (bar)


% Calculation of the transfer functions:
tf.s11 = (c1 - c2*params.field_round_trip_loss.*z)./(1 - c1*c2*params.field_round_trip_loss.*z);
% "in" to "through"
tf.s12 = -sqrt(params.power_coupling_1*params.power_coupling_2*params.field_round_trip_loss.*z)./(1 - c1*c2*params.field_round_trip_loss.*z);
% "add" to "through"
tf.s21 = -sqrt(params.power_coupling_1*params.power_coupling_2*params.field_round_trip_loss.*z)./(1 - c1*c2*params.field_round_trip_loss.*z);
% "in" to "drop"
tf.s22 = (c2 - c1*params.field_round_trip_loss.*z)./(1 - c1*c2*params.field_round_trip_loss.*z);
% "add" to "drop"

%--------------------------------------------------------------------------
% Calculation of some resonator properties
%--------------------------------------------------------------------------
tf.finesse_analytical_exact = pi/2/asin((1 - params.field_round_trip_loss*c1*c2)/(2*sqrt(params.field_round_trip_loss*c1*c2)));
% MRR finesse (analytical, exact)
  
tf.finesse_analytical_approximated = pi*sqrt(params.field_round_trip_loss*c1*c2)/(1 - params.field_round_trip_loss*c1*c2);
% MRR finesse (analytical, approximated)

tf.q_analytical_exact = (reference_frequency + params.centre_frequency)/params.fsr*tf.finesse_analytical_exact;
% Q factor (analytical, exact)

tf.q_analytical_approximated = (reference_frequency + params.centre_frequency)/params.fsr*tf.finesse_analytical_approximated;
% Q factor (analytical, approximated)

tf.extinction_ratio_s21 = -10*log10((1 - params.field_round_trip_loss*c1*c2)^2/(1 + params.field_round_trip_loss*c1*c2)^2);
% Extinction ratio in-to-drop, in dB

tf.extinction_ratio_s11 = 10*log10(((c1 + c2*params.field_round_trip_loss)*(1 - params.field_round_trip_loss*c1*c2))^2/((c1 - c2*params.field_round_trip_loss)*(1 + params.field_round_trip_loss*c1*c2))^2);
% Extinction ratio in-to-through, in dB



%--------------------------------------------------------------------------
% Quick visualisation of the transfer functions
%--------------------------------------------------------------------------

if vis.visualiser_status == 1
    
    figure('Name','MRR Transfer function: linear scale');
    subplot(2,1,1)
    plot(freq./1e9,abs(tf.s11).^2,'');
    hold on;
    plot(freq./1e9,abs(tf.s12).^2,'k');
    plot(freq./1e9,abs(tf.s21).^2,'--r');
    plot(freq./1e9,abs(tf.s22).^2,'--g');
    legend('"in" to "through"','"add" to "through"','"in" to "drop"','"add" to "drop"');
    grid;
    xlabel('frequency (GHz)');
    ylabel('attenuation');
    subplot(2,1,2);
    plot(freq./1e9,unwrap(-angle(tf.s11)),'');
    hold on;
    plot(freq./1e9,unwrap(-angle(tf.s12)),'k');
    plot(freq./1e9,unwrap(-angle(tf.s21)),'--r');
    plot(freq./1e9,unwrap(-angle(tf.s22)),'--g');
    legend('"in" to "through"','"add" to "through"','"in" to "drop"','"add" to "drop"');    
    grid;
    xlabel('frequency (GHz)');
    ylabel('phase (rad)');
    % Plot transfer function on linear scale.


    figure('Name','MRR transfer function: logarithmic scale');
    subplot(2,1,1);
    plot(freq./1e9,10*log10(abs(tf.s11).^2),'');
    hold on;
    plot(freq./1e9,10*log10(abs(tf.s12).^2),'k');
    plot(freq./1e9,10*log10(abs(tf.s21).^2),'--r');
    plot(freq./1e9,10*log10(abs(tf.s22).^2),'--g');
    legend('"in" to "through"','"add" to "through"','"in" to "drop"','"add" to "drop"');
    grid;
    xlabel('frequency (GHz)');
    ylabel('attenuation (dB)');
    subplot(2,1,2);
    plot(freq./1e9,unwrap(-angle(tf.s11)),'');
    hold on;
    plot(freq./1e9,unwrap(-angle(tf.s12)),'k');
    plot(freq./1e9,unwrap(-angle(tf.s21)),'--r');
    plot(freq./1e9,unwrap(-angle(tf.s22)),'--g');
    legend('"in" to "through"','"add" to "through"','"in" to "drop"','"add" to "drop"');
    grid;
    xlabel('frequency (GHz)');
    ylabel('phase (rad)');
    % Plot Bode diagram (logarithmic scale).
    
end

%--------------------------------------------------------------------------
% Save transfer function into file
%--------------------------------------------------------------------------
if vis.save_tf.status == 1
    Res = [freq/1.0e9;10*log10(abs(tf.s11).^2);unwrap(-angle(tf.s11));10*log10(abs(tf.s12).^2);unwrap(-angle(tf.s12));10*log10(abs(tf.s21).^2);unwrap(-angle(tf.s21));10*log10(abs(tf.s22).^2);unwrap(-angle(tf.s22))];
    fid = fopen(vis.save_tf.file_name,'w');
    fprintf(fid,'%s %s %s %s %s %s %s %s %s\n','frequency','s11_mag','s11_phase','s12_mag','s12_phase','s21_mag','s21_phase','s22_mag','s22_phase');
    fprintf(fid,'%s %s %s %s %s %s %s %s %s\n','-','in->through','in->through','add->through','add->through','in->drop','in->drop','add->drop','add->drop');
    fprintf(fid,'%s %s %s %s %s %s %s %s %s\n','(GHz)','(dB)','(rad)','(dB)','(rad)','(dB)','(rad)','(dB)','(rad)');
    fprintf(fid,'%f %f %f %f %f %f %f %f %f\n',Res);
    fclose(fid);   
end


    
end
