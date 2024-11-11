function [phase,fgd,gd,fgvd,gvd] = extract_dispersion_from_tf(freq,tf,units,save)
% Group delay and dispersion from optical filter complex transfer function
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates group-delay and dispersion from a given complex
% transfer function. The values of the transfer function are assumed to be
% provided at equally spaced frequency samples corresponding to a subset of
% the simulation frequency grid. 
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% save_dispersion.status = 0;
% save_dispersion.file_name = 'filter_tf.dat';
% [phase,fgd,gd,fgvd,gvd] = extract_dispersion_from_tf(frequency_array,tf,'si',save_dispersion);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%                       This is it.
%
% freq              frequency array at which the TF samples are specified,
%                       in Hz [real vector]
%                       
%                       This has to be a subset of the global 
%                       frequency_array
%
% tf                complex transfer function of the filter 
%                       [complex vector]
% 
%                       The usual conventions H(f)=|H(f)| exp(-i \phi(f))
%                       are used.
%
% units             system of units used to return the group delay and 
%                       dispersion values [string]
%                       
%                       units = 'eng'
%                           The GD is specified in ps and the GVD in ps/nm
%           
%                       units = 'si'
%                            The GD is specified in s and the GVD in s/m
%
% save                 transfer function saving parameters [structure]
%
%                       save.status 
%                           specify whether to save the transfer function
%                           in a file or not [0/1]
%
%                       save.file_name
%                           name of the file where the TF will be saved
%                           [string]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% phase                 phase of the filter transfer function, in rad
%                           [real vector]
%
% fgd                   frequency values at which the group delay is 
%                           calculated, in Hz [real vector]
%                           
% gd                    group delay, in the system of units specified 
%                           according to to the choice in the input 'units' 
%                           variable [real vector]
%
% fgvd                  frequency values at which the dispersion is 
%                           calculated, in Hz [real vector]
%
% gvd                   dispersion, in the system of units specified
%                           according to the choice in the input 'units' 
%                           variable [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT              essential physical constants [structure]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% df                    frequency samples separation, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global reference_frequency 
global df
global CONSTANT

phase = -unwrap(angle(tf));
% Calculate the phase of the transfer function. It is assumed the unwrap
% function is sufficiently robust.

gd = diff(phase)./(df*2*pi);
% Group delay
fgd = freq(1:length(freq)-1) + df/2;
% Frequencies at which the group delay is evaluated by finite differences.
% This corresponds to midpoints between the original frequency samples.

fgvd = fgd(1:length(fgd) - 1) + df/2;
% Frequencies at which the GVD will be evaluated
gvd = -(reference_frequency + fgvd).^2.*diff(gd)/df/CONSTANT.c;
% Group velocity dispersion

if strcmp(units,'eng')
    gd = gd/1.0e-12;
    % Convert the group delay from s to ps
    gvd = gvd*1.0e3;
    % Convert the dispersion from s/m to ps/nm
elseif strcmp(units,'si')
    % Do nothing; the values are already in SI units
else
    error('extract_dispersion_from_tf: output unit system not implemented.');
end
% Convert to engineering units if necessary

if save.status == 1
    freq_save = freq(2:length(freq) - 1);
    tf_save = tf(2:length(tf) - 1);
    phase_save = phase(2:length(phase) - 1);
    freq_gd_save = fgd(1:length(fgd) - 1);
    gd_save = gd(1:length(gd)-1);
    % Ensure all the saved arrays have the same length by truncating
    % them appropiately. Hopefully we should have enough samples.    
    fid = fopen(save.file_name,'w');
    fprintf(fid,'%s %s %s %s %s %s %s\n','frequency','tf','phase','freq_gd','gd','freq_gvd','gvd');
    if strcmp(units,'eng')
        fprintf(fid,'%s %s %s %s %s %s %s\n','(GHz)','(dB)','(rad)','(GHz)','(ps)','(GHz)','(ps/nm)');
        Res = [freq_save/1.0e9;10*log10(abs(tf_save).^2);phase_save;freq_gd_save/1.0e9;gd_save;fgvd/1.0e9;gvd];
    elseif strcmp(units,'si')
        fprintf(fid,'%s %s %s %s %s %s %s\n','(Hz)','(dB)','(rad)','(Hz)','(s)','(Hz)','(s/m)');
        Res=[freq_save;10*log10(abs(tf_save).^2);phase_save;freq_gd_save;gd_save;fgvd;gvd];
    end
    fprintf(fid,'%e %e %e %e %e %e %e\n',Res);
    fclose(fid);   
end
% Optionally save the full transfer function:
% frequency magnitude phase group-delay dispersion
% in a text file

end