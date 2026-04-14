function sig = set_osnr(sig,osnr_db,noise_bw,noise_pol)
% Add ASE noise corresponding to a given OSNR to a noise-free optical signal
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function is to be used to set the optical signal-to-noise ratio 
% (OSNR) of an initially noise free signal to a specified value.
% Observe that the input signal should be noise free for this function to
% make sense.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% osnr_target_db = 20;
% sig = set_osnr(sig,osnr_target_db,12.5e9,2);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input noise-free optical signal
%                       [optical signal structure]
%
% osnr_db           value of the target OSNR, in dB [real scalar]
%
%                       This value is specified in the chosen noise_bw.
%
% noise_bw          noise bandwith used to specify the OSNR, in Hz [real
%                       scalar]
%
%                       This bandwidth is typically chosen to be equal to 
%                       12.5 GHz (corresponding to 0.1 nm at 1550 nm)
%
% noise_pol         number of noise polarisations [string]
%
%                       noise_pol = '1pol'
%                           noise is only added to the -x polarisation
%                           The input signal should be polarised along -x 
%                           for this option to make sense.  
%
%                       noise_pol = '2pol'
%                           noise is added to both -x and -y polarisations
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               noisy output optical signal [optical signal structure]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% dt                time samples separation, in s [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global dt

nsamples = length(sig.x);
% Number of samples in the input signal

pav = char_opt_average_power(sig);
% Calculate average power of the signal

n0 = pav/(noise_pol*noise_bw*10^(osnr_db/10));
% Calculate the noise spectral density (per polarisation) for a desired 
% OSNR value  

ase_variance = n0/dt/2;
% ASE variance per quadrature and per polarisation
% Total noise power per polarisation is: n0*sample_rate = n0/dt

if noise_pol == 1  
    sig.x = sig.x + sqrt(ase_variance)*randn(1,nsamples) + 1i*sqrt(ase_variance)*randn(1,nsamples);
    % Noise is added to the -x polarisation
    sig.y = sig.y;
    % No noise is added to the -y polarisation. The signal should be
    % polarised along -x anyway for this calculation to make sense.   
    
elseif noise_pol == 2  
    sig.x = sig.x + sqrt(ase_variance)*randn(1,nsamples) + 1i*sqrt(ase_variance)*randn(1,nsamples);
    sig.y = sig.y + sqrt(ase_variance)*randn(1,nsamples) + 1i*sqrt(ase_variance)*randn(1,nsamples);    

else
    error('set_osnr: number of noise degrees of freedom (polarisation) should either be 1 or 2');

end