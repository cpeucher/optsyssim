function sig = opt_nlse_scalar_basic(sig,params,numparams)
% Basic scalar nonlinear Schroedinger equation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function solves a basic scalar nonlinear Schroedinger equation. 
% Fibre dispersion with any order of beta, nonlinear Kerr effect, and fibre 
% loss are considered. High order effects such as self-steepening and 
% stimulated Raman scattering are ignored. The nonlinear coefficient is
% assumed independent of frequency.
% The split-step calculation uses an adaptive step size determination based
% on a maximum allowed nonlinear phase rotation. A symmetric split-step
% scheme is implemented. No thrill.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_fibre.nonlinear_coefficient = 1.3e-3;% nonlinear coefficient, in 1/W/m
% params_fibre.dispersion = 16;% dispersion, in ps/nm/km
% params_fibre.dispersion_slope = 0.058;% dispersion slope, in ps/nm2/km
% params_fibre.dispersion_curvature = 0;% dispersion curvature, in ps/nm3/km
% params_fibre.dispersion_spec_frequency = reference_frequency;
% params_fibre.loss = 0.2;% loss, in dB/km
% params_fibre.length = 10e3;% fibre length, in m
% numparams_fibre.max_step_size = 1;% maximum step size, in m
% numparams_fibre.max_phase_shift = 1e-3;% maximum nonlinear phase shift, in radians
% params_fibre.beta_coefficients = conv_disp_d_beta([params_fibre.dispersion params_fibre.dispersion_slope params_fibre.dispersion_curvature],'to_beta','eng','si',params_fibre.dispersion_spec_frequency);
% params_fibre.loss_alpha = conv_loss_lin_log(params_fibre.loss);
% sig = opt_nlse_scalar_basic(sig,params_fibre,numparams_fibre);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input optical signal [optical signal structure]
%
%                       The signal should be polarised along the 
%                       -x direction.
%
% params            fibre parameters [structure]
%
%                       params.length       
%                           fibre length, in m [real scalar]
%
%                       params.loss_alpha   
%                           fibre attenuation, in 1/m [real vector]
%
%                           If params.loss_alpha consists of a single 
%                           element, it is assumed that the loss is 
%                           frequency independent. 
%                           Otherwise params.loss_alpha should have 
%                           length(sig.x) elements and represents the
%                           frequency dependent attenuation.
%
%                       params.beta_coefficients  
%                           fibre dispersion coefficients in SI units
%                           [real vector]
%
%                           [beta_1, beta_2, beta_3 ... beta_n] are terms 
%                           of the Taylor expansion of the propagation 
%                           constant as obtained by fitting a dispersion 
%                           versus wavelength curve or from direct
%                           conversion from [D S C] using the appropriate 
%                           functions.
%                           beta_1 is only present in the array for 
%                           consistency and its value is not relevant. 
%
%                       params_fibre.dispersion_spec_frequency
%                           frequency at which the dispersion parameters
%                           are provided, in Hz [real scalar]
%
%                       params_fibre.nonlinear_coefficient
%                           fibre nonlinear coefficient in 1/W/m 
%                           [real scalar]
%                   
% numparams         solver parameters [structure]
% 
%                       numparams.max_step_size
%                           maximum split-step size, in m [real scalar]
%   
%                       numparams.max_phase_shift
%                           maximum nonlinear phase shift, in radians 
%                           [real scalar]  
%
%                       The step length is determined according to the most
%                       constraining parameter.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output optical signal [optical signal structure]
% 
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% frequency_array       relative frequency samples, in Hz [real vector]
%
% reference_frequency   reference frequency, in Hz [real scalar]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global frequency_array 
global reference_frequency


nsamples = length(sig.x);
% Number of samples in the signal

% -------------------------------------------------------------------------
% Checking the inputs
% -------------------------------------------------------------------------
if sig.y ~= 0
    error('opt_nlse_scalar_basic (input): input signal should be polarised along -x.');
end
% Check that the signal is polarised along the -x direction

if length(params.nonlinear_coefficient) ~= 1
    error('opt_nlse_scalar_basic (input): only one value of nonlinear coefficient is expected.');
end
% Check that only one value is provided for the nonlinear coefficient

nalpha = length(params.loss_alpha);
if nalpha == 1    
    
        disp('opt_nlse_scalar_basic (input): loss is frequency independent.');
   
    elseif nalpha == nsamples
   
        disp('opt_nlse_scalar_basic (input): loss is frequency dependent.');
    
else
    error('opt_nlse_scalar_basic (input): wrong number of elements in the input array params.loss_alpha. Should be 1 or nsamples');
end
% Check that the input array params.loss_alpha contains one element (in case 
% the fibre loss is frequency independent, or Nsamples elements (in case 
% the fibre loss is frequency dependent).


% -------------------------------------------------------------------------
% Calculation of linear operator depending on the nature of the input loss
% and dispersion parameters
% -------------------------------------------------------------------------
angular_frequency_array = 2*pi*frequency_array;
% Calculate angular frequencies array omega=2*pi*f for convenience

wref  = 2*pi*(params.dispersion_spec_frequency - reference_frequency);
% Convert reference frequency for the beta_i's to angular frequency
% Observe at this point that:
% angular_frequency_array is a relative angular frequency array with '
% respect to the reference frequency of the simulation
% wref is now also a relative angular frequency with respect to the
% reference frequency of the simulation.

lin_op = zeros(1,nsamples);
% Preallocate linear operator
for idisp = 2:length(params.beta_coefficients)
    lin_op = lin_op + params.beta_coefficients(idisp)/factorial(idisp)*(angular_frequency_array - wref).^idisp;
end
% Dispersion term of the linear operator
lin_op = -1i*lin_op - params.loss_alpha/2;
% Total linear operator, including loss


ssig = sig.x;
% Scalar input signal is now only the -x polarisation component.

% -------------------------------------------------------------------------
% Isolate simple cases that do not require split-stepping
% -------------------------------------------------------------------------
if params.nonlinear_coefficient == 0
    % If no nonlinearity is present we just apply the linear operator.
    
    disp('opt_nlse_scalar_basic (run): the medium is actually linear.');
    
    ssig = fft(ssig);
    % Now ssig is in the frequency domain
    ssig = ssig.*fftshift(exp(lin_op*params.length));
    % Apply the linear operator
    ssig = ifft(ssig);
    % Back to the time domain
    % End of calculation in the case params.nonlinear_coefficient = 0
    
elseif params.beta_coefficients == 0 & length(params.loss_alpha) == 1
    % If no dispersion is present and the loss is frequency independent, we 
    % just apply the nonlinear operator.
    
        disp('opt_nlse_scalar_basic (run): the medium is actually not dispersive.');
        disp('opt_nlse_scalar_basic (run): loss is frequency independent.');
    
 
    if params.loss_alpha == 0
        % params.loss_alpha =0 case is distinguished for the calculation of the
        % effective length
        
        disp('opt_nlse_scalar_basic (run): the medium is actually not lossy.');
        
        leff = params.length;
    else
        leff = (1-exp(-params.loss_alpha*params.length))/params.loss_alpha;
        % Fibre effective length
    end
    % End of calculation of effective length
    
    ssig = ssig.*exp(-1i*params.nonlinear_coefficient*leff*abs(ssig).^2)*exp(-params.loss_alpha*params.length/2);
    % Signal at the output of the nonlinear fibre
    
else
% -------------------------------------------------------------------------
% Else we need to split step.
% -------------------------------------------------------------------------
   
    disp('opt_nlse_scalar_basic (run): split-step calculation running.');
    
    z = 0;
    % Initialise current position in the fibre
    
    nstep = 0;
    % Initialise current step

    ssig = fft(ssig)/nsamples;
    % The signal is now in the frequency domain and in fft order

    while z < params.length
        
        dz = min(numparams.max_step_size,numparams.max_phase_shift/(params.nonlinear_coefficient*max(abs(ssig).^2)));
        % First we determine the next step size, which is the minimum of
        % the 2 values determined by the maximum step size and the maximum
        % nonlinear phase shift.     

        dz = min(dz,params.length - z);
        % The step size needs to be adjusted when we reach the end of the
        % fibre.
        
        ssig = ifft(ssig.*fftshift(exp(lin_op*dz/2)))*nsamples;
        % Apply linear operator over first half of the step size and
        % convert back to the time domain

        ssig = ssig.*exp(-1i*params.nonlinear_coefficient*dz*abs(ssig).^2);
        % Apply nonlinear operator over the step size and remain in the
        % time domain   

        ssig = fft(ssig).*fftshift(exp(lin_op*dz/2))/nsamples;
        % Apply linear operator over second half of the step size
        % the signal is now in the frequency domain and in fft
        % order.      

        z = z + dz;
        % Increment position in fibre
        
        nstep = nstep + 1;
        % Increment number of completed steps
        
    end
    % End of while loop over fibre length.
    % The signal is in the frequency domain and in fft order when it exits
    % the loop.
    
    ssig = ifft(ssig)*nsamples;
    % Back to the time domain.
    
end
    
sig.x = ssig;
% Back to the 2 polarisation signal structure.

end