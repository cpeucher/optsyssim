function tf = opt_tf_fbg(freq,params,numparams)
% Transfer function of fibre Bragg gratings
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer functions of fibre Bragg grating
% filters (uniform, chirped, apodised, with or without average refractive 
% index variations) from coupled mode equations with the transfer matrix 
% formalism.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_fbg.length = 0.03;                 % Grating length, in m.
% params_fbg.centre_frequency = 193.1e12;   % Centre frequency of the grating, in Hz.
% params_fbg.n0 = 1.45;                     % Effective index of the fibre.
% params_fbg.dn = 1.0e-4;                   % Refractive index change.
% params_fbg.m = 0;                         % Control of average index change.
% params_fbg.chirp_rate = 0;                % Laser linear chirp rate, in 1/m.
% params_fbg.apodisation.type = 'uniform';     
% params_fbg.apodisation.eta = 0;
% params_fbg.apodisation.fwhm = 0;
% params_fbg.apodisation.profile = 0;
% numparams_fbg.nsections = 1000;           % Number of uniform grating sections.
% params_fbg.phase_shift = zeros(1,numparams.nsections);% Position of phase shifts.
% freq = frequency_array + reference_frequency;
% tf = opt_tf_fbg(freq,params_fbg,numparams_fbg);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              frequencies at which the transfer function will be
%                       evaluated, in Hz [real vector]
%
%                       Those are absolute frequencies.
%
% params            grating physical parameters [structure]
%
%                       params.length
%                           length of the grating, in m [real scalar]
%
%                       params.centre_frequency
%                           centre frequency of the grating passband / 
%                           stopband, in Hz [real scalar]
%
%                       params.dn
%                           amplitude of the refractive index change, no 
%                           unit, [real scalar]
%
%                       params.m
%                           parameter controlling the average refractive 
%                           index, no unit [real scalar]
%                           params.m is in the interval [0 1].
%                           If params.m = 0, there is no increase of the 
%                           average index.
%
%                       params.chirp_rate
%                           chirp rate of the grating, in 1/m [real scalar]
%                           Only linear chirp is currently supported.
%
%                       params.apodisation= 
%                           data about the apodisation profile [structure]
%
%                           params.apodisation.type
%                               type of apodisation [string]

%                               params.apodisation.type = 'uniform';
%                               params.apodisation.type = 'raised_cosine';
%                               params.apodisation.type = 'gaussian';
%                               params.apodisation.type = 'gaussian_relative';
%                               params.apodisation.type = 'quarter_cosine';
%                               params.apodisation.type = 'raised_cosine_flat_top';
%                               params.apodisation.type = 'blackman';
%                               params.apodisation.type = 'hamming';
%                               params.apodisation.type = 'custom';
%
%                           params.apodisation.fwhm 
%                               FWHM of the apodisation profile, in m
%                               [real scalar]
%
%                           params.apodisation.eta
%                               eta parameter for some apodisation profiles
%                               [real scalar]
%                           
%                           params.apodisation.profile
%                               apodisation profile in case of custom
%                               profile [real vector]
%
%                       params.phase_shift
%                           eventual phase shifts of the grating 
%                               [real vector]
%
%                           The positions of the phase shifts, expressed in 
%                           terms of section number, have to be calculated
%                           by the user.
%
%                           params.phase_shift(index)= dphi indicates 
%                           that a phase shift of dphi radians will be 
%                           present at the output of the uniform grating
%                           section indexed by index.
%
%                           In the absence of phase shifts, one would use:
%                           params.phase_shift = zeros(1,numparams.nsections);
%
% numparams         numerical parameters [structure]
%
%                       numparams.nsections
%                           number of uniform sections in which the grating
%                           is divided [integer scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                transfer functions of the grating [structure]
%
%                       tf.transmission
%                           transfer function in transmission [real vector]
%
%                       tf.reflection
%                           transfer function in reflection [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT          essential physical constants [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global CONSTANT


T{1,1} = ones(1,length(freq));
T{1,2} = zeros(1,length(freq));
T{2,1} = T{1,2};
T{2,2} = T{1,1};
% Initialise the transfer matrix. Without a grating structure, it is unit
% matrix.


lambda0 = CONSTANT.c/params.centre_frequency/2/(params.n0 + params.m*params.dn);
% Pitch at the centre (or more precisely where the apodisation function 
% takes its maximum value, which coincides with the middle of the grating
% with the classic symmetrical apodisation function
% Lambda0=CONSTANT.c/params.centre_frequency/2/params.n0; is corrected by the
% change of average refractive index that can be introduce if the parameter
% params.m is different from zero.


secparams.dz = params.length/numparams.nsections;
% Grating step size, in m

secparams.phi = 0;
% Initialise the phase of the grating for the 0th section. Taken
% arbitrarily equal to 0.
secparams.beta_b = 0;
% Initialise the betaB=pi/Lambda parameter for the 0th section. The value
% for the first section will be calculate in the following loop. Need to be
% defined here for compatibility with the phase continuity equation within
% the following loop.


for isection = 1:numparams.nsections
    % Start loop over grating sections 
    
    z = -params.length/2+(isection - 1)*secparams.dz;
    % Increment position along the grating
    
    
    secparams.neff = params.n0 + params.m*params.dn*opt_fbg_apodisation(z,params);
    % Average effective index of the present section
    
    secparams.dneff = params.dn*opt_fbg_apodisation(z,params);
    % Amplitude of the effective index variations of the present section
    
    secparams.lambda_pitch = lambda0 + params.chirp_rate*z;
    % Pitch of the present grating section, in 1/m
    
    secparams.beta_b = pi/secparams.lambda_pitch;
    % betaB= pi/Lambda of the present section
    
    tf_sec = opt_fbg_transfer_matrix(freq,secparams);
    % Compute transfer matrix of the section
    
    secparams.phi = secparams.phi + 2*secparams.beta_b*secparams.dz + params.phase_shift(isection);
    % Update the phase of the grating at the section output boundary
    % from its value at the section input boundary.
    % secparams.beta_b is pi/Lambda of the present section.
    % params.phase_shift is an eventual phase shift at the output of the
    % present section.
    
    tf_tmp = T;
    % Save the intermediate transfer matrix up to this section
    
    T{1,1} = tf_sec{1,1}.*tf_tmp{1,1} + tf_sec{1,2}.*tf_tmp{2,1};
    T{1,2} = tf_sec{1,1}.*tf_tmp{1,2} + tf_sec{1,2}.*tf_tmp{2,2};
    T{2,1} = tf_sec{2,1}.*tf_tmp{1,1} + tf_sec{2,2}.*tf_tmp{2,1};
    T{2,2} = tf_sec{2,1}.*tf_tmp{1,2} + tf_sec{2,2}.*tf_tmp{2,2};
    % Product of the tranfer matrices of the section (left) with the
    % intermediate transfer matrix (right)
    
end
% End of loop over grating sections


tf.reflection = -T{2,1}./T{2,2};
% Grating transfer function in reflection

tf.transmission = T{1,1} - (T{2,1}.*T{1,2})./T{2,2};
% Grating transfer function in transmission

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Section scattering matrix
% -------------------------------------------------------------------------
function ttf = opt_fbg_transfer_matrix(freq,secparams)

global CONSTANT

delta_beta = secparams.neff*2*pi*freq/CONSTANT.c - secparams.beta_b;
% \Delta\beta

lambda_bragg = 2*secparams.neff*secparams.lambda_pitch;
% Local Bragg wavelength
% \lambda_B= 2 n_eff \Lambda

kappa = pi*secparams.dneff/lambda_bragg;
% \kappa = \pi * \delta n_eff \lambda_B

gamma = sqrt(kappa*kappa - delta_beta.*delta_beta);
% \gamma^2 = \kappa^2 + \Delta\beta^2

exp11 = exp(-1i*secparams.beta_b*secparams.dz);
exp12 = exp11*exp(-1i*secparams.phi);
% Phase terms.

ttf{1,1} = (cosh(gamma*secparams.dz) - 1i*delta_beta./gamma.*sinh(gamma*secparams.dz)).*exp11;
ttf{1,2} = -(1i*kappa./gamma.*sinh(gamma*secparams.dz)).*exp12;
% T1 = conj(T{1,2});
% T2  =conj(T{1,1});
ttf{2,1} = (1i*kappa./gamma.*sinh(gamma*secparams.dz)).*conj(exp12);
ttf{2,2} = (cosh(gamma*secparams.dz) + 1i*delta_beta./gamma.*sinh(gamma*secparams.dz)).*conj(exp11);
% Elements of the transfer matrix.

end
% End of opt_fbg_transfer_matrix function.
% -------------------------------------------------------------------------
% End of scattering matrix
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Grating apodisation profile
% -------------------------------------------------------------------------
function alpha = opt_fbg_apodisation(z,params)
% Grating apodisation profile

switch params.apodisation.type
    
    case 'uniform'
        
        alpha = ones(1,length(z));
        
    case 'raised_cosine'
        
        alpha = 0.5*(1 + cos(pi*z/params.apodisation.fwhm));
        
    case 'gaussian'
        % Gaussian apodisation profile expressed in term of the width of
        % the Gaussian index profile
        
        alpha = exp(-(4*log(2)*z.*z)/(params.apodisation.fwhm*params.apodisation.fwhm));
        
        
    case 'gaussian_relative'
        % Alternative normalised expression for a Gaussian apodisation
        % profile
        % Included by lazyness for compatibility with the definitions 
        % of some of the other apodisation profiles expressed in terms of a
        % normalised (to the grating length) coordinate and where all other
        % relevant  physical parameters are reduce to a single eta value.
        
        alpha = exp(-params.apodisation.eta*z.*z/(params.length*params.length));
        
    case 'quarter_cosine'
        
        alpha = ones(1,length(z));
        % Initialise the appodisation profile to one
        
        alpha(abs(z) >= params.length*(0.5-params.apodisation.eta)) = cos(pi/2*(abs(z(abs(z)>= params.length*(0.5-params.apodisation.eta))) + params.length*(params.apodisation.eta-0.5))/(params.apodisation.eta*params.length));
        % if (abs(z)<(L*(0.5-eta))) 
        %   alpha=1;
        % else 
        %   z_prime=(abs(z)+L*(eta-0.5))/(eta*L);
        %   alpha=cos(pi*z_prime/2);        
        
    case 'raised_cosine_flat_top'  
        
        alpha = ones(1,length(z));
        % Initialise the appodisation profile to one
        
        alpha(abs(z) >= params.length*(0.5-params.apodisation.eta))=0.5*(1+cos(pi/2*(abs(z(abs(z)>= params.length*(0.5-params.apodisation.eta)))+params.length*(params.apodisation.eta-0.5))/(params.apodisation.eta*params.length)));
        % if (abs(z)<(L*(0.5-eta))) 
        %   alpha=1;
        % else 
        %   z_prime=(abs(z)+L*(eta-0.5))/(eta*L);
        %   alpha=0.5*(1+cos(pi*z_prime/2));   

        
    case 'blackman'
        % Blackman apodisation profile
        
        alpha = (1 + (1 + params.apodisation.eta)*cos(2*PI*z/params.length)+params.apodisation.eta*cos(4*PI*z/params.length))/(2+2*params.apodisation.eta);
        
    case 'hamming'
        % Hamming apodisation profile
        
        alpha = (1 + params.apodisation.eta*cos(2*pi*z/params.length))/(1 + params.apodisation.eta);  

        
    case 'custom'
        % Custom apodisation profile. The values of the normalised 
        % apodisation function are defined as an input vector.        
        
        alpha = params.apodisation.profile;         
   
        
    otherwise
        % Apodisation profile not defined.
        
        error('opt_fbg_apodisation: apodisation type not implemented.');
        
end


end
% -------------------------------------------------------------------------
% End of grating apodisation profile
% -------------------------------------------------------------------------
