function dml_response = calc_dml_frequency_response(ibias,params,freq)
% Calculation of small signal frequency response of a DML
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements textbook calculations of the small signal
% frequency response of a directly modulated laser.
% The definitions and notations correspond to the ones of my note:
% C. Peucheret,``Direct current modulation of semiconductor lasers,''
% Lecture note in 34130 Introduction to Optical Communication, Department
% of Photonics Engineering, Technical University of Denmark, version 2011.
% They are a mix of notations from versious texts, including:
% L. A Coldren and S. W. Corzine, ``Diode lasers and photonic integrated
% circuits,'' Wiley, New York, 1995.
% G. P. Agrawal and N. K. Dutta, ``Semiconductor lasers,'' 2nd edition, 
% Kluwer Academic Publishers, Dordrecht, 1993.
% K. Petermann, ``Laser Diode Modulation and Noise,'' Kluwer Academic 
% Publishers, Dordrecht, 1988.
% A. Yariv, ``Optical Electronics in Modern Communications,'' 5th edition, 
% Oxford University Press, New York, 1997.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_dml.tau_p = 2.6e-12;          % photon lifetime, in s
% params_dml.tau_c = 3.17e-9;          % carrier lifetime, in s
% params_dml.n_0 = 2.0e24;             % carrier density at transparency, in 1/m^3
% params_dml.sigma_g = 3.34e-20;       % gain cross section, in m^2
% params_dml.n_g = 4;                  % group effective index
% params_dml.Gamma = 0.2408;           % confinement factor
% params_dml.V = 3.6e-17;              % active volume, in m^3
% params_dml.epsilon_nl = 2.0e-23;     % gain suppression factor
% params_dml.alpha = 6;                % linewidth enhancement factor
% params_dml.beta = 1.0e-3;            % spontaneous emission factor
% params_dml.eta_0 = 0.2;              % differential quantum efficiency
% params_dml.emission_frequency = reference_frequency; % emission frequency, in Hz
% ibias = 20e-3;                       % bias current, in A
% freq = [0:0.1:20]*1e9;               % frequency range, in Hz
% dml_response = calc_dml_frequency_response(ibias,params_dml,freq);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% ibias             bias current, in A [real scalar]
%
% params            laser physical parameters [stucture]      
%
%                       params.tau_p
%                           photon lifetime, in s
%
%                       params.tau_c
%                           carrier lifetime, in s
%
%                       params.n_0
%                           carrier density at transparency, in 1/m^3
%
%                       params.sigma_g
%                           gain cross section, in m^2
%
%                       params.n_g
%                           group index
%                           
%                       params.Gamma
%                           confinement factor
%
%                       params.V
%                           active volume, in m^3
%
%                       params.epsilon_nl
%                           gain suppression factor, in m^3
%
%                       params.alpha
%                           linewidth enhancement factor
%
%                       params.beta
%                           spontaneous emission factor
%               
%                       params.eta_0
%                           differential quantum efficiency
%               
%                       params.emission_frequency
%                           emission frequency, in Hz
%
% freq              frequency values, at which the response will be 
%                       calculated in Hz [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% dml_response      frequency response and other calculated parameters
%                   	[structure]
%
%                       dml_response.Ith
%                           threshold current, in A
%
%                       dml_response.SlopeEfficiency
%                           slope efficiency, in W/A
%
%                       dml_response.Pout
%                           output power, in W
%
%                       dml_response.fR
%                           relaxation oscillations frequency, in Hz
%
%                       dml_response.tR
%                           relaxation oscillations period, in s
%
%                       dml_response.tDamp
%                           relaxation oscillations damping time constant, 
%                           in s
%
%                       dml_response.f3dB
%                           modulation 3-dB bandwidth, in Hz
%
%                       dml_response.f3dB_approx
%                           approximate modulation 3-dB bandwidth, in Hz
%
%                       dml_response.omega_0
%                           normalisation angular frequency for plotting
%                           normalised frequency responses, in rad/s
%
%                       dml_response.gamma
%                           normalised damping parameter, in 1/s
%
%                       dml_response.H
%                           small signal frequency response in linear scale
%
%                       dml_response.HdB
%                           small signal frequency responce in logarithmic 
%                           scale (10*log10(H))
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT          essential physical constants [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global CONSTANT


% -------------------------------------------------------------------------
% Calculation of laser threshold
% -------------------------------------------------------------------------
v_g = CONSTANT.c/params.n_g;
% Group velocity

G_N = (params.Gamma*v_g*params.sigma_g)/params.V;
% Ratio between the net rate of stimulated emission and the excess carrier
% number G=G_N (N - n_0)

n_0 = params.n_0*params.V;
% Carrier number calculated from the carrier density

N_th = n_0 + 1/(G_N *params.tau_p);
% Clamped number of carriers

dml_response.Ith = CONSTANT.q*N_th/params.tau_c;
% Laser threshold

fprintf(1,'\n\n\n%s','Directly modulated laser analytical calculations:');
fprintf(1,'\n%s','-------------------------------------------------');
% Prepare results display on screen
fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','Ith=',dml_response.Ith/1.0e-3,'mA','Threshold current');
% Display laser threshold on screen


% -------------------------------------------------------------------------
% Calculation of slope efficiency
% -------------------------------------------------------------------------

dml_response.SlopeEfficiency = CONSTANT.h*params.emission_frequency*params.eta_0/2/CONSTANT.q;
% Laser slope efficiency

fprintf(1,'\n%s\t\t%3.4f\t%s\t\t\t%s','DP/DI=',dml_response.SlopeEfficiency,'W/A','Slope efficiency');
% Display laser threshold on screen


% -------------------------------------------------------------------------
% Calculation of output power
% -------------------------------------------------------------------------

if ibias > dml_response.Ith
    % The bias current is above the lasing threshold.   
    dml_response.Pout = dml_response.SlopeEfficiency*(ibias - dml_response.Ith);    
else
    % The bias current is below the lasing threshold.
    dml_response.Pout = 0;    
end

fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','Ibias=',ibias/1.0e-3,'mA','Bias current');
fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','Pout=',dml_response.Pout/1.0e-3,'mW','Output power');
% Display laser output power on screen



% -------------------------------------------------------------------------
% Calculation of small signal frequency response
% -------------------------------------------------------------------------
G_b = 1/params.tau_p;
P_b = (params.tau_p/CONSTANT.q)*(ibias-dml_response.Ith);

gamma_PP = params.epsilon_nl*G_b*P_b;
gamma_PN = -G_N*P_b;
gamma_NP = G_b;
gamma_NN = G_N*P_b + 1/params.tau_c;

gamma_R = (gamma_PP + gamma_NN)/2;
Omega_R = sqrt(-gamma_PN*gamma_NP - 0.25*(gamma_NN - gamma_PP)^2);


dml_response.fR = Omega_R/(2*pi);
% Relaxation oscillations frequency
dml_response.tR = 1/dml_response.fR;
% Relaxation oscillations period

dml_response.tDamp = 1/gamma_R;
% Damping parameter 1/gamma_R

fprintf(1,'\n%s\t\t\t%3.2f\t%s\t\t\t%s','fR=',dml_response.fR/1.0e9,'GHz','Relaxation oscillations frequency');
% Display relaxation oscillations frequency on screen
fprintf(1,'\n%s\t\t\t%3.2f\t%s\t\t\t%s','TR=',dml_response.tR/1.0e-12,'ps','Relaxation oscillations period');
% Display corresponding relaxation oscillation period on screen
fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','tDamp=',dml_response.tDamp/1.0e-12,'ps','Relaxation oscillations damping time constant');
% Display damping parameter on screen


dml_response.f3dB = sqrt(Omega_R^2 - gamma_R^2 + 2*sqrt(gamma_R^4 + Omega_R^4 + Omega_R^2*gamma_R^2))/(2*pi);
% Laser 3 dB bandwidth

fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','f3dB=',dml_response.f3dB/1.0e9,'GHz','Modulation 3-dB bandwidth');
% Display 3 dB bandwidth on screen

dml_response.f3dB_approx=sqrt(3)*Omega_R/(2*pi);
% Approximate expression for the 3 dB bandwidth

fprintf(1,'\n%s\t\t%3.2f\t%s\t\t\t%s','f3dB ~',dml_response.f3dB_approx/1.0e9,'GHz','Approximate modulation 3-dB bandwidth');
% Display approximate 3 dB bandwidth on screen

dml_response.omega_0 = sqrt(Omega_R^2 + gamma_R^2);
% Normalised frequency for the modulation transfer function
dml_response.gamma = gamma_R/dml_response.omega_0;
% Normalised damping parameter

omega_m = freq*2*pi;
% Convert the input modulation frequencies array to angular frequencies

dml_response.H = 1./(1 - (omega_m/dml_response.omega_0).^2 + 2*1i*dml_response.gamma*(omega_m/dml_response.omega_0));
% Calculate small signal frequency response

dml_response.HdB = 10*log10(abs(dml_response.H));
% Convert to dB scale

end
