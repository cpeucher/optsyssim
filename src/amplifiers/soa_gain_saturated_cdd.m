function hN = soa_gain_saturated_cdd(g0,psat,pin)
% Exponential gain due to carrier density depletion under steady state
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function returns the exponential gain hN of an SOA under saturation
% when the effects of spectral hole burning (SHB) and carrier heating (CH)
% are neglected.
% Only carrier density depletion (CDD) is accounted for.
% 
% It solves the implicit equation: 
% psat*(g0 - hN)/pin = exp(hN) - 1
%
% The gain is then:
% G = exp(hN)
%
% The function is used to return an initial value for the carrier density 
% depletion differential equation in the model: 
% D. Cassioli, S. Scotti, and A. Mecozzi, "A time-domain computer simulator 
% of the nonlinear response of semiconductor optical amplifiers," 
% IEEE Journal of Quantum Electronics 36, 1072 (2000) 
% DOI: 10.1109/3.863960
% In this case pin is taken as the average power of the input to the SOA.
%
% It is also used to determine the exponential gain obtained with a CW 
% input whose power is equal to the mean power of the input signal in the 
% intermediate to large signal Saleh model:
% A.A.M. Saleh, "Modeling of nonlinearity in semiconductor optical 
% amplifiers" in IEEE Global Telecommunications Conference (Dallas, TX, 
% USA, 1989), pp. 665–670.
% DOI: 10.1109/GLOCOM.1989.64051
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% g0                small-signal exponential gain [real scalar]
%
%                       The SOA small-signal gain is:
%
%                           G0 = exp(g0)
%
% psat              saturation power, in W [real scalar]
%
%                        This is linked to the input saturation power
%                        Psat,in (SOA input power at which the gain is 
%                        half of its small signal value) according to
%       
%                        Psat = (G0 - 2)/(2*log(2)) Psat,in
% 
%                        and to the output saturation power Psat,out 
%                        (SOA output power at which the gain is half of
%                        its small signal value) according to
%
%                        Psat = (G0 - 2)/(G0*log(2))*Psat,out
%
%                        (Psat, Psat,in, Psat,out in W; G0 in linear
%                           units; log(2) is natural (neperian) logarithm)
%
% pin               input power, in W [real scalar]    
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% hN                exponential gain under pin [real scalar]
% 
%                       The SOA gain, neglecting SHB and SH is then:
%
%                           G = exp(hN)
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

func = @(hn) (psat*(g0 - hn) + pin*(1 - exp(hn)));
% Define implicit function that needs to be solved to find h_N under steady
% state. 
    
hN = fzero(func,0);

end