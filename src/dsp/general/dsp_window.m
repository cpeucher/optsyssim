function w = dsp_window(params)
% Window functions for DSP applications
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates some standard window functions for digital
% signal processing applications.
% The following functions are presently implemented:
% Kaiser        beta = pi*alpha                 symmetric (default)  
% Hann          no parameter                    symmetric, periodic
% Hamming       no parameter                    symmetric, periodic
% Bartlett      no parameter                    symmetric (default)
% Rectangular   no parameter
% These functions are typically implemented in the Signal Processing
% Toolbox of Matlab and the definitions used here should be fully compliant
% with those of this toolbox. These functions are to be used as substitutes
% if one has no access to this toolbox, which is my present situation.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_window.type = 'kaiser';%'hamming';'hann';'bartlett';'rectangular';
% params_window.length = 15;
% params_window.symmetry = 'periodic';%'symmetric';
% params_window.beta = 2.4*pi;% For Kaiser window.
% w = dsp_window(params_window); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            parameters of window function [structure]
%
%                       params.length
%                           number of points of the window [integer scalar]
%
%                       params.type
%                           window type [string]
%
%                           Supported types are: 
%                               'kaiser':           Kaiser window
%                               'hann':             Hann window
%                               'hamming':          Hamming window
%                               'bartlett':         Bartlett window
%                               'rectangular':      rectangular window
%
%                       params.symmetry
%                           symmetry type of the window for some windows
%                               [string]
%
%                               'symmetric':        symmetric
%                               'periodic':         periodic (for FFT
%                                                       purpose)
%
%                           These definitions are compatible with those 
%                           used in the 'Signal processing toolbox' of 
%                           Matlab.
%
%                       params.beta
%                           beta parameter (for Kaiser window)
%                               [real scalar]
%                               beta=pi*alpha
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% w                 window function [real vector]
%
%                       The value of the window is within [0,1]. 
%
%                       Note that depending on odd/even number of points
%                       the value of 1 may not be reached. 
%                       The value of 0 is clearly not reached for all types
%                       of windows.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

switch  params.type
    % Switch over type of window
    
    
    %----------------------------------------------------------------------
    % Kaiser window
    %----------------------------------------------------------------------
    case 'kaiser'
        % Kaiser window
        
        nn = [0:1:params.length-1];        
        w = besseli(0,params.beta*sqrt(1 - (2*nn/(params.length - 1) - 1).^2))/besseli(0,params.beta);  
        % Calculate window of length params.length
        
    %----------------------------------------------------------------------
    % Hanning window
    %----------------------------------------------------------------------    
    case 'hann'
        % Hann window        
        
        switch params.symmetry
            
            case 'periodic'
                % Suitable for DFT/FFT calculations
                
                nn = [0:1:params.length];                
                w = 0.5 - 0.5*cos(2*pi*nn/(params.length));
                % Calculate window of length params.length + 1               
                w = w(1:params.length);
                % Truncate the calculated window to the first params.length
                % points               
                
            case 'symmetric'
                
                nn = [0:1:params.length-1];                
                w = 0.5 - 0.5*cos(2*pi*nn/(params.length-1));
                % Calculate window of length params.length               
                
            otherwise
                
               error('dsp_window: symmetry type not defined for Hann window.');                
                
        end
        % End of options for Hann window       
        
    %----------------------------------------------------------------------
    % Hamming window
    %----------------------------------------------------------------------    
    case 'hamming'
        % Hamming window
        
        switch params.symmetry
            
            case 'periodic'
                % Suitable for DFT/FFT calculations
                
                nn = [0:1:params.length];                
                w = 0.54 - 0.46*cos(2*pi*nn/(params.length));
                % Calculate window of length params.length + 1               
                w = w(1:params.length);
                % Truncate the calculated window to the first params.length
                % points               
                
            case 'symmetric'
                
                nn = [0:1:params.length-1];                
                w = 0.54 - 0.46*cos(2*pi*nn/(params.length-1));
                % Calculate window of length params.length               
                
            otherwise
                
               error('dsp_window: symmetry type not defined for Hamming window.');                
                
        end
        % End of options for Hamming window               
        
    %----------------------------------------------------------------------
    % Bartlett window
    %----------------------------------------------------------------------    
    case 'bartlett'
        % Bartlett window
        
        
    check_even = mod(params.length,2);
    % To check if the window length is odd or even. 
    % This is known not to be robust for large numbers, but this should be
    % fine for the cases of interest here. Otherwise look for some of the
    % isodd or iseven functions proposed on Matlab file exchange.  
    
    switch check_even 
        
        case 0
            % The window length is even
            nn1 = [1:1:params.length/2];
            nn2 = [params.length/2 + 1:1:params.length];
            ww1 = 2*(nn1 - 1)/(params.length - 1);
            ww2 = 2*(params.length - nn2)/(params.length - 1);            
            
        case 1
            % The window length is odd
            nn1 = [0:1:(params.length - 1)/2];
            ww1 = 2*nn1/((params.length - 1));
            nn2 = [(params.length-1)/2 + 1:1:params.length - 1];
            ww2 = 2 - 2*nn2/((params.length - 1));           
    end
    % End of switch over even/odd window size
    
    w = [ww1 ww2];
    % Assemble the window from its domain definitions 
    
    
    %----------------------------------------------------------------------
    % None
    %----------------------------------------------------------------------
    case 'rectangular'
        % Rectangular window
        w = ones(1,params.length);       
        
        
    otherwise
        
        error('dsp_window: window type not defined.');
        
        
end
% End of switch over window type.


end