function [sig,delays] = opt_otdm_mux(sig,symbol_rate,multiplexing_factor,prbs_order,mode)
% OTDM multiplexer
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a serial PRBS preserving optical time
% multiplexer.
% Equation for the delays in the serial case is obtained from:
% Pengfei Xing, Msc thesis, Research Center COM, Technical University of
% Denmark, Oct. 2002.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% otdm_mux_multiplication_factor = 4;
% otdm_mux_mode = 'serial';% 'serial_random';'parallel';
% [sig,otdm_mux_delays] = opt_otdm_multiplexer(sig,symbol_rate,...
%       otdm_mux_multiplication_factor,params_prbs.order,otdm_mux_mode); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig                   input optical signal at the base rate
%                           [optical signal structure]
%
% symbol_rate           symbol rate of the input signal [real scalar]
%
% multiplexing_factor   bit rate multiplication factor
%                           If multiplexing_factor = 2, the bit rate 
%                               will be doubled. 
%                            If multiplexing_factor = 4 it will be 
%                               quadrupled, etc.
%
% prbs_order            order of the PRBS used to modulate the input 
%                           signal at the base rate [integer scalar]
%                           The PRBS will be preserved.
%
% mode                  structure of the multiplexer [string]
%                             
%                           mode = 'serial'
%                               PRBS preserving serial multiplexer
%                               structure
%
%                           mode = 'serial_random'
%                               serial multiplexer with random delays
%
%                           mode = 'parallel': parallel multiplexer
%                               structure.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% sig               output time multiplexed optical signal
%                       [optical signal structure]
%
% delays            values of the delays in the optical multiplexer
%                       [real vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

pshift = 0;
% PRBS shift. Set to zero. For compatibility with possible more general
% model.

if multiplexing_factor == 1
    % Do nothing if the bit rate multiplication factor is equal to 1
    delays = [];
    % Create an empty delays vector.
    
else    
    
    if strcmp(mode,'serial') || strcmp(mode,'serial_random')
        % Serial multiplexer.
        
        nstages = log2(multiplexing_factor);
        % Number of multiplexer stages.
        
        delays = zeros(1,nstages);
        % Initialise array for storing the values of delays.
        
        [sig_up,sig] = opt_coupler_2x2(opt_nosig,sig,'lin',0.5);
        % Input coupler.
        
        for kstage = 1:nstages
            
            if strcmp(mode,'serial')
                % Serial PRBS-preserving multiplexer.
                delay_stage = (pshift*(2^prbs_order - 1) + (2^prbs_order-1)/2^kstage)/symbol_rate;
                
            elseif strcmp(mode,'serial_random')
                % Serial multiplexer with random delays.
                delay_stage = 1/(2^kstage*symbol_rate) + randi([1,100])/(2^(kstage - 1)*symbol_rate);                
            end
            % End of delay calculation.
            
            [sig_up,delay_stage_actual] = opt_delay(sig_up,'off',delay_stage);
            % Implement delay in upper arm of the stage with no phase shift
            % in the delay.
            
            [sig_up,sig] = opt_coupler_2x2(sig_up,sig,'lin',0.5);
            % Coupler after multiplexing stage.
            
            delays(kstage) = delay_stage;
            % Add delay to the delays array.
        end        
    
        
    elseif strcmp(mode,'parallel')
        
        error('opt_otdm_multiplexer: parallel multiplexer structure not implemented yet.');        
        
    else
        error('opt_otdm_multiplexer: multiplexer structure not implemented.');
    end
    
end

end
% -------------------------------------------------------------------------
% End of function
% -------------------------------------------------------------------------