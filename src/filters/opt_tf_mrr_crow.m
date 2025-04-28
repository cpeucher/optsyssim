function M = opt_tf_mrr_crow(freq,params)
% Coupled microring resonator (CROW) chain matrix 
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements chain matrix computations of coupled
% micro-ring resonators of the "coupled-resonator optical waveguides"
% (CROW) type.
% The coupled ring resonators consist of:
% * 1 straight in-through waveguide
% * nrings ring resonators
% * Possibly one straight drop-add waveguide.
% For instance, for nrings = 3
%
%      drop <-  ---
%                o
%                o
%                o
%        in ->  ---  -> through
%
% All rings can have different properties (radii, effective refractive
% indices, coupling coefficients, loss, eventual phase shifts).
% The following frequency dependences are accounted for:
% * waveguide loss (in the ring)
% * effective indices (of the rings)
% * coupling coefficients (between the rings as well as to/from straight
% waveguides.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_crow.nrings = 3;
% params_crow.power_coupling = [0.45 0.5 0.3 0.45]';
% params_crow.radius = [10e-6 20e-6 10e-6];
% params_crow.phase_shift = [0 0 0];
% params_crow.neff = [2.54 2.54 2.54]';
% params_crow.loss_log = [1e-2 1e-2 1e-2]';
% M = opt_tf_mrr_crow(frequency_array,params_crow);
% tf_through = squeeze(-M(1,1,:)./M(1,2,:)).';
% tf_drop = squeeze(M(2,1,:) - M(1,1,:).*M(2,2,:)./M(1,2,:)).';
% tf_all_pass = squeeze((M(1,1,:) - M(2,1,:))./(M(2,2,:) - M(1,2,:))).';
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% freq              relative frequencies at which the transfer functions
%                       will be evaluated [real vector]
%
% params            parameters of the coupled MRR [structure]
%
%                       params.nrings
%                           number of microrings [integer scalar]
%
%                       params.power_coupling
%                           matrix containing the power coupling 
%                           coefficients (cross power coupling) versus 
%                           frequency of the params.nrings + 1 couplers. 
%
%                           If size(params.power_coupling) = [1 1],
%                               it is assumed that all params.nrings + 1 
%                               couplers have the same power coupling
%                               coefficient, which is independent of 
%                               frequency.
%
%                           If size(params.power_coupling) = [1 length(freq)],
%                               it is assumed that all params.nrings + 1 
%                               couplers have the same power coupling 
%                               coeffient, which is frequency dependent.
%
%                           If size(params.power_coupling) = [params.nrings+1 1],
%                               it is assumed that the params.nrings + 1 
%                               coupling coefficients, which may take 
%                               different values, are not
%                               frequency-dependent.
%
%                           If size(params.power_coupling) = [params.nrings+1 length(freq)],
%                               the lines of params.power_coupling are the
%                               frequency-dependent coupling coefficients 
%                               of the params.nrings + 1 couplers.
%
%                           Note that, in order to model an all pass structure,
%                           one should have
%                           params.power_coupling(params.nrings + 1,:) = 0;
%
%                       params.neff
%                           matrix containing the effective indices
%                           versus frequency of the waveguides of the 
%                           params.nrings rings. 
%
%                           If size(params.neff) = [1 1],
%                               it is assumed that all params.nrings rings
%                               have the same effective indices, which
%                               is independent of frequency (hence is equal
%                               to the group index).
%
%                           If size(params.neff) = [1 length(freq)],
%                               it is assumed that all params.nrings rings
%                               have the same effective indices, which
%                               are frequency dependent.
%
%                           If size(params.neff) = [params.nrings 1],
%                               it is assumed that the effective indices of 
%                               the params.nrings, which may take different 
%                               values, are not frequency-dependent.
%
%                           If size(params.neff) = [params.nrings length(freq)],
%                               the lines of params.neff are the frequency-
%                               dependent effective indices of the 
%                               waveguides of the params.nrings rings.
%
%                       params.loss_log
%                           matrix containing the loss of the waveguides of 
%                           each ring, in dB/m.
%
%                           If size(params.loss_log) = [1 1],
%                               it is assumed that all params.nrings rings
%                               have the same loss, which is independent of 
%                               frequency.
%
%                           If size(params.loss_log) = [1 length(freq)],
%                               it is assumed that all params.nrings rings
%                               have the loss, which is frequency
%                               dependent.
%
%                           If size(params.loss_log) = [params.nrings 1]
%                               it is assumed that the loss of the
%                               params.nrings, which may take different 
%                               values, is not frequency-dependent.
%
%                           If size(params.loss_log) = [params.nrings length(freq)],
%                               the lines of params.loss_log are the
%                               frequency-dependent loss values of
%                               the waveguides of the params.nrings rings.
%
%                       params.radius
%                           vector containing the radii of the microrings, 
%                           in m.
%
%                           The input variable is a vector of length
%                           params.nrings.
%
%                           If a single value is provided, it is assumed
%                           it is the same for all the rings.
%
%                       params.phase_shift
%                           vector containing eventual phase shifts
%                           applied to the rings, in rad.
%
%                           The input variable is a vector of length
%                           params.nrings. 
%
%                           If a single value is provided, it is assumed
%                           it is the same for all the rings.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% M                 "chain matrix" of the structure
%
%                       For an add-drop structure, the through and drop 
%                           transfer functions are:
%                           tf_through = -M{1,1}./M{1,2};
%                           tf_drop = M{2,1}-M{1,1}.*M{2,2}./M{1,2};
%                           respectively.
%
%                       For an all-pass structure, the through transfer
%                           function is: 
%                           tf_through = (M{1,1}-M{2,1})/(M{2,2}-M{1,2});
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
% Check params.power_coupling input
% -------------------------------------------------------------------------
dim_power_coupling = size(params.power_coupling);
% Dimension of the power_coupling matrix
% dim_power_coupling(1) is the number of lines, corresponding to the number 
% of power_coupling versus frequency vectors.
% dim_power_coupling(2) is the number of columns, that should match the 
% length of the freq vector.

if dim_power_coupling(1)*dim_power_coupling(2) == 1
    % First we test if params.power_coupling is a single value. If this is
    % the case, it means that all couplers have the same power_coupling, 
    % which is furthermore frequency independent.
    % We create a matrix of dimension NRings+1 x length(Freq) that stores
    % the same value of power_coupling for all couplers and all frequencies.
    params.power_coupling = params.power_coupling*ones(params.nrings + 1,length(freq));
elseif dim_power_coupling(2) == 1
    % Values of power_coupling are given only for one frequency. Therefore 
    % the power coupling coefficient does not depend on frequency.
    if dim_power_coupling(1) ~= params.nrings + 1
        % The number of values of power_coupling does not match the number 
        % of couplers.
        % We have a problem.
        error('opt_tf_mrr_crow: the number of lines in the power coupling matrix does not match the number of couplers.');    
    else
        % The number of values of power_coupling matches the number of
        % couplers.
        params.power_coupling = repmat(params.power_coupling,[1 length(freq)]);
        % We repeat the power coupling coefficients for all frequencies.
    end
elseif dim_power_coupling(1) == 1
    % Only one power_coupling vs frequency array is provided. 
    % It is therefore assumed that allcouplers  have the same values of 
    % power_coupling vs frequency.
    if dim_power_coupling(2) ~= length(freq)
        % The number of columns in the power_coupling array does not match 
        % the number of frequencies.
        % We have a problem. 
        error('opt_tf_mrr_crow: the number of columns in the power coupling matrix does not match the length of the frequency vector.');
    else
        % The number of columns in the power_coupling array matches the 
        % length of the frequency array.
        params.power_coupling = repmat(params.power_coupling,[params.nrings+1 1]);
        % We repeat power_coupling for all couplers.
    end
elseif dim_power_coupling(1) ~= params.nrings+1 || dim_power_coupling(2)~=length(freq)
    % The number of elements in the power_coupling matrix seems to be
    % unrelated to the number of rings and/or length of frequency array.
    % We have a problem.
    error('opt_tf_mrr_crow: the size of the power coupling matrix does not match the number of rings or length of frequency vector');
end 

%--------------------------------------------------------------------------
% Check params.neff input
%--------------------------------------------------------------------------
dim_neff = size(params.neff);
% Dimension of the neff matrix
% dim_neff(1) is the number of lines, corresponding to the number of neff
% versus frequency vectors.
% dim_neff(2) is the number of columns, that should match the length of the
% Freq vector.

if dim_neff(1)*dim_neff(2) == 1
    % First we test if params.neff is a single value. If this is the case,
    % it means that all rings have the same neff, which is furthermore
    % frequency independent (i.e. dispersion is not taken into account).
    % We create a matrix of dimension NRings x length(Freq) that stores the
    % same value of neff for all rings and all frequencies.
    params.neff = params.neff*ones(params.nrings,length(freq));
elseif dim_neff(2) == 1
    % Values of neff are given only for one frequency. Therefore it is 
    % assumed the waveguides are not dispersive.
    if dim_neff(1) ~= params.nrings
        % The number of values of neff does not match the number of rings.
        % We have a problem.
        error('opt_tf_mrr_crow: the number of lines in the effective index matrix does not match the number of rings.');    
    else
        % The number of values of neff matches the number of rings.
        params.neff = repmat(params.neff,[1 length(freq)]);
        % We repeat the index for all frequencies.
    end
elseif dim_neff(1) == 1
    % Only one neff array is provided. It is therefore assumed that all
    % rings have the same value of neff vs frequency.
    if dim_neff(2) ~= length(freq)
        % The number of columns in the neff array does not match the number 
        % of frequencies.
        % We have a problem. 
        error('opt_tf_mrr_crow: the number of columns in the effective index matrix does not match the length of the frequency vector.');
    else
        % The number of columns in the neff array matches the number 
        % of frequencies.
        params.neff = repmat(params.neff,[params.nrings 1]);
        % We repeat neff for all rings.
    end
elseif dim_neff(1)~=params.nrings || dim_neff(2)~=length(freq)
    % The number of elements in the neff matrix seems to be unrelated to
    % the number of rings and/or length of frequency array.
    % We have a problem.
    error('opt_tf_mrr_crow: the size of the effective index matrix does not match the number of rings or length of frequency vector');
end    



%--------------------------------------------------------------------------
% Check params.loss_log input
%--------------------------------------------------------------------------
dim_loss_log = size(params.loss_log);
% Dimension of the loss_log matrix
% dim_loss_log(1) is the number of lines, corresponding to the number of 
% loss_log versus frequency vectors.
% dim_loss_log(2) is the number of columns, that should match the length of 
% the Freq vector.

if dim_loss_log(1)*dim_loss_log(2) == 1
    % First we test if params.loss_log is a single value. If this is the 
    % case, it means that all rings have the same loss_log, which is 
    % furthermore frequency independent.
    % We create a matrix of dimension NRings x length(Freq) that stores the
    % same value of loss_log for all rings and all frequencies.
    params.loss_log = params.loss_log*ones(params.nrings,length(freq));
elseif dim_loss_log(2) == 1
    % Values of loss_log are given only for one frequency. Therefore it is 
    % assumed waveguide loss is not frequency dependent.
    if dim_loss_log(1) ~= params.nrings
        % The number of values of loss_log does not match the number of 
        % rings.
        % We have a problem.
        error('opt_tf_mrr_crow: the number of lines in the loss matrix does not match the number of rings.');    
    else
        % The number of values of loss_log matches the number of rings.
        params.loss_log = repmat(params.loss_log,[1 length(freq)]);
        % We repeat the index for all frequencies.
    end
elseif dim_loss_log(1) == 1
    % Only one loss_log array is provided. It is therefore assumed that all
    % rings have the same value of loss_log vs frequency.
    if dim_loss_log(2) ~= length(freq)
        % The number of columns in the loss_log array does not match the 
        % number of frequencies.
        % We have a problem. 
        error('opt_tf_mrr_crow: the number of columns in the loss matrix does not match the length of the frequency vector.');
    else
        % The number of columns in the loss_log array matches the number 
        % of frequencies.
        params.loss_log = repmat(params.loss_log,[params.nrings 1]);
        % We repeat loss_log for all rings.
    end
elseif dim_loss_log(1)~=params.nrings || dim_loss_log(2)~=length(freq)
    % The number of elements in the loss_log matrix seems to be unrelated to
    % the number of rings and/or length of frequency array.
    % We have a problem.
    error('opt_tf_mrr_crow: the size of the loss matrix does not match the number of rings or length of frequency vector');
end 

%--------------------------------------------------------------------------
% Check params.radius input
%--------------------------------------------------------------------------
if length(params.radius) == 1
    % No, the radii are not frequency dependent...
    % Only one single value of radius is provided. Then it is the
    % same for all rings.
    params.radius = params.radius*ones(1,params.nrings);
elseif length(params.radius) ~= params.nrings
    error('opt_tf_mrr_crow: the dimension of the ring radii array does not match the number of rings.');    
end    

%--------------------------------------------------------------------------
% Check params.phase_shift input
%--------------------------------------------------------------------------
if length(params.phase_shift) == 1
    % So far the option to have frequency dependent phase shifts is not
    % considered.
    % Only one single value of phase shift is provided. Then it is the
    % same for all rings.
    params.phase_shift = params.phase_shift*ones(1,params.nrings);
elseif length(params.phase_shift) ~= params.nrings
    error('opt_tf_mrr_crow: the dimension of the ring phase shift array does not match the number of rings.');    
end 


%--------------------------------------------------------------------------
% The inputs have been checked. We can move on.
%--------------------------------------------------------------------------

k = -1i*sqrt(params.power_coupling);
t = sqrt(1 - params.power_coupling);
% Cross and through field coupling coefficient of the couplers

alpha = params.loss_log*log(10)/10;
% Convert loss from dB/m to alpha coefficient for all ring waveguides

M(1,1,:) = -t(1,:)./k(1,:);
M(1,2,:) = 1./k(1,:);
M(2,1,:) = -1./k(1,:);
M(2,2,:) = conj(t(1,:))./k(1,:);
% Input coupler matrix


for iring = 1:params.nrings  
        
    betac = 2*pi*params.neff(iring,:).*freq/CONSTANT.c - 1i*alpha(iring,:)/2;
    % Complex propagation constant in the ring, including loss
    
    Q(1,1,:) = zeros(1,length(freq));
    Q(1,2,:) = exp(-1i*betac*pi*params.radius(iring)).*exp(-1i*params.phase_shift(iring));
    Q(2,1,:) = exp(1i*betac*pi*params.radius(iring));
    Q(2,2,:) = zeros(1,length(freq));
    % Q matrix of the ring
    
    M = prod_mm(Q,M);
    % Product of matrix of ring by matrix of coupler   
    
    if isequal(k(iring + 1,:),zeros(1,length(freq)))        
        % Check that the next coupler does not have a cross coupling 
        % coefficient equal to 0
        % If this is the case the coupled-ring resonator is an all-pass
        % structure and we need to exit the loop.
        break;
    else
        % Else we carry on and move to the next coupler
        
        P(1,1,:) = -t(iring + 1,:)./k(iring + 1,:);
        P(1,2,:) = 1./k(iring + 1,:);
        P(2,1,:) = -1./k(iring + 1,:);
        P(2,2,:) = conj(t(iring + 1,:))./k(iring + 1,:);
        % P matrice of the iring + 1 coupler
        
        M = prod_mm(P,M);
        % New transfer matrix after irings
        
    end
    
end
% End of loop over rings

end