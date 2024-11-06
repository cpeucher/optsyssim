function out = conv_disp_d_beta(in,mode,units_in,units_out,freq)
% Conversion of dispersion between D, S, C and beta_n
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function converts dispersion parameters to and from the customary
% engineering parameters of D, S and dS/dlambda to beta_2, beta_3, beta_4
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% dispersion = 16;              % dispersion, in ps/nm/km
% dispersion_slope = 0.058;     % dispersion slope, in ps/nm2/km
% dispersion_curvature = 0;     % dispersion curvature, in ps/nm3/km
% dispersion_spec_frequency = CONSTANT.c/1550e-9;
% beta = conv_disp_d_beta([dispersion dispersion_slope dispersion_curvature],'to_beta','eng','si',dispersion_spec_frequency);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% in                input dispersion [real vector]
%
%                       if mode ='to_beta':
%
%                           in =[D dD d2D]
%                           where:
%                           D is the dispersion per unit length
%                           dD is the first order derivative of D with
%                           respect to wavelength (dispersion slope)
%                           d2D is the second order derivative of D with 
%                           respect to wavelength (dispersion curvature)
%
%                           The units depend on the units_in parameters
%                           units_in         'si'         'eng'
%                           D                s/m^2        ps/nm/km
%                           dD               s/m^3        ps/nm^2/km
%                           d2D              s/m^4        ps/nm^3/km
%
%                           Conversion rule:
%                           1 ps/nm^n/km= 10^(9n-15) s/m^(n+1)
%
%                           If the input array has less than 3 elements, it 
%                           is assumed it consists only of the [D] or 
%                           [D dD] and the missing elemnts are arbitrarily
%                           set to 0.
%
%                       if mode= 'from_beta':
%
%                           in= [beta_1 beta_2 beta_3 beta_4]
%                           where the beta's are the coefficients of the 
%                           Taylor expansion of the propagation constant 
%                           around the frequency defined by the input 
%                           parameter freq.
%
%                           The units depend on the units_in parameters
%                           units_in         'si'         'eng'
%                           beta_1           s/m          ps/km
%                           beta_2           s^2/m        ps^2/km
%                           beta_3           s^3/m        ps^3/km
%                           beta_4           s^4/m        ps^4/km
%
%                           Conversion rule: 
%                           1 s^n/m = 10^(12n+3) ps^n/km
%
%                           If the input array has less than 4 elements, 
%                           it is assumed it consists only of 
%                           [beta_1 beta_2 beta_3] or [beta_1 beta_2] or 
%                           [beta_1] and the missing elements are 
%                           arbitrariy set to 0.
%
%                           Note that the value of beta_1 is not meaningful
%                           and is not used in any of the calculations. 
%                           It is just here to follow our usual convention 
%                           on arrays containing beta's according to which
%                           the element of index i is equal to 
%                           beta_i. Setting beta_1 to 0 is just fine.
%
% mode              conversion mode [string]
%
%                       mode = 'to_beta'
%                           conversion from [D dD d2D] to
%                           [beta_1 beta_2 beta_3 beta_4]
%
%                       mode  ='from_beta'
%                           conversion from [beta_1 beta_2 beta_3 beta_4] 
%                           to [D dD d2D]
%
% units_in          specifies whether the input vector is provided in 'si'
%                       or 'eng' units according to the tables above 
%                       [string]
%
%                       units_in = 'si'
%                       units_in = 'eng'
%
% units_out         specifies whether the output vector is provided in 'si' 
%                       or 'eng' units according to the tables above 
%                       [string]
%
%                       units_out = 'si'
%                       units_out = 'eng'
%
% freq              frequency at which the values of [D dD d2D] or
%                       beta_1 beta_2 beta_3 beta_4] are calculated, in Hz
%                       [real scalar]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% out               converted dispersion parameters [real vector]
%
%                       if mode = 'from_beta':
%
%                           out = [D dD/dlambda d2D/dlambda2]
%
%                       if mode = 'to_beta':
%
%                           out = [beta_1 beta_2 beta_3 beta_4]
%
%                       Obviously beta_1 is not relevant but is just here
%                       so that the element of index i is equal to beta_i.
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% CONSTANT              essential physical constants [structure]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global CONSTANT


% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
if strcmp(mode,'to_beta')
    nelements = length(in);
    % Number of terms in the input array.
    % if keep_log.Verbose
        if nelements > 3
            disp('conv_disp_d_beta: elements of input array D over d^2D/dlambda^2 will be ignored.');
        elseif nelements < 3
            fprintf(1,'%s%i%s\n','conv_disp_d_beta: only ', nelements,' elements provided in input array D.');
        end
        % Display notification if number of elements in the input array is
        % different from 3.
    % end
    
    for ii = nelements+1:3
        in(ii) = 0;
    end
    % The input array is completed up to d2D in case fewer elements are
    % provided as input. In this way the general equations below for
    % beta_2, beta_3 and beta_4 can be used.
    
    if strcmp(units_in,'eng')
        % if keep_log.verbose
            disp('conv_disp_d_beta: values in input array are assumed to be expressed in ps/nm^i/km.');
            disp('conv_disp_d_beta: converting input array to SI units of s/m^(i+1).');
        % end
        for ii = 1:3
            in(ii) = in(ii)*10^(9*ii - 15);
        end        
    end
    % Converts the units of the betas from engineering to SI units if
    % needed.
    out = zeros(1,4);
    % Preallocate output array that will contain [beta_1 beta_2 beta_3
    % beta_4].
    lambda = CONSTANT.c/freq;
    % Wavelength in m.
    out(1) = 0;
    % beta_1 in s/m
    % beta_1 is arbitrarily set to 0. Only here to ensure that
    % out(i)=beta_i, i.e. easier to deal with the indices in the arrays
    % containing the beta's.  
    out(2) = -lambda^2/(2*pi*CONSTANT.c)*in(1);
    % beta_2 in s^2/m
    out(3) = (lambda^2/(2*pi*CONSTANT.c))^2*(2/lambda*in(1)+in(2));
    % beta_3 in s^3/m
    out(4) = -(lambda^4/(2*pi*CONSTANT.c)^3)*(6*in(1)+6*lambda*in(2)+lambda^2*in(3));
    % beta_4 in s^4/m

    if strcmp(units_out,'eng')
        % if keep_log.verbose
            disp('conv_disp_d_beta: converting output array to Eng units of ps^i/km.');
        % end
        for ii = 1:4
            out(ii) = out(ii)*10^(12*ii + 3);
        end
    end
    % Convert the output array from SI to engineering units if needed. 
    
    
elseif strcmp(mode,'from_beta')           
    nelements = length(in);
    % Number of beta terms in the input array.   
    % if keep_log.verbose
        if nelements > 4
            disp('conv_disp_d_beta: elements of input array beta over beta_4 will be ignored.');
        elseif nelements < 4
            fprintf(1,'%s%i%s\n','conv_disp_d_beta: only ', nelements,' elements provided in input array beta.');
        end
        % Display notification if number of elements in the input array is
        % different from 4.
    % end
    
    for ii = nelements + 1:4
        in(ii) = 0;
    end
    % The beta array is completed up to beta_4 in case fewer elements are
    % provided as input. In this way the general equations below for D, dD
    % and d2d can still be used
    
    if strcmp(units_in,'eng')
        % if keep_log.verbose
            disp('conv_disp_d_beta: values of input array beta_i are assumed to be expressed in ps^i/km.');
            disp('conv_disp_d_beta: converting input array beta_i to SI units of s^i/m.');
        % end
        for ii = 1:4
            in(ii) = in(ii)/10^(12*ii + 3);
        end        
    end
    % Converts the units of the betas from engineering to SI units if
    % needed.
    out = zeros(1,3);
    % Preallocate output array that will contain [D dD d2D]
    
    out(1) = -2*pi*freq^2/CONSTANT.c*in(2);
    % D from beta_2; D in s/m^2    
    out(2) = 4*pi*freq^3/CONSTANT.c^2*(in(2) + pi*freq*in(3));
    % S=dD/dlambda from beta_2 and beta_3; dD in s/m^3    
    out(3) = -4*pi*freq^4/CONSTANT.c^3*(3*in(2) + 6*pi*freq*in(3) + 2*pi^2*freq^2*in(4));
    % dS/dlambda=d^2D/dlambda^2 from beta_2, beta_3 and beta_4
    % d2D in s/m^4    

    if strcmp(units_out,'eng')
        % if keep_log.verbose
            disp('conv_disp_d_beta: converting output array to Eng units of ps/nm^i/km.');
        % end
        for ii = 1:3
            out(ii) = out(ii)/10^(9*ii - 15);
        end
    end
    % Convert the output array from SI to engineering units if needed.  
    
    
else
    error('conv_disp_d_beta: dispersion conversion mode not implemented.');
end



end
