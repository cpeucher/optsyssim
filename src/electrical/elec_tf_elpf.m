function tf = elec_tf_elpf(params,freq)
% Transfer functions of electrical low-pass filters
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function calculates the transfer functions of electrical low pass 
% filters of various types.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% params_elpf.type = 'bessel';%'butterworth';'gaussian;'none';'rc';'rectangular';'raised_cosine';'root_raised_cosine';
% params_elpf.order = 4;
% params_elpf.f3dB = 0.75*symbol_rate;
% params_elpf.roll_off = 0.15;
% params_elpf.symbol_rate = symbol_rate;
% params_elpf.hold_time = 1/symbol_rate/2;
% params.samples_per_symbol = nsamples_per_symbol;
% tf = elec_tf_elpf(params_elpf,frequency_array);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% params            low-pass filter parameters [structure]
%       
%                       params.type
%                           type of the filter [string]
%
%                           params.type = 'none'
%                           params.type = 'bessel'
%                           params,type = 'butterworth';
%                           params.type = 'gaussian'
%                           params.type = 'raised_cosine'
%                           params.type = 'rc'
%                           params.type = 'rectangular'
%                           params.type = 'root_raised_cosine'
%                           params.type = 'zero_order_hold'
%
%                       params.order
%                           filter order (for 'bessel', 'butterworth' and 
%                           'gaussian') [integer scalar]
%
%                       params.f3dB 
%                           cut-off frequency, in Hz [real scalar]
%
%                       params.roll_off
%                           roll-off factor (for 'raised_cosine' and
%                           'root_raised_cosine') [real scalar]
%
%                       params.symbol_rate
%                           symbole rate (for 'raised_cosine' and
%                           'root_raised_cosine') [real scalar]
%
%                       params.hold_time
%                           hold time (for 'zero_order_hold') [real scalar]
% 
%                       params.samples_per_symbol = nsamples_per_symbol
%                           number of samples (for 'zero_order_hold') 
%                           [integer scalar]
%
% freq              vector of frequency values at which the transfer
%                       function is calculated [real vector]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% tf                electrical filter complex transfer function 
%                       [complex vector]
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

nsamples = length(freq);
% Number of frequency values at which the transfer function should be
% calculated

switch params.type

    case 'none'

        % In case a filter module is present in a simulation set-up but one
        % whishes to use no actual filtering

        tf = ones(1,nsamples);

    case 'bessel'

        % Bessel filters of various orders

        switch params.order
            case 2
                p = 1i*freq*1.36165412871613/params.f3dB;
                tf = 3./(3 + 3*p + p.^2);
            case 3
                p = 1i*freq*1.75567236868121/params.f3dB;
                tf = 15./(15 + 15*p + 6*p.^2 + p.^3);
            case 4
                p = 1i*freq*2.11391767490422/params.f3dB;
                tf = 105./(105 + 105*p + 45*p.^2 + 10*p.^3 + p.^4);
            case 5
                p = 1i*freq*2.42741070215263/params.f3dB;
                tf = 945./(945 + 945*p + 420*p.^2 + 105*p.^3 + 15*p.^4 + p.^5);
            case 6
                p = 1i*freq*2.70339506120292/params.f3dB;
                tf = 10395./(10395 + 10395*p + 4725*p.^2 + 1260*p.^3 + 210*p.^4 + 21*p.^5 + p.^6);
            case 7
                p = 1i*freq*2.95172214703872/params.f3dB;
                tf = 135135./(135135 + 135135*p + 62370*p.^2 + 17325*p.^3 + 3150*p.^4 + 378*p.^5 + 28*p.^6 + p.^7);
            case 8
                p = 1i*freq*3.17961723751065/params.f3dB;
                tf = 2027025./(2027025 + 2027025*p + 945945*p.^2 + 270270*p.^3 + 51975*p.^4 + 6930*p.^5 + 630*p.^6 + 36*p.^7 + p.^8);
            case 9
                p = 1i*freq*3.39169313891166/params.f3dB;
                tf = 34459425./(34459425 + 34459425*p + 16216200*p.^2 + 4729725*p.^3 + 945945*p.^4 + 135135*p.^5 + 13860*p.^6 + 990*p.^7 + 45*p.^8 + p.^9);
            case 10
                p = 1i*freq*3.59098059456916/params.f3dB;
                tf = 654729075./(654729075 + 654729075*p + 310134825*p.^2 + 91891800*p.^3 + 18918900*p.^4 + 2837835*p.^5 + 315315*p.^6 + 25740*p.^7 + 1485*p.^8 + 55*p.^9 + p.^10);
        end

    case 'butterworth'

        % Butterworth filters of various orders

        p = zeros(1,params.order + 1);
        % The coefficients of the Butterwoth polynomial will be stored in p

        p(1) = 1;
        % Term of degree 0

        for ii = 2:params.order+1
            % Recurrence relation to compute the coefficient of the
            % Butterworth polynomial

            p(ii) = cos((ii - 2)*pi/2/params.order)/sin((ii - 1)*pi/2/params.order)*p(ii - 1);
            % Term of degree ii - 1

        end

        p = fliplr(p);
        % Matlab polynomial arrays go from coefficient of degree n to
        % coefficient of degree 0

        tf = 1./polyval(p,1i*freq/params.f3dB);
        % Evaluate response

    case 'gaussian'

        % Gaussian filter

        tf = exp(-log(sqrt(2))*(freq/params.f3dB).^(2*params.order));

    case 'rc'

        % RC filter

        tf = 1./(1 + 1i*(freq/params.f3dB));

    case 'rectangular'

        % Rectangular filter.

        tf = zeros(1,nsamples);
        tf(abs(freq) <= params.f3dB) = 1;

    case {'raised_cosine','root_raised_cosine'}

        % Raised cosine or root raised cosine filters

        ts = 1/params.symbol_rate;
        % Reciprocal of the symbol rate

        freq_range1 = freq(abs(freq) <= (1 - params.roll_off)/(2*ts));
        freq_range1i = freq(abs(freq) > (1 - params.roll_off)/(2*ts));
        freq_range2p = freq_range1i(freq_range1i <= (1 + params.roll_off)/(2*ts) & freq_range1i > 0);
        freq_range2m = freq_range1i(freq_range1i >= -(1 + params.roll_off)/(2*ts) & freq_range1i < 0);
        freq_range3p = freq(freq > (1 + params.roll_off)/(2*ts));
        freq_range3m = freq(freq <- (1 + params.roll_off)/(2*ts));
        % Divide the input frequency range in subdomains for both positive and
        % negative frequencies for piecewise definiton of the filters

        tf1= ts*ones(1,length(freq_range1));
        tf2p = ts/2*(1 + cos((pi*ts/params.roll_off)*(abs(freq_range2p) - ((1 - params.roll_off)/(2*ts)))));
        tf2m = ts/2*(1 + cos((pi*ts/params.roll_off)*(abs(freq_range2m) - ((1 - params.roll_off)/(2*ts)))));
        tf3p = zeros(1,length(freq_range3p));
        tf3m = zeros(1,length(freq_range3m));
        % Calculate the raised cosine transfer function over those domains

        tfrc = [tf3m tf2m tf1 tf2p tf3p];
        % Reassemble the raised cosine transfer function over those domains

        if strcmp(params.type,'raised_cosine')
            % Raised cosine filter

            tf = tfrc/max(tfrc);
            % The output is the (normalised) raised cosine transfer 
            % function

        elseif strcmp(params.type,'root_raised_cosine')
            % Root raised cosine filter

            tf = sqrt(tfrc)/max(sqrt(tfrc));
            % The output is the (normalised) root raised cosine transfer
            % function
        end
        % End of calculation of transfer function for raised cosine and root
        % raised cosine filters


    case 'zero_order_hold'

        % Zero-order hold filter

        tf = params.samples_per_symbol*func_sinc(pi*freq*params.hold_time).*exp(-1i*pi*freq*params.hold_time);
        % This will result in some oscillations in the output signal. Strictly
        % speaking one should use a proper windowing function
        % Should be followed by a low-pass filter that will take care of
        % suppressing these oscillations.

    otherwise
        error('elec_tf_elpf: the requested low-pass filter transfer function is not implemented.');

end
% End of switch over filter type

end