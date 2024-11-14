function [fflag,message] = isoptical(sig)
% Test if a signal is optical
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function tests if a signal is optical, according to the definition
% used in optsyssim.
% The criteria are:
% - the signal is a structure
% - The structure possesses the fields 'x' and 'y'
% - The fields 'x' and 'y' have the same length as time_array
% If any of these conditions is not verified, sig cannot be an optical
% signal structure.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% fflag = isoptical(sig);
% [fflag,message] = isoptical(sig);
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% sig               input signal to test
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% fflag             result of the test [logical]
%
% message           reason why the test has failed [string]  
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% time_array            time samples, in s [real vector]
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

global time_array


if ~isstruct(sig)

    fflag = false;
    message = 'isoptical: signal is not a structure';

else

    if ~and(isfield(sig,'x'),isfield(sig,'y'))

        fflag = false;
        message = 'isoptical: field x or y missing in input structure';

    else 

        nsamples_x = length(sig.x);
        nsamples_y = length(sig.y);
        nsamples = length(time_array);

        if (nsamples_x ~= nsamples) || (nsamples_y ~= nsamples)

            fflag = false;
            message = 'isoptical: field x or y does not have the expected length';

        else

            fflag = true;
            message = 'isoptical: this could well be an optical signal';

        end

    end
end

end