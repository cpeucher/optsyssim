function streams_out = serial_parallel(stream_in,nstreams)
% Serial-to-parallel conversion
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function implements a digital serial-to-parallel converter.
% For instance:
% in = [1 2 3 4 5 6 7 8 9 10 11 12]
% nstreams = 3
% out(1,:) = [1 4 7 10]
% out(2,:) = [2 5 8 11]
% out(3,:) = [3 6 9 12]
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% streams_out = serial_parallel(bits_bin,nstreams); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% stream_in         input stream [vector]
%
% nstreams          number of outputs streams [integer]
%
%                       The length of the input stream should be an 
%                       integer multiple of the number of output streams.
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% streams_out       output streams [matrix]
%
%                       The data in the ith output is accessible 
%                       in streams_out(i,:).
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if mod(length(stream_in),nstreams) ~= 0
    error('serial_parallel: the number of symbols does not match the number of outputs.');
end

for iout = 1:nstreams
    streams_out(iout,:) = stream_in(iout:nstreams:length(stream_in));
end
% Serial to parallel conversion

end