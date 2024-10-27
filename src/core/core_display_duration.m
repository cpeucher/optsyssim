function core_display_duration(start_time,end_time)
% Display execution time of a simulation
%
% -------------------------------------------------------------------------
% DESCRIPTION:
% -------------------------------------------------------------------------
% This function displays the duration between two specified events in days
% hours, minutes, seconds.
%
% -------------------------------------------------------------------------
% FUNCTION CALL:
% -------------------------------------------------------------------------
% start_time = datetime("now");
% end_time = datetime("now");
% core_display_duration(start_time,end_time); 
%
% -------------------------------------------------------------------------
% INPUTS:
% -------------------------------------------------------------------------
% start_time        starting date [datetime data type]
%
% end_time          ending date [datetime data type]
%
% -------------------------------------------------------------------------
% OUTPUTS:
% -------------------------------------------------------------------------
% 
%
% -------------------------------------------------------------------------
% GLOBAL:
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

duration = seconds(end_time - start_time);
% Calculate duration in seconds between the start and stop of the
% simulation

ndays = floor(duration/86400);
nhours = floor((duration - ndays*86400)/3600);
nmin = floor((duration - ndays*86400 - nhours*3600)/60);
nsec = duration - ndays*86400 - nhours*3600 - nmin*60;
% Convert duration to number of days, minutes, hours and seconds 

fprintf('\n\n%s%s\n\n','Simulation ended on ',end_time);
% Display end date of the simulation
fprintf('%s%i%s%i%s%i%s%3.2f%s\n\n\n','Duration: ',ndays,' days ',nhours,' hours ',nmin,' minutes ',nsec,' seconds');
% Display duration of the simulation

end