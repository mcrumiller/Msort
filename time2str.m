function timestr = time2str(seconds)
%--------------------------------------------------------------------------
% time2str.m - Given a number in seconds, converts the time to hh:mm:ss
% format.  Returns a string.
%
% Usage: time2str();
%
% Input:  
% Output: 
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
timestr=sprintf('%02.0f:%02.0f:%02.0f',floor(seconds/3600),floor(mod(seconds,3600)/60),round(mod(seconds,60)));