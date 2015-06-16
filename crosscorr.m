function [t,n]=crosscorr(x1, x2, stepsize, maxlag)
%CROSSCORR   generate cross-correlation function of two spike trains
%   [t,n] = CROSSCORR(x1,x2,stepsize,maxlag) generates a cross-correlation
%   vector of two spike trains x1 and x2, in increments of stepsize, and
%   running from -maxlag to maxlag.
%
%   If x1==x2, then CROSSCORR returns the autocorrelation function of x1.
%   
%   Input:  x1                * vector of spikes for spiketrain 1
%           x2                * vector of spikes for spiketrain 2
%           stepsize          * timestep increment
%           maxlag            * maximum temporal deviation
%   Output: t                 * vector of time values (lags)
%           n                 * value of cross correlation function at each
%                               time point (lag)
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%--------------------------------------------------------------------------
if(~exist('bin_width','var')), stepsize=1e-2; end; % 10 ms bins
if(~exist('lag','var')), maxlag=0.5; end % in seconds

% zero-pad vectors to same length, trim extraneous zeros
loc=find(x1~=0,1,'last'); x1=x1(1:loc);
loc=find(x2~=0,1,'last'); x2=x2(1:loc);
len1=length(x1); len2=length(x2);
len=max(len1,len2);
x1(end+1:len)=0; x2(end+1:len)=0;

% Calculate all pairs of differences between the two sets
s1 = repmat(x1,len,1);
s2 = repmat(x2',1,len);
diffs=s2-s1;
diffs=reshape(diffs,1,[]);
diffs(diffs<-maxlag | diffs>maxlag)=[];

% update: num points definted by number of points we have
num_pts=ceil(maxlag/stepsize);
t=(-stepsize*num_pts):stepsize:(stepsize*num_pts);
n=histc(diffs,t);