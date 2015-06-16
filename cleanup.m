function locs = cleanup(waveforms,std_dev)
%--------------------------------------------------------------------------
% cleanup.m - Finds large outliers, and returns an index vector to
% legitimate waveforms, with outliers removed.
%
% Usage: locs = cleanup(waveforms);
%
% Input:  waveforms                 * CxWxT matrix of waveforms
% Output: locs                      * 1xW logical vector of valid waveforms
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
waveforms=expand_waveforms(waveforms);
waveforms=zscore(waveforms);
locs=~any(waveforms>abs(std_dev) | waveforms<-abs(std_dev),2)';