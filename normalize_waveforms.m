function waveforms = normalize_waveforms(waveforms)
%--------------------------------------------------------------------------
% normalize_waveforms.m - Normalizes a set of waveforms, so all vectors are
% unit length.
%
% Usage: waveforms = normalize_waveforms(waveforms);
%
% Input:  waveforms              * CxWxT matrix of waveforms
%                                * C channels, W waveforms, T timepoints
% Output: waveforms              * normalized waveforms
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
len=sqrt(sum(sum(waveforms.^2,3),1));
len = repmat(len,[size(waveforms,1) 1 size(waveforms,3)]);
waveforms = waveforms./len;