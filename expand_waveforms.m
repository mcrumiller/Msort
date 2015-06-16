function waveforms = expand_waveforms(waveforms)
%--------------------------------------------------------------------------
% expand_waveforms.m - Expands waveforms from a [C,W,T] matrix into a
% [W,T*C] matrix--essentially concatenating the waveforms from each
% channel. Dimensions:
%    -C Channels
%    -W Waveforms
%    -T time points per waveform
%
% Usage: waveforms = expand_waveforms(waveforms)
%
% Input:  waveforms                  * [C,W,T] matrix of waveforms.
% Output: waveforms                  * [C,W*T] matrix of waveforms.
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[C,~,T]=size(waveforms);
waveforms=permute(waveforms,[3 1 2]);
waveforms=reshape(waveforms,T*C,[])';