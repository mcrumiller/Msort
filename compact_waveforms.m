function waveforms = compact_waveforms(waveforms,num_channels)
%--------------------------------------------------------------------------
% compact_waveforms.m - Takes a [W,C*T] set of concatenated waveforms and
% expands them back into their full [C,W,T] form.
%
% Usage: waveforms = compact_waveforms(waveforms,num_channels);
%
% Input:  waveforms                     * [W,C*T] matrix of waveforms
%         num_channels                  * number of separate channels
% Output: waveforms                     * [C,W,T] matrix of waveforms
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[W,len]=size(waveforms); C=num_channels;
if(mod(len,C)~=0), error('matrix size does not match up with num_channels'); end

C=num_channels; T=len/C;

waveforms = reshape(waveforms',[T C W]);
waveforms = permute(waveforms, [2 3 1]);