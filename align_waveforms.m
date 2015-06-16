function [waveforms, timestamps, t] = align_waveforms(waveforms, timestamps, t, search_width, trigger_ch)
%ALIGN_WAVEFORMS   align waveforms to peak
%   [waveforms, timestamps, t] = ALIGN_WAVEFORMS(waveforms, timestamps, t)
%   aligns a set of waveforms to their peaks and re-assigns the spike time
%   to the newly-aligned peak.
%
%   ALIGN_PEAKS(...,search_width) only aligns peaks within a window of
%   search_width time points.
%
%   ALIGN_PEAKS(...,trigger_ch) aligns peaks to the channels described in
%   trigger_ch.
%
%   Usage: [waveforms, timestamps, t] = align_waveforms(waveforms, timestamps, t);
%
%   Input:  waveforms            * CmWxT matrix of waveform values
%           timestamps           * CxW matrix of waveform timestamps
%           t                    * 1xT vector corresponding to time points
%           trigger_ch           * 1xT list of channels that triggered each
%                                  waveform.
%   Output: waveforms            * CmWxT matrix of aligned waveform values
%           timestamps           * CxW matrix of waveform timestamps
%           t                    * 1xT vector corresponding to time points
%
%   Written by Marshall Crumiller
%   email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
if(~exist('search_width','var') || isempty(search_width))
    search_width=3;
end

[C,W,T] = size(waveforms);
dt=mean(diff(t));

% determine where the peak is supposed to be
peakloc=find_index(t,0);
search_ind=(peakloc-search_width:peakloc+search_width);
peakloc=search_width+1; % midpoint of search_ind

% zoom in on current search space
wfs = waveforms(:,:,search_ind);

[~,min_locs]=max(abs(wfs),[],3);
min_locs=min_locs-peakloc;

% align by channel triggers, if we have them
if(exist('trigger_ch','var') && ~isempty(trigger_ch))
    % determine distance to peakloc for this channel
    locs=sub2ind(size(min_locs),trigger_ch,1:W);
else
    % determine closest peak to the peakloc
    [~,offset_locs]=min(abs(min_locs),[],1);
    locs=sub2ind(size(min_locs),offset_locs,1:W);
end

offsets=min_locs(locs);

max_offset=max(abs(offsets));

% generate new index offset matrix
ind=repmat(1:T,W,1);
ind=ind+repmat(offsets',1,T);
ind=ind(:,max_offset+1:end-max_offset);

wfs=zeros(C,W,size(ind,2),'single');
t=t(max_offset+1:end-max_offset);
for i = 1:W
    wfs(:,i,:)=waveforms(:,i,ind(i,:));
end
waveforms=wfs; clear wfs;

timestamps=timestamps+offsets*dt;

% re-sort
[timestamps,ind]=sort(timestamps);
waveforms=waveforms(:,ind,:);