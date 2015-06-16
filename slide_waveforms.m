function slide_waveforms(handles,direction)
%SLIDE_WAVEFORMS   Move waveforms left or right to align peak
%   SLIDE_WAVEFORMS Shifts all waveforms to the left or to the right until the maximum peak on the current channel is
%   aligned with t=0.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms palette;

% we need to look for peaks to the right of the current peak, of opposite
% polarity
idx=getappdata(handles.output,'idx');
selected_waveforms=getappdata(handles.output,'selected_waveforms');
cluster_number=getappdata(handles.output,'cluster_number');
if(isempty(selected_waveforms))
    selected_waveforms=idx==cluster_number;
end

% get selected channels
selected_channels=getappdata(handles.output,'selected_channels');
if(sum(selected_channels)~=1 && length(selected_channels)>1)
    set(handles.status_label,'String','This action can only be performed if with a single channel selected.');
    set(handles.status_light,'BackgroundColor',palette.orange);
    return;
end
if(~any(selected_channels)),selected_channels=~selected_channels; end

t=getappdata(handles.output,'t');
zeroloc=find(t==0);

% move waveforms left
if(direction==-1)
    wfs = squeeze(waveforms(selected_channels,selected_waveforms,zeroloc:end));
    t=t(zeroloc:end);
    %peak_sin = sign(mean(wfs(:,1)));
% move waveforms right
else
    wfs = squeeze(waveforms(selected_channels,selected_waveforms,1:zeroloc));
    t=t(1:zeroloc);
    %peak_sin = sign(mean(wfs(:,end)));
end

%if(peak_sin==1), wfs=-wfs; end
wfs=abs(wfs);

% find peak, set a threshold of 1
[peaks,peaklocs] = find_peaks(wfs,std(wfs(:)));
if(isempty(peaks)),return; end
[~,locs]=max(peaks,[],2);
inds=sub2ind(size(peaklocs),1:size(peaklocs,1),locs');
peaklocs=peaklocs(inds);

timestamp_locs=peaklocs~=0;
peaklocs=peaklocs(timestamp_locs);
offset=t(peaklocs);
timestamps=getappdata(handles.output,'timestamps');
s=find(selected_waveforms);
timestamps(s(timestamp_locs))=timestamps(s(timestamp_locs))+offset;
setappdata(handles.output,'timestamps',timestamps);

% update trigger channels for these waveforms
trigger_ch=getappdata(handles.output,'trigger_ch');
trigger_ch(s(timestamp_locs))=find(selected_channels);
setappdata(handles.output,'trigger_ch',trigger_ch);

setappdata(handles.output,'selected_waveforms',[]);
capture_waveforms(handles);
show_featurespanel(handles);