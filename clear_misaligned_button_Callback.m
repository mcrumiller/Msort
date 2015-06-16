function clear_misaligned_button_Callback(~, ~, handles)
%CLEAR_MISALIGNED_BUTTON_CALLBACK   Send waveforms that are not peak-aligned to selected to channels to junk
%   [] = CLEAR_MISALIGNED_BUTTON_CALLBACK() 
%   
%   Input:  
%   Output: 
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms;
%waveforms=getappdata(handles.output,'waveforms');

% note: This isn't the actual channel #, just the index
selected_channels=getappdata(handles.output,'selected_channels');
if(~any(selected_channels)),selected_channels=~selected_channels; end
current_channel=find(selected_channels);

idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);

selected_axes=getappdata(handles.output,'selected_axes');
cluster_number=find(selected_axes)-1;
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);

if(isempty(cluster_number))
    cluster_locs=1:size(waveforms,2);
else
    c=false(1,size(waveforms,2));
    for i = 1:length(cluster_number)
        c(idx==cluster_number(i))=true;
    end
    cluster_locs=find(c);
end

t=getappdata(handles.output,'t');
center_loc=find_index(t,0);

waveforms_temp=waveforms(current_channel,cluster_locs,:);
minlocs=zeros(size(waveforms_temp,2),length(current_channel));
maxlocs=zeros(size(waveforms_temp,2),length(current_channel));

% find bad peaks
for i = 1:length(current_channel)
    [~,minlocs(:,i)]=min(waveforms_temp(i,:,:),[],3);
    [~,maxlocs(:,i)]=max(waveforms_temp(i,:,:),[],3);
end
badlocs=all(minlocs~=center_loc,2) & all(maxlocs~=center_loc,2);
idx(cluster_locs(badlocs))=0;
setappdata(handles.output,'idx',idx);

% update display
fix_clusters(handles);
show_featurespanel(handles);