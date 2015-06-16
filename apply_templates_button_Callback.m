function apply_templates_button_Callback()
%APPLY_TEMPLATES_BUTTON_CALLBACK   Sort spikes according to current templates
%   [] = APPLY_TEMPLATES_BUTTON_CALLBACK() 
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms palette;
set(handles.status_light,'backgroundcolor',palette.red);
set(handles.status_label,'String','Matching templates...'); drawnow;

% match the original waveforms
%waveforms=getappdata(handles.output,'waveforms');
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);
templates=generate_templates(waveforms,idx);

% If certain channels are selected, only match along those dimensions
selected_channels=find(getappdata(handles.output,'selected_channels'));
if(~isempty(selected_channels))
    templates=templates(:,selected_channels,:);
    waveforms_temp=waveforms(selected_channels,:,:);
else
    waveforms_temp=waveforms;
end

[idx,dists] = match_spike_templates(templates, waveforms_temp);

% Attempt to throw out spikes with low confidence
%{
if(length(dists)>1)
    idx2=kmeans(dists,2);
    group1=mean(dists(idx2==1)); group2=mean(dists(idx2==2));
    if(group1<group2),idx(idx2==2)=0; else idx(idx2==1)=0; end
    bad_dists=dists>(mean(dists)+2*std(dists));
    idx(bad_dists)=0;
end
%}

setappdata(handles.output,'idx',idx);
fix_clusters(handles);
set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor',palette.green);

% display features
show_featurespanel(handles);