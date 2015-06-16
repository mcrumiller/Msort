function show_featurespanel(handles)
%SHOW_FEATURESPANEL   Display features panel when user clicks on "Features" button
%   SHOW_FEATURESPANEL(handles) Loads a view of the current feature spanel.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global palette;
setappdata(handles.output,'panel_selected','features');
global features;
clear_panel(handles.features_panel);
timestamps=getappdata(handles.output,'timestamps');
if(isempty(timestamps))
    set(handles.status_label,'String','No waveforms at this threshold.');
    set(handles.status_light,'BackgroundColor',palette.orange);
    return;
end
fix_feature_axis_locations(handles);

if(isempty(features))
    calculate_features(handles);
end

% Update waveforms and clusters
generate_merge_buttons(handles);
plot_mean_waveforms(handles);
set(handles.feature_control_panel,'visible','on');
cluster_number=getappdata(handles.output,'cluster_number');
num_clusters=getappdata(handles.output,'num_clusters');
if(isempty(cluster_number) || cluster_number(1)<0 || cluster_number(1)>num_clusters-1)
    cluster_number=0;
    setappdata(handles.output,'cluster_number',cluster_number);
end

% set axis frame
num_clusters=getappdata(handles.output,'num_clusters');
selected_axes=false(1,num_clusters); selected_axes(cluster_number+1)=true;
frame_handles=getappdata(handles.output,'frame_handles');
set(frame_handles(selected_axes),'visible','on');
setappdata(handles.output,'selected_axes',selected_axes);
plot_features(handles);
show_all_waveforms(handles);

fix_feature_axis_locations(handles);
fix_featurepanel_display(handles);

set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor',palette.green);