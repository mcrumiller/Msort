function recluster_button_Callback(~, ~, handles)
%RECLUSTER_BUTTON_CALLBACK   Cluster neurons based on selected waveform features
%   RECLUSTER_BUTTON_CALLBACK Automatically clusters neurons based on the features selected out by feature_mask.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global features palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Clustering data...');
drawnow;

idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);

feature_mask=get_feature_mask(handles);
features_temp=features(:,feature_mask);

% if we want fast sort, only use a subset of waveforms
fastsort=get(handles.fastsort_checkbox,'Value');
if(fastsort)
    max_wfs=str2double(get(handles.fastsort_field,'String'));
    if(max_wfs>size(features,1))
        fastsort=false;
    else
        ind=unique(round(linspace(1,size(features,1),max_wfs)));
        features_temp=features_temp(ind,:);
    end
end

% KlustaKwik Algorithm
cluster_type=get(handles.cluster_menu,'String');
cluster_type=cluster_type{get(handles.cluster_menu,'Value')};

if(strcmp(cluster_type,'KlustaKwik'))
    %idx_tmp = KlustaKwik(features_temp,'-MinClusters 2 -MaxClusters 12');
    idx_tmp = KlustaKwik_java(features_temp,'-MinClusters 2 -MaxClusters 12','KlustaKwik.jar');
elseif(strcmp(cluster_type,'K-Medoids'))
    % K Medoids
    idx_tmp=cell(1,8); energy=zeros(1,8); energy(1)=inf;
    for i = 2:8
        [idx_tmp{i}, energy(i)] = kmedoids(double(features_temp'),i);
    end
    [~,best_loc]=min(energy);
    idx_tmp=idx_tmp{best_loc};
elseif(strcmp(cluster_type,'DBSCAN'))
    idx_tmp=dbscan(features_temp,5,(max(max(features_temp))-min(min(features_temp)))/10);
    idx_tmp(idx_tmp==-1)=0;
end

% note: cluster 0 is noise
u=unique(idx_tmp); u(u==0)=[];
num_clusters = length(u)+1;
setappdata(handles.output,'num_clusters',num_clusters);

cluster_colors=get_cluster_colors(handles);
setappdata(handles.output,'cluster_colors',cluster_colors);

idx=zeros(1,size(features,1),'uint8');
if(fastsort)
    idx(ind)=idx_tmp;
    setappdata(handles.output,'idx',idx);
    apply_templates_button_Callback([],[],handles);
else
    idx=idx_tmp;
    setappdata(handles.output,'idx',idx);
    % update plots
    show_featurespanel(handles);
end

% play beep
if(get(handles.sound_checkbox,'Value')), stop_alert; end

set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);