function split_button_Callback(~, ~, handles)
%SPLIT_BUTTON_CALLBACK   Execute automatic clustering of spike subsets
%   SPLIT_BUTTON_CALLBACK runs RECLUSTER() on a subset of neurons selected by the user.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global features palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Reclustering...'); drawnow;

feature_mask=get_feature_mask(handles);
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);
selected_cell=getappdata(handles.output,'cluster_number');

% get waveforms of selected neurons
locs=false(1,length(idx));
for i = 1:length(selected_cell)
    locs(idx==selected_cell(i))=true;
end

% note: cluster 0 is noise
u=unique(idx);
num_clusters=length(u);
if(~any(u==0)), num_clusters=num_clusters+1; end

features_temp=features(locs,feature_mask);

cluster_type=get(handles.cluster_menu,'String');
cluster_type=cluster_type{get(handles.cluster_menu,'Value')};

% KlustaKwik
if(strcmp(cluster_type,'KlustaKwik'))
    %idx_2 = KlustaKwik(features_temp,'-MinClusters 2 -MaxClusters 6');
    idx_2 = KlustaKwik_java(features_temp,'-MinClusters 2 -MaxClusters 12','KlustaKwik.jar');
% K Medoids
elseif(strcmp(cluster_type,'K-Medoids'))
    idx_2=cell(1,3); energy=zeros(1,3);
    for i = 1:3
        [idx_2{i}, energy(i)] = kmedoids(double(features_temp'),i+1);
    end
    [~,best_loc]=min(energy);
    idx_2=idx_2{best_loc};
elseif(strcmp(cluster_type,'DBSCAN'))
    idx_2=dbscan(features_temp',5,(max(max(features_temp))-min(min(features_temp)))/10);
    idx_2(idx_2==-1)=0;
end

locs=find(locs);
num_new_cells=length(unique(idx_2));
if(num_new_cells==1)
    set(handles.status_label,'String','No new clusters found.');
    set(handles.status_light,'BackgroundColor',palette.green);
    return;
end

for i = 1:num_new_cells
    idx(locs(idx_2==i)) = num_clusters+i-1;
end
setappdata(handles.output,'idx',idx);

% fix clusters
fix_clusters(handles);
num_clusters=getappdata(handles.output,'num_clusters');
cluster_number=num_clusters-num_new_cells:num_clusters-1;
setappdata(handles.output,'cluster_number',cluster_number);
selected_axes=false(1,num_clusters);
selected_axes(num_clusters-num_new_cells+1:end)=true;
setappdata(handles.output,'selected_axes',selected_axes);

show_featurespanel(handles);

if(length(selected_cell)==1)
    cluster_string=sprintf('Cluster %g',selected_cell);
else
    cluster_string=sprintf('Clusters %g',selected_cell(1));
    for i = 2:length(selected_cell)-1
        cluster_string=sprintf('%s, %g',cluster_string,selected_cell(i));
    end
    cluster_string=sprintf('%s, and %g',cluster_string,selected_cell(end));
end

set(handles.merge_button,'Enable','on');

% play beep
if(get(handles.sound_checkbox,'Value')), stop_alert; end
set(handles.status_label,'String',sprintf('%s split into %g new clusters.',cluster_string,num_new_cells));
set(handles.status_light,'BackgroundColor',palette.green);