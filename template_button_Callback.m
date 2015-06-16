function template_button_Callback(~, ~, handles)
%TEMPLATE_BUTTON_CALLBACK   Use preexisting clusters from previous recording as templates for sorting
%   TEMPLATE_BUTTON_CALLBACK uses clusters from a separate recording and uses template matching. All waveforms are
%   matched to current templates. It is up to the user to create new clusters, should they arise.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms palette;
pathname=getappdata(handles.output,'pathname');
[filename,pathname] = uigetfile('*.plx','Select a .plx file.',pathname);

if(~filename)
    set(handles.status_label,'String','Please load a .plx file');
    set(handles.status_light,'BackgroundColor',palette.green);
    return
end

% parse filename
expr = regexp(filename,'(?<basename>.*)\.plx','names');
basename = expr.basename;
%setappdata(handles.output,'basename',basename);
%setappdata(handles.output,'pathname',pathname);

% look for PLX file
mat_file = sprintf('%s%s.mat',pathname,basename);
% make sure mat file exists
if(~exist(mat_file,'file'))
    set(handles.status_light,'BackgroundColor',palette.red);
    set(handles.status_label,'String','Error: associated .mat file not found.');
end

% make sure .mat file contains what we need (waveforms)
group_num=getappdata(handles.output,'group_num');

names=whos('-file',mat_file);
names={names.name};
varname=sprintf('idx_%g',group_num);
if(~any(strcmp(names,varname)))
    set(handles.status_light,'BackgroundColor',palette.red);
    set(handles.status_label,'String','Error: selected .plx file has no associated clusters');
    return;
end

% ask to load same threshold levels
reply=questdlg('Would you like to apply the same threshold levels used in this file? (Note: you cannot Undo)','Load Thresholds','Yes','No','Yes');
new_thresh=false;
if(strcmp(reply,'Yes'))
    new_thresh=true;
    % grab xml file
    set(handles.status_light,'backgroundcolor',palette.red);
    set(handles.status_label,'String','Grabbing thresholds...'); drawnow;
    xml_file=sprintf('%s%s.xml',pathname,basename);
    exp = xml2struct(xml_file);
    disabled=getappdata(handles.output,'disabled');
    thresholds=zeros(2,length(disabled));
    for i = 1:length(disabled)
        if(disabled(i)), continue;
        else
            if(length(disabled)==1)
                if(isfield(exp.experiment.shank{group_num}.channel,'threshold2'))
                    thresholds(1,i) = str2double(exp.experiment.shank{group_num}.channel.threshold.Text);
                    thresholds(2,i) = str2double(exp.experiment.shank{group_num}.channel.threshold2.Text);
                else
                    t=str2double(exp.experiment.shank{group_num}.channel.threshold.Text);
                    if(t>0), thresholds(1,i)=t;
                    else thresholds(2,i)=t;
                    end
                end
            else
                if(isfield(exp.experiment.shank{group_num}.channel{i},'threshold2'))
                    thresholds(1,i) = str2double(exp.experiment.shank{group_num}.channel{i}.threshold.Text);
                    thresholds(2,i) = str2double(exp.experiment.shank{group_num}.channel{i}.threshold2.Text);
                else
                    t=str2double(exp.experiment.shank{group_num}.channel{i}.threshold.Text);
                    if(t>0), thresholds(1,i)=t;
                    else thresholds(2,i)=t;
                    end
                end
            end
        end
    end
    thresholds(:,disabled)=[];
    setappdata(handles.output,'thresholds',thresholds);
    
    % apply new thresholds
    set(handles.status_label,'String','Applying thresholds...'); drawnow;
    run_threshold_spikes(handles);
    
    % calculate features
    calculate_features(handles);
end

% load waveforms and clusters
load(mat_file,varname);
eval(sprintf('idx=%s;',varname));
varname=sprintf('original_waveforms_%g',group_num);
load(mat_file,varname);
eval(sprintf('waveforms_temp=%s;',varname));
varname=sprintf('t_%g',group_num);
load(mat_file,varname);
eval(sprintf('t_template=%s;',varname));
varname=sprintf('timestamps_%g',group_num);
load(mat_file,varname);
eval(sprintf('timestamps=%s;',varname));


% interpolate waveforms if we have to
if(get(handles.interp_waveforms_checkbox,'Value'))
    %trigger_ch=getappdata(handles.output,'trigger_ch');
    interp_factor=round(str2double(get(handles.interp_factor_field,'String')));
    [w,t_template]=interp_waveforms(waveforms_temp,t_template,interp_factor); %#ok<NODEF>
    [waveforms_temp,~,t_template]=align_waveforms(w,timestamps,t_template,interp_factor-1);
end

% determine mean waveforms
if(get(handles.scale_checkbox,'Value'))
    waveforms_temp=normalize_waveforms(waveforms_temp);
end
templates=generate_templates(waveforms_temp,idx); %#ok<NODEF>
clear waveforms_temp idx;

% If waveforms are different lengths, truncate one or the other for wf
% matching
t=getappdata(handles.output,'t');

% truncate either templates or waveforms for matching
t_center=find_index(t,0); t_template_center=find_index(t_template,0);
% truncate start
if(t_center<t_template_center)
    % truncate template
    offset=t_template_center-t_center;
    t_template(1:offset)=[]; templates(:,:,1:offset)=[];
elseif(t_center>t_template_center)
    % truncate t
    offset=t_center-t_template_center;
    t(1:offset)=[]; waveforms(:,:,1:offset)=[];
end

% truncate end
if(length(t)<length(t_template))
    % truncate template
    offset=length(t_template)-length(t);
    templates(:,:,end-offset+1:end)=[];
elseif(length(t)>length(t_template))
    % truncate t
    offset=length(t)-length(t_template);
    waveforms(:,:,end-offset+1:end)=[];
end

% now in position to match templates
set(handles.status_label,'String','Matching Templates...');

% if certain channels are selected, only template sort along those channels
% If certain channels are selected, only match along those dimensions
selected_channels=find(getappdata(handles.output,'selected_channels'));
if(~isempty(selected_channels))
    templates=templates(:,selected_channels,:);
    waveforms=waveforms(selected_channels,:,:);
end

[idx,dists] = match_spike_templates(templates, waveforms);

% Attempt to throw out spikes with low confidence
try
    idx2=kmeans(dists',2);
    group1=mean(dists(idx2==1)); group2=mean(dists(idx2==2));
    if(group1<group2),idx(idx2==2)=0; else idx(idx2==1)=0; end
    bad_dists=dists>(mean(dists)+2*std(dists));
    idx(bad_dists)=0;
catch
end

if(new_thresh), setappdata(handles.output,'idx_old',idx); end
setappdata(handles.output,'idx',idx);
fix_clusters(handles);
capture_waveforms(handles);

set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor',palette.green);

% display features
show_featurespanel(handles);