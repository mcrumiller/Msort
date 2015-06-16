function load_group(group_num, handles)
%LOAD_GROUP   Load data for a particular channel group from M-file
%   LOAD_GROUP() Loads all of the required data for Msort to operate from the selected matlab .m file. This function may
%   take a while with large recordings. Additionally, it determines whether or not data has been already been converted
%   from the associated .plx file and, if not, asks the user to convert.
%
%   This function also generates many GUI-related elements, which are dependent on properties of the loaded data.
%
%   Input: group_num..............channel group to load
%          handles................vector of Msort gui handles
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
clearvars -global -except palette;
global waveforms features data_v Ts palette; %#ok<NUSED>
set(handles.status_label,'String',sprintf('Loading group %g.  Please wait...',group_num));

set(handles.status_light,'BackgroundColor',palette.red); drawnow;
setappdata(handles.output,'group_num',group_num);
setappdata(handles.output,'original_timestamps',[]);
setappdata(handles.output,'timestamps',[]);
setappdata(handles.output,'idx',[]);

% swap visible group panels
set(handles.group_panel,'visible','off'); drawnow;
set(handles.group_cmd_panel,'visible','on');

% Group panel
group_panel=[handles.voltage_plot_button, handles.return_button, handles.rethreshold_button, handles.features_button, handles.stats_button];
set(handles.group_cmd_panel,'Title',sprintf('Group %g',group_num));

% Clusters Panel
clusters_panel=[handles.cluster_menu, handles.merge_button, handles.split_button, handles.restore_wfs_button,...
    handles.recluster_button handles.channel_sort_button handles.cluster_options_button handles.template_button handles.apply_templates_button];

% Parameters Panel
params_panel=[handles.cluster_options_button, handles.scale_checkbox, handles.param_button];

% Cleanup Panel
cleanup_panel = [handles.display_cleanup_button, handles.clear_misaligned_button handles.move_misaligned_button handles.errorbar_field, handles.errorbar_slider, handles.std_dev_label handles.cleanup_button handles.ISI_move_button];

set([group_panel clusters_panel params_panel cleanup_panel],'enable','off','visible','off');

setappdata(handles.output,'show_display_tools',false);
setappdata(handles.output,'big_mode',false);

% delete associated data, reset all panels
setappdata(handles.output,'timestamps',[]);
setappdata(handles.output,'idx',[]);
setappdata(handles.output,'normalized_waveforms',[]);
setappdata(handles.output,'t',[]);
clear_panel(handles.features_panel);
drawnow;

% load necessary data
plx_file = getappdata(handles.output,'plx_file');

% get channels from file
xml_file=getappdata(handles.output,'xml_file');
exp=xml2struct(xml_file);
setappdata(handles.output,'exp',exp);
if(~iscell(exp.experiment.shank)),exp.experiment.shank={exp.experiment.shank}; end
num_channels = length(exp.experiment.shank{group_num}.channel);
channels = zeros(1,num_channels);
disabled = false(1,num_channels);
thresholds = zeros(2,num_channels);
for i = 1:num_channels
    if(num_channels==1)
        channels(i)=str2double(exp.experiment.shank{group_num}.channel.Attributes.num);
    else
        channels(i) = str2double(exp.experiment.shank{group_num}.channel{i}.Attributes.num);
    end
    
    % should have two thresholds as of 2013.02.14; top is positive, bottom is negative
    if(num_channels==1)
        if(isfield(exp.experiment.shank{group_num}.channel,'threshold2'))
            thresholds(1,i)=str2double(exp.experiment.shank{group_num}.channel.threshold.Text);
            thresholds(2,i)=str2double(exp.experiment.shank{group_num}.channel.threshold.Text);
        else
            t=str2double(exp.experiment.shank{group_num}.channel.threshold.Text);
            if(t>0), thresholds(1,i)=t;
            else thresholds(2,i)=t;
            end
        end
        if(strcmp(exp.experiment.shank{group_num}.channel.enabled.Text,'false'))
            disabled(i)=true;
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
        if(strcmp(exp.experiment.shank{group_num}.channel{i}.enabled.Text,'false'))
        disabled(i)=true;
        end
    end
    
end

channels(disabled)=[]; thresholds(:,disabled)=[];
num_channels=length(channels);
setappdata(handles.output,'num_channels',num_channels);
setappdata(handles.output,'disabled',disabled);

% determine if our .mat file contains all processed channels
mat_file=getappdata(handles.output,'mat_file');
vars=whos('-file',mat_file);
vars={vars.name}; expr=regexp(vars,'data_(?<channel>\d*$)','names');
locs=cellfun(@(x) ~isempty(x),expr); vars_data=expr(locs);
num_channels_completed=length(vars_data);
channels_completed=zeros(1,num_channels_completed);
for i = 1:num_channels_completed
    channels_completed(i) = str2double(vars_data{i}.channel);
end
% compare to our loaded .xml file
channels_needed=channels;
found=0;
for i = 1:length(channels_needed)
    if(any(channels_completed==channels_needed(i)))
        found=found+1;
    end
end
% request plx2mat to user
if(found<length(channels_needed))
    reply = questdlg(sprintf('%g of %g total channels have been processed. Continue with the rest? This may take a while.',found,length(channels_needed)),'PLX -> MAT Processing','Yes','No','Yes');
    %choice = 'Yes';
    switch reply
        case 'Yes'
            % get electrode type
            electrode_type=exp.experiment.electrode.Text;
            switch plx_file(end-2:end)
                case 'plx'
                    plx2mat(plx_file,electrode_type,group_num);
                case 'abf'
                    abf2mat(plx_file);
            end
            % play beep
            if(get(handles.sound_checkbox,'Value')), stop_alert; end
            load_group(group_num,handles);
            return;
        case 'No'
            set(handles.group_cmd_panel,'visible','off');
            set(handles.group_panel,'visible','on');
            return;
    end
end

% determine what we can load
vars=whos('-file',mat_file);
names=cell(1,length(vars));
for i=1:length(vars)
    names{i}=vars(i).name;
end

% generate list of variables we need
varlist = {'Ts',...
    sprintf('idx_%g',group_num),...
    sprintf('original_timestamps_%g',group_num),...
    sprintf('num_PCs_%g',group_num),...
    sprintf('trigger_ch_%g',group_num);
    };
for i = 1:length(channels)
    varlist{end+1}=sprintf('data_%g',channels(i)); %#ok<AGROW>
end
std_list=cell(1,num_channels);
for i = 1:length(channels)
    std_list{i}=sprintf('std_%g',channels(i));
end

% determine which variables are in there
vars_exist=false(1,length(varlist));
for i = 1:length(varlist)
    if(any(strcmp(varlist{i},names)))
        vars_exist(i)=true;
    end
end

% load the data
varlist=varlist(vars_exist);

h=waitbar(0,'Loading variables','name','Loading Variables','color',palette.background);
set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
for i = 1:length(varlist)
    waitbar(i/length(varlist),h);
    load(mat_file,varlist{i});
end
delete(h);

if(~isempty(std_list))
    load(mat_file,std_list{:});
end

% join data channels
if(exist(sprintf('data_%g',channels(1)),'var'))
    s=sprintf('[data_%g''',channels(1));
    s2=sprintf('[std_%g',channels(1));
else
    s='[';
    s2='[';
end
for i = 2:num_channels
    if(exist(sprintf('data_%g',channels(i)),'var'))
        s=sprintf('%s data_%g''',s,channels(i));
        s2=sprintf('%s std_%g',s2,channels(i));
    end
end
s=sprintf('%s]',s);
s2=sprintf('%s]',s2);
cmd=sprintf('data_v = %s;',s); eval(cmd);
cmd=sprintf('std_devs = %s;',s2); eval(cmd);

% grab data squared
%data_squared=data_v.^2;

% save data into standard variables
%standard_names = {'Ts','idx','timestamps','num_PCs','trigger_ch','data_squared'};
standard_names = {'Ts','idx','timestamps','num_PCs','trigger_ch'};
vars_exist(end-length(channels)+1:end)=[]; %vars_exist(end+1)=true;
varlist(end-length(channels)+1:end)=[]; %varlist{end+1}='data';
%varlist{end+1}='data_squared';
standard_names=standard_names(vars_exist);

for i = 1:length(varlist)
    eval(sprintf('setappdata(handles.output,standard_names{i},%s);',varlist{i}));
end
setappdata(handles.output,'std_devs',std_devs);

setappdata(handles.output,'t_before',str2double(get(handles.t_before_field,'String')));
setappdata(handles.output,'t_after',str2double(get(handles.t_after_field,'String')));


% set channel colors
num_channels=length(channels);
channel_colors=get_channel_colors(num_channels);
setappdata(handles.output,'channel_colors',channel_colors);

% If we have num_PCs and normalized_waves, set those boxes as well
if(any(strcmp(standard_names,'normalized_waves')))
    eval(sprintf('normalized_waves=normalized_waves_%g;',group_num));
    if(normalized_waves), set(handles.scale_checkbox,'Value',get(handles.scale_checkbox,'Max'));
    else set(handles.scale_checkbox,'Value',get(handles.scale_checkbox,'Min'));
    end
end

% grab waveforms using timestamps
varname=sprintf('original_timestamps_%g',group_num);
if(~exist(varname,'var'))
    setappdata(handles.output,'timestamps',[]);
else
    eval(sprintf('setappdata(handles.output,''timestamps'',original_timestamps_%g);',group_num));
end
capture_waveforms(handles);

% set idx & num clusters, if exists
if(any(strcmp(standard_names,'idx')))
    eval(sprintf('idx=idx_%g;',group_num));
    if(any(idx))
        u=unique([0 idx]);
        num_clusters=length(u);
        setappdata(handles.output,'idx',idx);
    else
        num_clusters=1;
    end
else
    num_clusters=1;
end
setappdata(handles.output,'num_clusters',num_clusters);
setappdata(handles.output,'cluster_number',num_clusters-1);
selected_axes=false(1,num_clusters);
selected_axes(num_clusters)=true;
setappdata(handles.output,'selected_axes',selected_axes);

cluster_colors=get_cluster_colors(handles);
setappdata(handles.output,'cluster_colors',cluster_colors);

% update application
setappdata(handles.output,'channels',channels);
setappdata(handles.output,'thresholds',thresholds);
setappdata(handles.output,'display_cleanup_tools',false);

% re-enable buttons and GUI components
channel_colors=get_channel_colors(num_channels);
setappdata(handles.output,'channel_colors',channel_colors);

% populate the clusters panel with relevant data, depends on # channels
channel_tag_pos=get(handles.channel_tag,'position');
x_offset=channel_tag_pos(1)+channel_tag_pos(3); margin_right=channel_tag_pos(1);
ch_tag_width=(1-x_offset-margin_right)/num_channels;
ch_tag_height=channel_tag_pos(4);

% create channel tags and UI controls
cluster_ch_tags=zeros(1,num_channels);
num_PCs_field=zeros(1,num_channels);
num_PCs_slider=zeros(1,num_channels);
peak_checkbox=zeros(1,num_channels);
valley_checkbox=zeros(1,num_channels);
slope_checkbox=zeros(1,num_channels);
energy_checkbox=zeros(1,num_channels);
amplitude_checkbox=zeros(1,num_channels);
for i = 1:num_channels
    x=x_offset+(ch_tag_width)*(i-1);
    
    % channel tag
    y=channel_tag_pos(2);
    cluster_ch_tags(i)=uicontrol('Style','text','fontsize',16,'fontweight','bold',...
        'units','normalized','position',[x y ch_tag_width ch_tag_height],'HorizontalAlignment','right','String',...
        sprintf('Ch %g',channels(i)),'Parent',handles.cluster_panel,'BackgroundColor',palette.background,'ForegroundColor',channel_colors(i,:));
    
    % num PCs
    xPC=x+ch_tag_width;
    y=get(handles.num_PCs_tag,'position'); y=y(2);
    field_width=.05; slider_width=.03;
    sliderMax=3; sliderMin=0;
    num_PCs_field(i)=uicontrol('Style','Edit','fontsize',12,'units','normalized',...
        'position',[xPC-field_width-slider_width y field_width ch_tag_height],...
        'Parent',handles.cluster_panel,'String','3','BackgroundColor','white');
    num_PCs_slider(i)=uicontrol('Style','slider','units','normalized',...
        'position',[xPC-slider_width y slider_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',3,'Min',sliderMin,'Max',sliderMax,'Sliderstep',[1 1]/(sliderMax-sliderMin));
    
    % width for checkbox
    checkbox_width=.02;
    
    % Peak
    y=get(handles.peak_tag,'position'); y=y(2);
    peak_checkbox(i)=uicontrol('Style','Checkbox','units','normalized',...
        'position',[xPC-checkbox_width y checkbox_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',0,'fontsize',14,'BackgroundColor',palette.background);
    
    % Valley
    y=get(handles.valley_tag,'position'); y=y(2);
    valley_checkbox(i)=uicontrol('Style','Checkbox','units','normalized',...
        'position',[xPC-checkbox_width y checkbox_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',0,'fontsize',14,'BackgroundColor',palette.background);
    
    % Slope
    y=get(handles.slope_tag,'position'); y=y(2);
    slope_checkbox(i)=uicontrol('Style','Checkbox','units','normalized',...
        'position',[xPC-checkbox_width y checkbox_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',0,'fontsize',14,'BackgroundColor',palette.background);
    
    % Energy
    y=get(handles.energy_tag,'position'); y=y(2);
    energy_checkbox(i)=uicontrol('Style','Checkbox','units','normalized',...
        'position',[xPC-checkbox_width y checkbox_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',0,'fontsize',14,'BackgroundColor',palette.background);
    
    % amplitude
    y=get(handles.amplitude_tag,'position'); y=y(2);
    amplitude_checkbox(i)=uicontrol('Style','Checkbox','units','normalized',...
        'position',[xPC-checkbox_width y checkbox_width ch_tag_height],...
        'Parent',handles.cluster_panel,'Value',0,'fontsize',14,'BackgroundColor',palette.background);
end

% set handles
handles.num_PCs_field=num_PCs_field;
handles.num_PCs_slider=num_PCs_slider;
handles.peak_checkbox=peak_checkbox;
handles.valley_checkbox=valley_checkbox;
handles.slope_checkbox=slope_checkbox;
handles.energy_checkbox=energy_checkbox;
handles.amplitude_checkbox=amplitude_checkbox;

% makes for easy deleting
handles.cluster_controls=[num_PCs_field num_PCs_slider peak_checkbox valley_checkbox slope_checkbox energy_checkbox amplitude_checkbox];

% Create channel tags/checkboxes on the bottom left
ch_tags = zeros(1,num_channels);
ch_colorboxes = zeros(1,num_channels);
margin_y=.05; margin_x=.1;
ch_width=(1-3*margin_x)/2;
ch_height=(1-(num_channels+1)*margin_y)/max(8,num_channels);
for i = 1:num_channels
    y=1-(margin_y+ch_height)*i;
    ch_tags(i)=uicontrol('Style','Checkbox','units','Normalized','String',sprintf('%g',channels(i)),...
        'position',[margin_x y ch_width ch_height],'Parent',handles.channel_panel,...
        'BackgroundColor',palette.background,'ForegroundColor',palette.text);
    
    ch_colorboxes(i)=uicontrol('Style','text','units','Normalized','String','',...
        'position',[2*margin_x+ch_width y ch_width ch_height],'Parent',handles.channel_panel,...
        'BackgroundColor',channel_colors(i,:));
end
handles.ch_tags=ch_tags;
handles.ch_colorboxes=ch_colorboxes;
% update handles structure so we can pass it to the callback functions
guidata(handles.output,handles);

% Update call back functions
for i = 1:num_channels
    set(num_PCs_field(i),'CallBack',{@num_PCs_field_Callback,num_PCs_slider(i),handles});
    set(num_PCs_slider(i),'CallBack',{@num_PCs_slider_Callback,num_PCs_field(i),handles});
    set(peak_checkbox(i),'CallBack',{@cluster_checkbox_Callback,peak_checkbox,handles});
    set(valley_checkbox(i),'CallBack',{@cluster_checkbox_Callback,valley_checkbox,handles});
    set(slope_checkbox(i),'CallBack',{@cluster_checkbox_Callback,slope_checkbox,handles});
    set(energy_checkbox(i),'CallBack',{@cluster_checkbox_Callback,energy_checkbox,handles});
    set(amplitude_checkbox(i),'Callback',{@cluster_checkbox_Callback,amplitude_checkbox,handles});
    set(ch_tags(i),'CallBack',{@ch_tags_Callback,handles});
end

fix_feature_menu(handles);

selected_channels=false(1,num_channels);
setappdata(handles.output,'selected_channels',selected_channels);

% create dummy waveform axes
ax_waveforms=zeros(1,num_channels);
for i = 1:num_channels
    ax_waveforms(i)=axes('parent',handles.features_panel,'color','w','visible','off',...
        'xcolor',palette.text,'ycolor',palette.text,'xtick',[],'ytick',[],'clipping','off','linewidth',1.5); %#ok<LAXES>
    hold(ax_waveforms(i),'all');
end
setappdata(handles.output,'ax_waveforms',ax_waveforms);

ax_features=axes('position',[.225 0.04 1-.225 1-.04],'parent',handles.features_panel,...
    'color',palette.ax_background,'visible','off','linewidth',1.5,'box','on','xcolor',palette.text,'ycolor',palette.text);
set(ax_features,'ButtonDownFcn',{@select_features,handles});

hold(ax_features,'all');
setappdata(handles.output,'ax_features',ax_features);

% restore visibility
set([group_panel clusters_panel params_panel cleanup_panel handles.undo_button],'enable','on','visible','on');

% Update application data
setappdata(handles.output,'selected_features',[1 8 1]);

% Display voltage plot
setappdata(handles.output,'voltage_start',1);
plot_voltages(handles);

% set keyboard handling
set(handles.output,'WindowKeyPressFcn',{@key_press,handles},'WindowKeyReleaseFcn',@ctrl_release);

% update hidden clusters
fix_clusters(handles);
setappdata(handles.output,'clusters_hidden',true);
set(handles.hide_clusters_label,'visible','off');

% set status
set(handles.status_light,'BackgroundColor',palette.green);
set(handles.status_label,'String','Ready');