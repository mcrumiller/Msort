function varargout = Msort(varargin)
% MSORT MATLAB code for Msort.fig
%      MSORT, by itself, creates a new MSORT or raises the existing
%      singleton*.
%
%      H = MSORT returns the handle to a new MSORT or the handle to
%      the existing singleton*.
%
%      MSORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSORT.M with the given input arguments.
%
%      MSORT('Property','Value',...) creates a new MSORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Msort_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Msort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Msort

% Last Modified by GUIDE v2.5 17-Dec-2014 15:32:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Msort_OpeningFcn, ...
    'gui_OutputFcn',  @Msort_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Msort is made visible.
function Msort_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Msort (see VARARGIN)
p=strrep(path,'\','/');
contains_Plex=regexp(p,'Plexon/mexPlex', 'once');
if(isempty(contains_Plex))
    Plexon_folder = sprintf('%s\\Plexon',pwd);
    mexPlex_folder = sprintf('%s\\Plexon\\mexPlex',pwd);
    addpath(Plexon_folder,mexPlex_folder);
end

% determine default data directory
Msort_path=which('Msort.m');
p=regexp(Msort_path,strrep(sprintf('(?<path>.*%c).*.m',filesep),'\','\\'),'names');
param_file=sprintf('%sMsort.cfg',p.path);
params=load_params(param_file);
default_dir = params.default_file_dir;
setappdata(hObject,'default_dir',default_dir);

% Choose default command line output for Msort
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set_theme(handles); % set palette
warning off all;

% UIWAIT makes Msort wait for user response (see UIRESUME)
% uiwait(handles.fig_msort);



% --- set color theme
function set_theme(handles)
global palette;
theme = 'material';


% note: theme must implement the following values:
%   palette.background
%   palette.gridlines
%   palette.text
%   palette.header
%   palette.button
%   palette.button_text
%   palette.undo
%   palette.green
%   palette.red
%   palette.orange
%   palette.blue
%   palette.ax_background
switch(theme)
    case 'material'
        colors=import_colors(theme);
        palette.background=colors.light_grey;
        palette.gridlines=colors.grey;
        palette.text=colors.dark_grey;
        palette.header=colors.blue;
        palette.button=1-(1-colors.blue)/3;
        palette.button_text=colors.dark_grey;
        palette.undo=1-(1-colors.pink)/2;
        palette.view_button=1-(1-colors.indigo)/3;
        palette.green=colors.green;
        palette.red=colors.red;
        palette.orange=colors.orange;
        palette.blue=colors.indigo;
        palette.teal=colors.teal;
        palette.ax_background = 1-(1-colors.light_grey)/2;
        cluster_colors=[colors.red;colors.purple;colors.indigo;colors.teal;...
                        colors.brown;colors.blue_grey;colors.green;colors.orange;...
                        colors.pink;colors.green;colors.indigo;colors.amber;...
                        colors.cyan;colors.deep_purple;colors.light_green;...
                        colors.dark_grey;colors.pink];
    case 'dark'
        palette.background=[0 0 0];
        palette.gridlines=[.5 .5 .5];
        palette.text=[.7 .7 .7];
        palette.header=[.87 .87 .87];
        palette.button=[237 237 237]/255;
        palette.button_text=[.2 .2 .2];
        palette.undo=[236 214 214]/255;
        palette.view_button=[186 212 244]/255;
        palette.green=[0 .7 0];
        palette.red=[.7 0 0];
        palette.orange=[255 127 39]/255;
        palette.blue=[0 0 .7];
        palette.ax_background = [0 0 0];
        cluster_colors=distinguishable_colors(10,'k');
    case 'solarized'
        
end
% set background
set([handles.filename_panel handles.information_panel handles.features_panel handles.feature_control_panel handles.options_panel handles.cluster_panel ...
    handles.group_panel handles.group_cmd_panel handles.params_panel handles.clusters_panel handles.cleanup_panel...
    handles.interface_panel handles.channel_panel handles.status_panel],...
    'backgroundcolor',palette.background,'foregroundcolor',palette.header);
set([get(handles.group_panel,'children')' get(handles.group_cmd_panel,'children')' ...
    get(handles.clusters_panel,'children')' get(handles.cleanup_panel,'children')' ...
    get(handles.interface_panel,'children')' get(handles.information_panel,'children')'...
    get(handles.channel_panel,'children')' get(handles.options_panel,'children')' ...
    get(handles.filename_panel,'children')' get(handles.status_panel,'children')' ...
    get(handles.params_panel,'children')' get(handles.feature_control_panel,'children')'],...
    'backgroundcolor',palette.background,'foregroundcolor',palette.text);

% set buttons
set(findobj(handles.output,'style','pushbutton'),'backgroundcolor',palette.button,'foregroundcolor',palette.button_text);
set([handles.deselect_button handles.undo_button],'backgroundcolor',palette.undo,'foregroundcolor',palette.button_text);

set(handles.output,'color',palette.background);

% determine palette-based cluster colors
palette.cluster_colors=cluster_colors;


% --- Outputs from this function are returned to the command line.
function varargout = Msort_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
set(get(handles.output,'JavaFrame'),'Maximized',1);


% --- Executes on button press in group1_button.
function group1_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(1, handles);

% --- Executes on button press in group2_button.
function group2_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(2, handles);

% --- Executes on button press in group3_button.
function group3_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group3_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(3, handles);

% --- Executes on button press in group4_button.
function group4_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group4_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(4, handles);

% --- Executes on button press in group5_button.
function group5_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group5_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(5, handles);

% --- Executes on button press in group6_button.
function group6_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group6_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(6, handles);

% --- Executes on button press in group7_button.
function group7_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group7_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(7, handles);

% --- Executes on button press in group8_button.
function group8_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to group8_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_group(8, handles);


% --- Loads group data into the Features panel
function load_group(group_num, handles)
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



% controls the num_PCs slider & field
function num_PCs_slider_Callback(hObject,~,field,handles)
value=round(get(hObject,'Value'));

% if we're linked, make sure to update all of the values
if(value>get(hObject,'Max')), value=get(hObject,'Max'); end
if(value<get(hObject,'Min')), value=get(hObject,'Min'); end

if(get(handles.link_channels_checkbox,'Value'))
    set(handles.num_PCs_field,'String',sprintf('%2.0f',value));
    set(handles.num_PCs_slider,'Value',value);
else
    set(hObject,'Value',value);
    set(field,'String',sprintf('%2.0f',value));
end

% controls the num_PCs slider & field
function num_PCs_field_Callback(hObject,~,slider,handles)
value=str2double(get(hObject,'String'));
if(isnan(value)), value=1; end

if(value>get(slider,'Max')), value=get(hObject,'Max'); end
if(value<get(hObject,'Min')), value=get(hObject,'Min'); end

if(get(handles.link_channels_checkbox,'Value'))
    set(handles.num_PCs_field,'String',sprintf('%2.0f',value));
    set(handles.num_PCs_slider,'Value',value);
else
    set(slider,'Value',value);
    set(hObject,'String',sprintf('%2.0f',value));
end

function cluster_checkbox_Callback(hObject,~,row,handles)
linked=get(handles.link_channels_checkbox,'Value');
if(linked)
    set(row,'Value',get(hObject,'Value'));
end


% -- Generate plot of voltages
function plot_voltages(handles)
global Ts data_v data_spk data_t window_size y_offset voltage_plots spike_plots threshold_plots ax palette;
clear_panel(handles.features_panel);
set(handles.feature_control_panel,'visible','off');
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Plotting voltage traces...'); drawnow;

channels=getappdata(handles.output,'channels');
num_channels=length(channels);

% Display voltage plots
Fs=getappdata(handles.output,'Fs');
thresholds=getappdata(handles.output,'thresholds');

% Get rid of non-spikes from voltage trace
if(isempty(data_v))
    %if(get(handles.square_voltage_checkbox,'Value'))
    %    data_v=getappdata(handles.output,'data_squared')';
    %else
        data_v=(getappdata(handles.output,'data'))';
    %end
end

% create scrollbar
total_time=Ts(end)-Ts(1);
window_time=1; % in seconds
window_size=round(window_time*Fs);
if(~get(handles.square_voltage_checkbox,'Value'))
    std_devs = get(handles.voltage_range_slider,'Value');
    if(~isfinite(std_devs)),std_devs=15; end
    y_offset=repmat(num_channels-1:-1:0,window_size,1)*std_devs+std_devs;
else
    std_devs=15^2;
    y_offset=repmat(num_channels-1:-1:0,window_size,1)*std_devs+std_devs;
end


% create initial plot
start=1; stop=start+window_size-1;
t_start=Ts(start); t_stop=Ts(stop);
data_chunk=data_v(start:stop,:)+y_offset;

margin_left = .03; margin_right = .01;
margin_top = .02; margin_bot = .05;
plot_width = 1-(margin_left+margin_right);
plot_height = 1-(margin_top+margin_bot);

ax=axes('Parent',handles.features_panel,'position',[margin_left margin_bot plot_width plot_height],...
    'color',palette.ax_background,'ycolor',palette.text,'xcolor',palette.text,'box','on','xlim',[Ts(start) Ts(stop)],...
    'ylim',[0 (num_channels+.5)*std_devs],'ytick',(0:.5:num_channels-1)*std_devs,...
    'yticklabel',[],'linewidth',1.5);

xlabel('Time (s)','fontsize',16,'fontweight','bold','color',palette.text);
hold(ax,'all'); grid(ax,'on');
channel_colors=getappdata(handles.output,'channel_colors');
cluster_colors=getappdata(handles.output,'cluster_colors');

% make dim
timestamps=getappdata(handles.output,'timestamps');
if(~isempty(timestamps))
    set(ax,'ColorOrder',palette.gridlines);
else
    set(ax,'ColorOrder',channel_colors);
end


% plot data
voltage_plots=plot(Ts(start:stop),data_chunk(start:stop,:),'linesmooth','on');

% plot spike waveforms
if(~isempty(timestamps))
    idx=getappdata(handles.output,'idx');
    if(isempty(idx)),idx=zeros(1,length(timestamps),'uint8'); setappdata(handles.output,'idx',idx); end
    num_clusters=getappdata(handles.output,'num_clusters');
    
    t_before=getappdata(handles.output,'t_before');
    t_after=getappdata(handles.output,'t_after');
    
    % we need spike times and their indices into the data_chunk matrix
    locs=timestamps>t_start-t_after*1e-6 & timestamps<t_stop+t_before*1e-6;
    visible_timestamps=timestamps(locs)-Ts(start);
    visible_idx=idx(locs);
    
    % get vector and indexing for actual waveform shapes
    Fs=getappdata(handles.output,'Fs');
    pts_before=round(t_before*1e-6*Fs); pts_after=round(t_after*1e-6*Fs);
    t_ind=-pts_before:pts_after; t=t_ind/Fs;
    
    % we need a matrix of actual spike times, xC for each channel
    counts=zeros(1,num_clusters);
    [counts_exist,u]=count_uniques(visible_idx);
    counts(u+1)=counts_exist;
    max_spikenum=max(counts);
    if(max_spikenum==0),max_spikenum=1; end
    expanded_len=(pts_before+pts_after+1)*num_channels;
    t_template = [t nan]';

    
    % we have one spike matrix, each column is the same color
    max_len=max_spikenum*(expanded_len+num_channels);
    data_spk = nan(max_len,num_clusters);
    data_t = nan(max_len,num_clusters);

    for i = u+1
        spk_locs=round(visible_timestamps(visible_idx==i-1)*Fs)+1;
        
        x_offset=repmat(visible_timestamps(visible_idx==i-1),expanded_len+num_channels,1);
        t=repmat(t_template,num_channels,counts(i))+x_offset;
        
        % get waveform indices
        offset=repmat(t_ind',[1 counts(i)]);
        spk_locs2=repmat(spk_locs,[pts_before+pts_after+1 1]);
        
        spikes=spk_locs2+offset;
        spikes_allchannels=repmat(spikes,[1 1 num_channels]);
        
        
        offset=repmat((0:num_channels-1)*size(data_chunk,1),[counts(i) 1 pts_before+pts_after+1]);
        offset=permute(offset,[3 1 2]);
        
        spk_locs=permute(spikes_allchannels+offset,[1 3 2]);
        neglocs=spk_locs<1 | spk_locs>window_size*num_channels;
        spk_locs(neglocs)=1;
        wf_data=data_chunk(spk_locs);
        wf_data(neglocs)=nan;
        
        wf_data=reshape([wf_data;nan(1,num_channels,counts(i))],1,[]);
        t_data=reshape(t,1,[]);
        
        data_t(1:length(t_data),i)=t_data;
        data_spk(1:length(t_data),i)=wf_data;
    end
    
    set(ax,'ColorOrder',cluster_colors);
    spike_plots=plot(data_t,data_spk,'linesmooth','on');   
end


% plot thresholds
threshold_plots=zeros(2,num_channels);
if(get(handles.square_voltage_checkbox,'Value'))
    thresholds(2,:)=0;
    threshold_offset=repmat(num_channels-1:-1:0,2,1)*std_devs+std_devs;
else
    threshold_offset=repmat(num_channels-1:-1:0,2,1)*std_devs+std_devs;
end
zerolocs=thresholds==0;
thresholds=thresholds+threshold_offset; thresholds(zerolocs)=0;
bright_colors=1-(1-channel_colors)*.4;
for i = 1:num_channels
    if(thresholds(1,i)~=0)
        threshold_plots(1,i)=plot([Ts(1) Ts(end)],[thresholds(1,i) thresholds(1,i)],'--','color',bright_colors(i,:),'linewidth',1);
    end
    
    if(thresholds(2,i)~=0)
        threshold_plots(2,i)=plot([Ts(1) Ts(end)],[thresholds(2,i) thresholds(2,i)],'--','color',bright_colors(i,:),'linewidth',1);
    end
end
threshold_plots(zerolocs)=[];
uistack(threshold_plots(:),'bottom');

% Create channel tags
ch_tag_height=.03;
subplot_height=std_devs/(std_devs*(num_channels+1));
ch_tag_y = (num_channels-1:-1:0)*(subplot_height)+subplot_height/2-ch_tag_height+margin_bot+ch_tag_height/2+subplot_height/2;
ch_tag_width=margin_left-.013;
t=zeros(1,num_channels);
for i = 1:num_channels
    t(i)=uicontrol('Style','Text','String',sprintf('%g',channels(i)),'fontsize',20,...
        'fontweight','bold','ForegroundColor',channel_colors(i,:),'BackgroundColor',palette.background,'HorizontalAlignment','Right',...
        'position',[margin_left-ch_tag_width-.005 ch_tag_y(i) ch_tag_width ch_tag_height],...
        'Units','Normalized','Parent',handles.features_panel);
end

overlap=0.3;
minor_step=window_time*overlap/total_time; major_step=1/5;
slider=uicontrol('style','slider','Parent',handles.features_panel,'units','normalized','position',[margin_left margin_bot plot_width .02],...
    'backgroundcolor',palette.background,'SliderStep',[minor_step major_step],'Min',1,'Max',length(Ts)-window_size,'Value',1,...
    'CallBack',{@voltage_slider_Callback,handles});
setappdata(slider,'voltage_plots',voltage_plots);

voltage_start=getappdata(handles.output,'voltage_start');
if(~isempty(voltage_start))
    set(slider,'Value',voltage_start);
    voltage_slider_Callback(slider,[],handles);
end

setappdata(handles.output,'voltage_plot',ax);
setappdata(handles.output,'ch_voltage_tags',t);
setappdata(handles.output,'panel_selected','voltage');
setappdata(handles.output,'voltage_slider',slider);


% --- helper function for the voltage plot slider
function voltage_slider_Callback(hObject,~,handles)
global Ts data_v data_spk data_t window_size y_offset voltage_plots spike_plots ax;

% determine which range of data to display
start=max(floor(get(hObject,'Value')),get(hObject,'Min')); stop=start+window_size-1;
setappdata(handles.output,'voltage_start',start);
setappdata(handles.output,'voltage_stop',stop);

ind=start:stop;
if(stop>length(Ts)) % uh oh
    d=stop-length(Ts);
    start=start-d; stop=stop-d;
end

ysiz=size(y_offset,1); len=length(ind);
if(ysiz>len),y_offset=y_offset(1:len,:);
elseif(ysiz<len)
    extra=len-ysiz;
    y_offset(end+1:end+extra,:)=repmat(y_offset(1,:),extra,1);
end

% determine which data chunk to show
data_chunk=data_v(start:stop,:)+y_offset;
for i = 1:size(data_chunk,2)
    set(voltage_plots(i),'YData',data_chunk(:,i),'XData',Ts(start:stop));
end

% determine which spikes to show
timestamps=getappdata(handles.output,'timestamps');
num_channels=length(voltage_plots);
t_start=Ts(start); t_stop=Ts(stop);

if(~isempty(timestamps))
    idx=getappdata(handles.output,'idx');
    num_clusters=getappdata(handles.output,'num_clusters');
    
    t_before=getappdata(handles.output,'t_before');
    t_after=getappdata(handles.output,'t_after');
    
    % we need spike times and their indices into the data_chunk matrix
    locs=timestamps>t_start-t_after*1e-6 & timestamps<t_stop+t_before*1e-6;
    visible_timestamps=timestamps(locs)-Ts(start);
    %visible_timestamps=timestamps(locs);
    visible_idx=idx(locs);
    
    % get vector and indexing for actual waveform shapes
    Fs=getappdata(handles.output,'Fs');
    pts_before=round(t_before*1e-6*Fs); pts_after=round(t_after*1e-6*Fs);
    t_ind=-pts_before:pts_after; t=t_ind/Fs;
    
    % we need a matrix of actual spike times, xC for each channel
    counts=zeros(1,num_clusters);
    [counts_exist,u]=count_uniques(visible_idx);
    counts(u+1)=counts_exist;
    max_spikenum=max(counts);
    expanded_len=(pts_before+pts_after+1)*num_channels;
    t_template = [t nan]';

    % we have one spike matrix, each column is the same color
    max_len=max_spikenum*(expanded_len+num_channels);
    data_spk = nan(max_len,num_clusters);
    data_t = nan(max_len,num_clusters);

    for i = u+1
        % why do we need -1?
        spk_locs=round(visible_timestamps(visible_idx==i-1)*Fs)+1;
        
        x_offset=repmat(visible_timestamps(visible_idx==i-1),expanded_len+num_channels,1);
        t=repmat(t_template,num_channels,counts(i))+x_offset;
        
        % get waveform indices
        offset=repmat(t_ind',[1 counts(i)]);
        spk_locs2=repmat(spk_locs,[pts_before+pts_after+1 1]);
        
        spikes=spk_locs2+offset;
        spikes_allchannels=repmat(spikes,[1 1 num_channels]);
        
        
        offset=repmat((0:num_channels-1)*size(data_chunk,1),[counts(i) 1 pts_before+pts_after+1]);
        offset=permute(offset,[3 1 2]);
        
        spk_locs=permute(spikes_allchannels+offset,[1 3 2]);
        neglocs=spk_locs<1 | spk_locs>window_size*num_channels;
        spk_locs(neglocs)=1;

        wf_data=data_chunk(spk_locs);
        wf_data(neglocs)=nan;
        
        wf_data=reshape([wf_data;nan(1,num_channels,counts(i))],1,[]);
        t_data=reshape(t,1,[]);
        
        data_t(1:length(t_data),i)=t_data+Ts(start);
        data_spk(1:length(t_data),i)=wf_data;
    end
    
    for i = 1:length(spike_plots)
        set(spike_plots(i),'XData',data_t(:,i),'Ydata',data_spk(:,i));
    end
end

set(ax,'xlim',[Ts(start) Ts(stop)]);
drawnow;


% --- Executes on slider movement.
function numPCs_slider_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to numPCs_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
num_PCs = round(get(hObject,'Value'));
if(num_PCs<0)
    num_PCs=0;
elseif(num_PCs>8)
    num_PCs=8;
end

% sync up slider and field
set(handles.numPCs_field,'String',sprintf('%g',num_PCs));
set(hObject,'Value',num_PCs);
setappdata(handles.output,'num_PCs',num_PCs);

% --- Executes during object creation, after setting all properties.
function numPCs_slider_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to numPCs_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
minVal=0; maxVal=8; % these are defaults
dX=1/(maxVal-minVal);
set(hObject,'SliderStep',[dX dX]);

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function numPCs_field_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to numPCs_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numPCs_field as text
%        str2double(get(hObject,'String')) returns contents of numPCs_field as a double
num_PCs = round(str2double(get(hObject,'String')));
if(num_PCs<1)
    num_PCs=1;
elseif(num_PCs>8)
    num_PCs=8;
end
setappdata(handles.output,'num_PCs',num_PCs);

% sync up slider and field
set(hObject,'String',sprintf('%g',num_PCs));
set(handles.numPCs_slider,'Value',num_PCs);

% --- Executes during object creation, after setting all properties.
function numPCs_field_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to numPCs_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in loadfile_button.
function loadfile_item_Callback(~, ~, handles,filename,pathname)
% hObject    handle to loadfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Loading file...');

return_button_Callback([],[],handles);

if(~exist('filename','var') && ~exist('pathname','var'))
    default_dir=getappdata(handles.output,'default_dir');
    [filename,pathname] = uigetfile({'*.plx;*.abf'},'Select the a file.',default_dir);
end

if(~filename)
    set(handles.status_label,'String','Please load a .plx file');
    set(handles.status_light,'BackgroundColor',palette.green);
    return;
end

% parse filename
expr = regexp(filename,'(?<basename>.*)\.(plx|abf|mat)','names');
basename = expr.basename;
setappdata(handles.output,'basename',basename);
setappdata(handles.output,'pathname',pathname);

% look for PLX file
%plx_file = sprintf('%s%s.plx',pathname,basename);
plx_file=sprintf('%s%s',pathname,filename);
setappdata(handles.output,'plx_file',plx_file);

% look for XML file
xml_file = sprintf('%s%s.xml',pathname,basename);
if(~exist(xml_file,'file'))
    % generate XML file
    set(handles.status_label,'String','Generating associated .xml file...');
    switch(filename(end-2:end))
        case 'plx'
            generate_standard_Plexon_xml(plx_file);
        case 'abf'
            generate_standard_ABF_xml(plx_file);
    end
end
setappdata(handles.output,'xml_file',xml_file);

% look for MAT file
mat_file=sprintf('%s%s.mat',pathname,basename);
if(~exist(mat_file,'file'))
    save(mat_file,'-v7.3','');
end
setappdata(handles.output,'mat_file',mat_file);

% load struct into memory
exp = xml2struct(xml_file);
setappdata(handles.output,'exp',exp);


% get electrode type
setappdata(handles.output,'electrode',exp.experiment.electrode.Text);

% Populate Information fields
num_groups = length(exp.experiment.shank);
channel_groups = cell(1,num_groups);
group_nums = zeros(1,num_groups);
group_thresholds = cell(1,num_groups);

if(~iscell(exp.experiment.shank)),exp.experiment.shank={exp.experiment.shank}; end

for i = 1:num_groups
    % get group numbers
    group_nums(i) = str2double(exp.experiment.shank{i}.Attributes.num);
    num_channels = length(exp.experiment.shank{i}.channel);
    thresh = zeros(2,length(exp.experiment.shank{i}.channel));
    disabled=false(1,length(exp.experiment.shank{i}.channel));
    % get channel numbers
    
    % Note: this special case is due to the fact that xml2struct doesn't
    % place single items into cell arrays
    if(num_channels==1)
        channel_groups{i} = str2double(exp.experiment.shank{i}.channel.Attributes.num);
        
        % grab threshold
        if(isfield(exp.experiment.shank{group_nums(i)}.channel,'threshold2'))
            thresh(1,i) = str2double(exp.experiment.shank{group_nums(i)}.channel.threshold.Text);
            thresh(2,i) = str2double(exp.experiment.shank{group_nums(i)}.channel.threshold2.Text);
        else
            t=str2double(exp.experiment.shank{group_nums(i)}.channel.threshold.Text);
            if(t>0), thresh(1,i)=t;
            else thresh(2,i)=t;
            end
        end
        
        % check if channel is enabled
        if(strcmp(exp.experiment.shank{group_nums(i)}.channel.enabled.Text,'false'))
            disabled=true;
        end
    else
        % get channel number and threshold
        for c = 1:num_channels
            channel_groups{i}(c) = str2double(exp.experiment.shank{i}.channel{c}.Attributes.num);
            
            % grab both thresholds
            if(isfield(exp.experiment.shank{group_nums(i)}.channel{c},'threshold2'))
                thresh(1,i) = str2double(exp.experiment.shank{group_nums(i)}.channel{c}.threshold.Text);
                thresh(2,i) = str2double(exp.experiment.shank{group_nums(i)}.channel{c}.threshold2.Text);
            else
                t=str2double(exp.experiment.shank{group_nums(i)}.channel{c}.threshold.Text);
                if(t>0), thresh(1,i)=t;
                else thresh(2,i)=t;
                end
            end
            
            if(strcmp(exp.experiment.shank{i}.channel{c}.enabled.Text,'false'))
                disabled(c)=true;
            end
        end
    end
    channel_groups{i}(disabled)=[];
    thresh(:,disabled)=[];
    group_thresholds{i} = thresh;
end
Fs = str2double(exp.experiment.ADRate.Text);
num_total_channels = sum(cell2mat(cellfun(@length,channel_groups,'UniformOutput',false)));

% duration of experiment
exp_duration = str2double(exp.experiment.duration.Text);
setappdata(handles.output,'duration',exp_duration);

% save data to application
setappdata(handles.output,'num_groups',num_groups);
setappdata(handles.output,'channel_groups',channel_groups);
setappdata(handles.output,'group_nums',group_nums);
setappdata(handles.output,'group_thresholds',group_thresholds);
setappdata(handles.output,'Fs',Fs);
setappdata(handles.output,'exp_duration',exp_duration);
setappdata(handles.output,'num_total_channels',num_total_channels);

% update Information fields
set(handles.filename_label,'String',sprintf('%s%s.plx',pathname,basename));
set(handles.info_groups_label,'String',num2str(num_groups));
set(handles.info_channels_label,'String',num2str(num_total_channels));
set(handles.info_duration_label,'String',time2str(exp_duration));
set(handles.info_Fs_label,'String',sprintf('%g kHz',Fs/1e3));

% Enable the proper # of groups
strings=cell(1,num_groups+1); strings{1}='Select Group...';
for i = 1:num_groups, strings{i+1}=sprintf('Group %g',i); end
set(handles.group_menu,'String',strings);

% load event data
setappdata(handles.output,'events',[]);
load(mat_file,'events');
if(exist('events','var'))
    setappdata(handles.output,'events',events);
end
% note: gaps are in case user paused recording
setappdata(handles.output,'gaps',[]);
load(mat_file,'gaps');
if(exist('gaps','var'))
    setappdata(handles.output,'gaps',gaps);
end

set(handles.convert_button,'enable','on','visible','on');
set(handles.group_menu,'visible','on');
% Update status message
set(handles.status_label,'String','Select a group to process');


setappdata(handles.output,'panel_selected','none');
%load_group(1, handles);


% --- Executes on button press in return_button.
function return_button_Callback(~, ~, handles)
% hObject    handle to return_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear_panel(handles.features_panel);
clear_panel(handles.channel_panel);
set(handles.output,'WindowKeyPressFcn',@key_press,'WindowKeyReleaseFcn',@ctrl_release);
set(handles.group_cmd_panel,'visible','off');
set(handles.group_panel,'visible','on');
set(handles.feature_control_panel,'visible','off');
%set(handles.features_panel,'backgroundcolor',[.9412 .9412 .9412]);
set(handles.cluster_panel,'visible','off');
set(handles.group_menu,'Value',1);
if(isfield(handles,'cluster_controls'))
    delete(handles.cluster_controls(ishandle(handles.cluster_controls)));
end


% --- Executes on button press in rethreshold_button.
function rethreshold_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to rethreshold_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waveforms data_v Ts features palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Thresholding...');

group_num=getappdata(handles.output,'group_num');
        channels = getappdata(handles.output,'channels');
        exp=getappdata(handles.output,'exp');

reply=questdlg('Apply automatic Quiroga thresholding?','Automatic Threshold','Yes','No','Cancel','Yes');
%reply='No';
        
switch reply
    % Apply Quirogram thresholding of 5 * median{abs(data)/.6745}
    case 'Yes'
        set(handles.status_label,'String','Calculating Quiroga Thresholds...'); drawnow;
        %if(get(handles.square_voltage_checkbox,'Value'))
        %    data_tmp=getappdata(handles.output,'data_squared');
        %else
            data_tmp=data_v;
        %end
        M=5*median(abs(data_tmp)/.6745,1);
        thresholds(1,:)=M;
        thresholds(2,:)=-M;
    case 'No'
        clear_panel(handles.features_panel);
        set(handles.feature_control_panel,'visible','off');
        
        %if(get(handles.square_voltage_checkbox,'Value'))
        %    data_tmp=getappdata(handles.output,'data_squared');
        %else
            data_tmp=data_v;
        %end
        Fs=getappdata(handles.output,'Fs');
        plot_span=2;
        
        start=getappdata(handles.output,'voltage_start');
        if(isempty(start)),start=1; end
        stop=getappdata(handles.output,'voltage_stop');
        if(isempty(stop)),stop=start+plot_span*Fs-1; end
        ind=start:stop;
        
        % generate axes
        margin_left = .03; margin_right = .01;
        margin_top = .02; margin_bot = .05;
        plot_width = 1-(margin_left+margin_right);
        plot_height = 1-(margin_top+margin_bot);
        ax = axes('parent',handles.features_panel,'position',[margin_left margin_bot plot_width plot_height],...
            'ytick',[],'color',palette.ax_background,'xcolor',palette.text,'ycolor',palette.text,'linewidth',1.5,'box','on'); hold(ax,'on');
        
        % get thresholds
        num_channels=length(channels);
        thresholds=zeros(2,num_channels);
        if(get(handles.square_voltage_checkbox,'Value'))
            maxval=15^2; minval=-1;
        else
            maxval=15; minval=-15;
        end
        xgrid=linspace(Ts(1),Ts(ind(end)),9);
        ygrid=linspace(-double(maxval),double(maxval),11);
        channel_colors=getappdata(handles.output,'channel_colors');
        bright_colors=1-(1-channel_colors)*.4;
        for i = 1:num_channels
            cla(ax);
            plot(ax,Ts(ind),data_tmp(ind,i),'color',channel_colors(i,:)); axis([Ts(ind(1)) Ts(ind(end)) minval maxval]);
            xlabel('Time (s)','fontsize',16,'fontweight','bold','color',palette.text);
            minorgrid(xgrid,ygrid,ax);
            title_inside(sprintf('Channel %g',channels(i)),'fontsize',24,'fontname','Cambria','fontweight','bold','color',channel_colors(i,:));
            [~,y]=my_ginput(ax,bright_colors(i,:),false,handles.status_label);
            % user shift-clicked -- threshold above and below
            if(length(y)==2)
                thresholds(1,i)=abs(y(1));
                thresholds(2,i)=-abs(y(1));
                % only single threshold used
            else
                if(y>0),thresholds(1,i)=y;
                else thresholds(2,i)=y;
                end
            end
        end
        set(gcf,'pointer','arrow');
    otherwise
        return;
end

setappdata(handles.output,'thresholds',thresholds);
setappdata(handles.output,'normalized_waveforms',[]);
setappdata(handles.output,'idx_old',[]);
setappdata(handles.output,'num_clusters',[]);

% update struct
% generate channel list
disabled=getappdata(handles.output,'disabled');
ind=1;
if(~iscell(exp.experiment.shank)),exp.experiment.shank={exp.experiment.shank}; end
for i = 1:length(disabled)
    if(disabled(i)), continue;
    else
        if(length(disabled)==1)
            exp.experiment.shank{group_num}.channel.threshold.Text=sprintf('%g',thresholds(1,ind));
            exp.experiment.shank{group_num}.channel.threshold2.Text=sprintf('%g',thresholds(2,ind));
        else
            exp.experiment.shank{group_num}.channel{i}.threshold.Text=sprintf('%g',thresholds(1,ind));
            exp.experiment.shank{group_num}.channel{i}.threshold2.Text=sprintf('%g',thresholds(2,ind));
        end
    ind=ind+1;
    end
end
setappdata(handles.output,'exp',exp);

% remove old features and waveforms
waveforms=[]; features=[];
setappdata(handles.output,'normalized_waveforms',[]);

% remove clusters
setappdata(handles.output,'num_clusters',0);
setappdata(handles.output,'idx',[]);
plot_mean_waveforms(handles);
run_threshold_spikes(handles);
plot_voltages(handles);

% play beep
if(get(handles.sound_checkbox,'Value')), stop_alert; end
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);


% --- rethreshold spikes
function run_threshold_spikes(handles)
global waveforms features;
set(handles.status_label,'String','Applying thresholds...'); drawnow;

% rethreshold spikes
t_before=getappdata(handles.output,'t_before')*1e-6;
t_after=getappdata(handles.output,'t_after')*1e-6;

thresholds=getappdata(handles.output,'thresholds');

% threshold spikes
%if(get(handles.square_voltage_checkbox,'Value'))
%    data_squared=getappdata(handles.output,'data_squared');
%    [timestamps,trigger_ch] = M_threshold_spikes(handles, thresholds,t_before,t_after,data_squared);
%else
    [timestamps,trigger_ch] = M_threshold_spikes(handles, thresholds,t_before,t_after);
%end
if(isempty(timestamps))
    set(handles.status_label,'String','No Waveforms found; reduce threshold.');
    set(handles.status_light,'backgroundcolor',[255 178 0]/255);
    waveforms=[];
    setappdata(handles.output,'timestamps',[]);
    return;
end

% aligning may have aligned identical spikes; remove duplicates
setappdata(handles.output,'timestamps',timestamps);
setappdata(handles.output,'trigger_ch',trigger_ch);
capture_waveforms(handles);

num_waveforms=size(waveforms,2);
idx=zeros(1,num_waveforms,'uint8');
features=[];
setappdata(handles.output,'trigger_ch',trigger_ch);

setappdata(handles.output,'idx',idx);
setappdata(handles.output,'idx_old',idx);
fix_clusters(handles);


% --- calculates features from extracted waveforms
function calculate_features(handles,PCA_type)
global features palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Calculating Features...');
drawnow;

set(handles.status_label,'String','Extracting Features...'); drawnow;
t=getappdata(handles.output,'t');

% grab other features as well
num_PCs=get(handles.num_PCs_slider,'Max');
if(~iscell(num_PCs)),num_PCs={num_PCs}; end
params.PCA=max([num_PCs{:}]);
params.peaks=1;
params.valleys=1;
params.slope=1;
params.energy=1;
params.amplitude=1;

if(~exist('PCA_type','var')), PCA_type = 'all'; end
M_extract_features(t,params,handles,PCA_type);


% flatten features
features = permute(features,[3 1 2]);
features=reshape(features,size(features,1)*size(features,2),[])';

% Generate feature plots
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);



% --- Removes outliers from data
function remove_outliers(handles)
global features waveforms;

% find outliers and remove them
sigma_limit=str2double(get(handles.outliers_field,'String'));

%sigma_limit = 10; % 10 sigma limit
z = zscore(features,[],1);
ind = find(abs(z)>sigma_limit);
[I,~]=ind2sub(size(features),ind);
bad_locs = unique(I);

if(~isempty(I))
    % update features
    features(bad_locs,:)=[];
    
    % update interpolated data
    w=getappdata(handles.output,'waveforms_interp');
    w(:,bad_locs,:)=[];
    setappdata(handles.output,'waveforms_interp',w);
    w=getappdata(handles.output,'normalized_waveforms_interp');
    w(:,bad_locs,:)=[];
    setappdata(handles.output,'normalized_waveforms_interp',w);
    t=getappdata(handles.output,'timestamps_interp');
    t(bad_locs)=[];
    setappdata(handles.output,'timestamps_interp',t);
    
    % update regular data
    w=getappdata(handles.output,'original_waveforms');
    w(:,bad_locs,:)=[];
    setappdata(handles.output,'original_waveforms',w); clear w;
    w=getappdata(handles.output,'original_normalized_waveforms');
    w(:,bad_locs,:)=[];
    setappdata(handles.output,'original_normalized_waveforms',w); clear w;
    t=getappdata(handles.output,'original_timestamps');
    t(bad_locs)=[];
    setappdata(handles.output,'original_timestamps',t); clear t;
    
    % update current data
    waveforms(:,bad_locs,:)=[];
    w=getappdata(handles.output,'normalized_waveforms');
    w(:,bad_locs,:)=[];
    setappdata(handles.output,'normalized_waveforms',w);
    t=getappdata(handles.output,'timestamps');
    t(bad_locs)=[];
    setappdata(handles.output,'timestamps',t);
    
    % update trigger channels
    trigger_ch=getappdata(handles.output,'trigger_ch');
    trigger_ch(bad_locs)=[];
    setappdata(handles.output,'trigger_ch',trigger_ch);
    
    % update idx
    idx=getappdata(handles.output,'idx');
    idx_old=getappdata(handles.output,'idx_old');
    if(isempty(idx_old)), idx_old=zeros(1,length(idx),'uint8'); end
    idx(bad_locs)=[];
    idx_old(bad_locs)=[];
    
    setappdata(handles.output,'idx',idx);
    setappdata(handles.output,'idx_old',idx_old);
    
    features_button_Callback([],[],handles);
    set(handles.status_label,'String',sprintf('%g waveforms removed.',length(bad_locs)));
else
    set(handles.status_label,'String',sprintf('No waveforms removed.'));
    return;
end


% --- Executes on button press in features_button.
function features_button_Callback(~, ~, handles)
% hObject    handle to features_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
show_featurespanel(handles);


% --- plots features
function show_featurespanel(handles)
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



% -- Updates the features axis with features
function plot_features(handles)
global rotate features palette;
rotate = false;
if(isempty(features))
    calculate_features(handles);
end
timestamps=getappdata(handles.output,'timestamps');
features_temp=[features timestamps(:)];

set(handles.status_light,'backgroundcolor',palette.red);
set(handles.status_label,'String','Plotting features...'); drawnow;

% show features axis
ax_features=getappdata(handles.output,'ax_features');
feat_children=get(ax_features,'children'); if(iscell(feat_children)),feat_children=cell2mat(feat_children); end
delete(feat_children);
set(ax_features,'XTickMode','auto','YTickMode','auto','xlimmode','manual','ylimmode','manual','zlimmode','manual');

% Determine which features we have selected
selected_features=getappdata(handles.output,'selected_features');
if(isempty(selected_features))
    feature1=get(handles.feature1_menu,'Value');
    feature2=get(handles.feature2_menu,'Value');
    feature3=get(handles.feature3_menu,'Value');
    setappdata(handles.output,'selected_features',[feature1 feature2 feature3]);
else
    set(handles.feature1_menu,'Value',selected_features(1));
    set(handles.feature2_menu,'Value',selected_features(2));
    set(handles.feature3_menu,'Value',selected_features(3));
    feature1=selected_features(1);
    feature2=selected_features(2);
    feature3=selected_features(3);
end
feature1_list=get(handles.feature1_menu,'String');
feature2_list=get(handles.feature2_menu,'String');
feature3_list=get(handles.feature3_menu,'String');
cluster_colors=getappdata(handles.output,'cluster_colors');
num_clusters=getappdata(handles.output,'num_clusters');
feature_plots=zeros(1,num_clusters);

feature_names=getappdata(handles.output,'feature_names');

%feature_names=get(handles.feature1_menu,'String');


sphere_plots=zeros(1,num_clusters);
num_stds=get(handles.errorbar_slider,'Value');
wfs_in_sphere=false(size(features,1),num_clusters);
show_spheres=getappdata(handles.output,'display_cleanup_tools');
if(show_spheres),visible='on'; else visible='off'; end
idx=getappdata(handles.output,'idx');

if(feature2==length(feature2_list) || feature3==length(feature3_list))
    hist_mode=true;
    optimize_width=get(handles.hist_optimize_checkbox,'Value');
else
    hist_mode=false; hist_plots=[];
end
setappdata(handles.output,'hist_mode',hist_mode);
num_pts=size(features_temp,1);

% We have a 2D plot
if(feature3==1)
    % create histogram instead of scatterplot
    if(hist_mode)
        minFeat=min(features_temp(:,feature1)); maxFeat=max(features_temp(:,feature1));
        % optimize bin width
        if(~optimize_width)
            num_bins=round(20*log10(num_pts));
            edges=linspace(minFeat,maxFeat,num_bins);
            % double the X to make histogram
            x=[reshape(repmat(edges,2,1),1,[]) edges(end-1:-1:1)];
        end
        
        
        max_N=0;
        hist_plots=zeros(1,num_clusters);
        
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            feat1=features_temp(ind,feature1);
            if(optimize_width)
                [~,edges]=best_bin_width(feat1);
                %edge_diff=edges(2)-edges(1);
                %num_bins=ceil((maxFeat-minFeat)/edge_diff);
                %edges=linspace(minFeat,maxFeat,num_bins);
                % double the X to make histogram
                x=[reshape(repmat(edges,2,1),1,[]) edges(end-1:-1:1)];
            end
            N=histc(feat1,edges); N(end)=[]; max_N=max(max_N,max(N));
            N2=reshape(repmat(N',2,1),1,[]); N2=[0 N2 0]; %#ok<AGROW>
            hist_plots(i)=patch(x,[N2 zeros(1,length(N))],cluster_colors(i,:),'linesmooth','on','parent',ax_features,'edgecolor','none','FaceAlpha','flat');
        end
        h=ylabel(ax_features,'Count','color',palette.text,'fontsize',16,'fontweight','bold');
        span=edges(end)-edges(1);
        set(ax_features,'ButtonDownFcn',[],'xlim',[edges(1)-span*.01 edges(end)+span*.01],'ylim',[0 max_N*1.05]);
        alpha(hist_plots,.6);
    else
        MaxF=[-inf -inf]; MinF=[inf inf];
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1);
            feat2=features_temp(ind,feature2);
            feats=[feat1 feat2];
            MaxF=max([MaxF;max(feats,[],1)]);
            MinF=min([MinF;min(feats,[],1)]);
            
            % calculate ellipse
            feature_plots(i)=plot(ax_features,feat1,feat2,'.','markersize',1,'linesmooth','on','color',cluster_colors(i,:));
            if(i~=1 && size(feats,1)>2)
                [Ex,Ey]=ellipsoid_cov(mean(feats),num_stds,cov(feats),30);
                sphere_plots(i)=patch(Ex,Ey,cluster_colors(i,:),'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'parent',ax_features,'visible',visible);
                wfs_in_sphere(ind,i)=inhull(feats,get(sphere_plots(i),'Vertices'));
            end
        end
        alpha(sphere_plots(ishandle(sphere_plots) & sphere_plots~=0),0.2);
        set(ax_features,'xlim',[MinF(1) MaxF(1)],'ylim',[MinF(2) MaxF(2)],'ButtonDownFcn',{@select_features,handles});
        h=ylabel(ax_features,feature_names{feature2},'color',palette.text,'fontsize',16,'fontweight','bold');
    end
    view(ax_features,2); set(ax_features,'projection','orthographic');
    xlabel(ax_features,feature_names{feature1},'color',palette.text,'fontsize',16,'fontweight','bold');
    ylim=get(ax_features,'ylim');
    pos=get(h,'position'); set(h,'position',[pos(1) (ylim(2)-ylim(1))/3 pos(3)]);
    set(handles.output,'WindowButtonMotionFcn',[]);

% we have a 3D plot
else
    feature3=feature3-1;
    if(hist_mode)
        if(feature2==length(feature2_list)),feature2=feature1; end

        minFeat1=min(features_temp(:,feature1)); maxFeat1=max(features_temp(:,feature1));
        minFeat2=min(features_temp(:,feature2)); maxFeat2=max(features_temp(:,feature2));
        del = .001; % space between bars, relative to bar size
        if(~optimize_width)
            num_binsX=round(20*log10(num_pts)); num_binsY=num_binsX;
            edgesX=linspace(minFeat1,maxFeat1,num_binsX);
            edgesY=linspace(minFeat2,maxFeat2,num_binsY);
            %dX=edgesX(2)-edgesX(1); dY=edgesY(2)-edgesY(1);
            
            % Build x-coords for the eight corners of each bar.
            dX = edgesX(2)-edgesX(1);
            xx = [edgesX edgesX(end)+dX];
            xx = [xx(1:num_binsX)+del*dX; xx(2:num_binsX+1)-del*dX];
            xx = [reshape(repmat(xx(:)',2,1),4,num_binsX); NaN(1,num_binsX)];
            xx = [repmat(xx(:),1,4) NaN(5*num_binsX,1)];
            xx = repmat(xx,1,num_binsX);
            
            % Build y-coords for the eight corners of each bar.
            dY = edgesY(2)-edgesY(1);
            yy = [edgesY edgesY(end)+dY];
            yy = [yy(1:num_binsY)+del*dY; yy(2:num_binsY+1)-del*dY];
            yy = [reshape(repmat(yy(:)',2,1),4,num_binsY); NaN(1,num_binsY)];
            yy = [repmat(yy(:),1,4) NaN(5*num_binsY,1)];
            yy = repmat(yy',num_binsY,1);
        end
        
        max_N=0;
        hist_plots=zeros(1,num_clusters);
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2);
            if(optimize_width)
                [~,edgesX]=fdrule(feat1);
                [~,edgesY]=fdrule(feat2);
                %edgesX=best_bin_width(feat1); edgesY=best_bin_width(feat2);
                %edge_diffX=edgesX(2)-edgesX(1);
                %edge_diffY=edgesY(2)-edgesY(1);
                %num_binsX=ceil((maxFeat1-minFeat1)/edge_diffX);
                %num_binsY=ceil((maxFeat2-minFeat2)/edge_diffY);
                %edgesX=linspace(minFeat1,maxFeat1,num_binsX);
                %edgesY=linspace(minFeat2,maxFeat2,num_binsY);
                % Build x-coords for the eight corners of each bar.
                dX = edgesX(2)-edgesX(1);
                xx = [edgesX edgesX(end)+dX];
                xx = [xx(1:num_binsX)+del*dX; xx(2:num_binsX+1)-del*dX];
                xx = [reshape(repmat(xx(:)',2,1),4,num_binsX); NaN(1,num_binsX)];
                xx = [repmat(xx(:),1,4) NaN(5*num_binsX,1)];
                xx = repmat(xx,1,num_binsY);

                % Build y-coords for the eight corners of each bar.
                dY = edgesY(2)-edgesY(1);
                yy = [edgesY edgesY(end)+dY];
                yy = [yy(1:num_binsY)+del*dY; yy(2:num_binsY+1)-del*dY];
                yy = [reshape(repmat(yy(:)',2,1),4,num_binsY); NaN(1,num_binsY)];
                yy = [repmat(yy(:),1,4) NaN(5*num_binsY,1)];
                yy = repmat(yy',num_binsX,1);
            end
            n=hist3([feat1 feat2],{edgesX edgesY},'edgecolor','none','parent',ax_features);
            max_N=max(max_N,max(n(:)));
            % Build z-coords for the eight corners of each bar.
            
            zz = zeros(5*num_binsX, 5*num_binsY);
            n(n==0)=NaN;
            zz(5*(1:num_binsX)-3, 5*(1:num_binsY)-3) = n;
            zz(5*(1:num_binsX)-3, 5*(1:num_binsY)-2) = n;
            zz(5*(1:num_binsX)-2, 5*(1:num_binsY)-3) = n;
            zz(5*(1:num_binsX)-2, 5*(1:num_binsY)-2) = n;
            hist_plots(i)=surf(ax_features, xx, yy, zz,'facecolor',cluster_colors(i,:),'edgecolor','none','tag','hist3','linesmooth','off','facelighting','flat','Facealpha',0.6);
        end
        
        %alpha(hist_plots,.6);
        set(ax_features,'xlim',[edgesX(1) edgesX(end)],'ylim',[edgesY(1) edgesY(end)],'zlim',[0 max_N],'box','off');
    else
        MaxF=[-inf -inf -inf]; MinF=[inf inf inf];
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2); feat3=features_temp(ind,feature3);
            feats=[feat1 feat2 feat3];
            MaxF=max([MaxF;max(feats,[],1)]);
            MinF=min([MinF;min(feats,[],1)]);
            feature_plots(i)=plot3(ax_features,feat1,feat2,feat3,'.','markersize',1,'color',cluster_colors(i,:),'buttondownfcn',[],'hittest','off');
            
            % calculate ellipsoid
            if(size(feats,1)>2 && feature1 ~= feature2 && feature1 ~= feature3 && feature2 ~= feature3 && i~=1)
                [Ex,Ey,Ez]=ellipsoid_cov(mean(feats),num_stds,cov(feats),30);
                % sometimes we get a weird result with complex numbers with
                % zero imaginary component
                if(any(~isreal(Ex))),continue; end
                sphere_plots(i)=surf(Ex,Ey,Ez,'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'facelighting','gouraud','parent',ax_features,'visible',visible);
                [~,v]=surf2patch(Ex,Ey,Ez);
                wfs_in_sphere(ind,i)=inhull(feats,v);
            end
        end
        set(ax_features,'xlim',[MinF(1) MaxF(1)],'ylim',[MinF(2) MaxF(2)],'zlim',[MinF(3) MaxF(3)],'ButtonDownFcn',[]);
        zlabel(feature3_list{feature3},'fontsize',16,'fontweight','bold','color',palette.text);
        alpha(sphere_plots(ishandle(sphere_plots) & sphere_plots~=0),0.2);
    end
    
    xlabel(feature1_list{feature1},'fontsize',16,'fontweight','bold','color',palette.text);
    ylabel(feature2_list{feature2},'fontsize',16,'fontweight','bold','color',palette.text);
    
    view(ax_features,[15 45]);
    h=light('parent',ax_features);
    camlight(h,45,45);
    set(ax_features,'projection','perspective');
    
    h=rotate3d(ax_features);
    setappdata(handles.output,'rotate_handle',h);
    set(handles.output,'WindowButtonMotionFcn',{@enable_rotate,get(handles.features_panel,'position'),ax_features,h});
end
set(ax_features,'color',palette.ax_background,'xcolor',palette.gridlines,'ycolor',palette.gridlines,'fontweight','bold','zcolor',palette.gridlines);
grid(ax_features,'on');

setappdata(handles.output,'feature_plots',feature_plots);
setappdata(handles.output,'sphere_plots',sphere_plots);
setappdata(handles.output','wfs_in_sphere',wfs_in_sphere);
setappdata(handles.output,'hist_plots',hist_plots);

% highlight selected
highlight_selected_cluster(handles);
hide_clusters(handles);

set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor',palette.green);



function enable_rotate(f,~,featurepanel_pos,axHandle,rotate_handle)
global rotate;
currentPoint=get(f,'currentpoint');

axPos=get(axHandle,'position');
x_min=featurepanel_pos(1)+featurepanel_pos(3)*axPos(1);
x_max=featurepanel_pos(1)+featurepanel_pos(3)*(axPos(1)+axPos(3));
y_min=featurepanel_pos(2)+featurepanel_pos(4)*axPos(2);
y_max=featurepanel_pos(2)+featurepanel_pos(4)*(axPos(2)+axPos(4));

% we're outside
if(currentPoint(1)<x_min || currentPoint(1)>x_max || currentPoint(2)<y_min || currentPoint(2)>y_max)
    if(rotate)
        rotate=false;
        set(rotate_handle,'enable','off');
    end
    
% we're inside
else
    if(~rotate)
        rotate=true;
        set(rotate_handle,'enable','on');
    end
end



% --- hides clusters that aren't currently active
function hide_clusters(handles)
feature_plots=getappdata(handles.output,'feature_plots');

if(~isempty(feature_plots) && all(ishandle(feature_plots)))
    clusters_hidden=getappdata(handles.output,'clusters_hidden');
    selected_axes=getappdata(handles.output,'selected_axes');
    
    if(clusters_hidden)
        set(feature_plots(selected_axes & feature_plots~=0),'visible','on');
        set(feature_plots(~selected_axes & feature_plots~=0),'visible','off');
    else
        set(feature_plots(feature_plots~=0),'visible','on');
    end
    
    hist_plots=getappdata(handles.output,'hist_plots');
    if(any(hist_plots))
        valid_handles=ishandle(hist_plots) & hist_plots~=0;
        if(clusters_hidden)
            set(hist_plots(selected_axes & valid_handles),'visible','on');
            set(hist_plots(~selected_axes & valid_handles),'visible','off');
        else
            set(hist_plots(valid_handles),'visible','on');
        end
    end
end


% --- lights up dots in the selected cluster
function highlight_selected_cluster(handles)
feature_plots=getappdata(handles.output,'feature_plots');

if(~isempty(feature_plots) && all(ishandle(feature_plots)))
    selected_axes=getappdata(handles.output,'selected_axes');
    set(feature_plots(selected_axes & feature_plots~=0),'markersize',5);
    set(feature_plots(~selected_axes & feature_plots~=0),'markersize',1);
    
    % move selected to the top if we're in 2D view
    % Note: OpenGL rendering sorts by ZData, not by uistack
    twodee=get(handles.feature3_menu,'Value')==1;
    if(twodee)
        locs=find(selected_axes & feature_plots~=0);
        badlocs=find(~selected_axes & feature_plots~=0);
        for i = locs
            X=get(feature_plots(i),'XData');
            set(feature_plots(i),'ZData',ones(size(X))*eps);
        end
        for i = badlocs
            X=get(feature_plots(i),'XData');
            set(feature_plots(i),'ZData',zeros(size(X)));
        end
    end
    
    
    hist_plots=getappdata(handles.output,'hist_plots');
    
    if(any(hist_plots))
        valid_handles=ishandle(hist_plots) & hist_plots~=0;
        alpha(hist_plots(selected_axes & valid_handles),.9);
        alpha(hist_plots(~selected_axes & valid_handles),.6);
        if(twodee)
            locs=find(selected_axes & valid_handles);
            badlocs=find(~selected_axes & valid_handles);
            for i = locs
                X=get(hist_plots(i),'XData');
                set(hist_plots(i),'ZData',ones(size(X))*eps);
            end
            for i = badlocs
                X=get(hist_plots(i),'XData');
                set(hist_plots(i),'Zdata',zeros(size(X)));
            end
        end
    end
end


% -- Clears all the elements of a panel
function clear_panel(handle)
delete(get(handle,'Children'));
set(ancestor(handle,'figure'),'windowbuttonmotionfcn',[]);


% -- Fixes feature menu for plotting
function fix_feature_menu(handles)
% Determine which features we're going to calculate
num_channels=getappdata(handles.output,'num_channels');
channels=getappdata(handles.output,'channels');
num_PCs = get(handles.num_PCs_slider,'Max');
if(~iscell(num_PCs)), num_PCs={num_PCs}; end
num_PCs=max([num_PCs{:}]);
feature_list = cell(1,num_channels*num_PCs);
feature_names = cell(1,num_channels*num_PCs);

c=round(getappdata(handles.output,'channel_colors')*255);
ind=1;
for i = 1:num_channels
    for j = 1:num_PCs
        feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - PC %g</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i),j);
        feature_names{ind}=sprintf('%g - PC %g',channels(i),j);
        ind=ind+1;
    end
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Peak</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Peak',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Vall</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Valley',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Slope</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Slope',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Ener</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Energy',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Amp</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Amplitude',channels(i)); ind=ind+1;
end

% note: 'time' means spike time, 'hist' creates a 3D of the other features
% displayed
feature_list = [feature_list 'Time'];
feature_names = [feature_names 'Time'];

set(handles.feature1_menu,'String',feature_list);
set(handles.feature2_menu,'String',[feature_list(:);{'Hist'}]);
set(handles.feature3_menu,'String',[{'<none>'};feature_list(:);{'Hist'}]);
set(handles.feature1_menu,'Value',1);
if(num_channels>1),val=8; else val=2; end
set(handles.feature2_menu,'Value',val);
set(handles.feature3_menu,'Value',1);
setappdata(handles.output,'feature_names',feature_names);

% -- Plots mean waveforms in the right panel
function plot_mean_waveforms(handles)
global waveforms palette;
idx=getappdata(handles.output,'idx');
num_clusters=getappdata(handles.output,'num_clusters');
if(isempty(num_clusters) || num_clusters==0), return; end

t=getappdata(handles.output,'t');

set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Plotting waveforms...');

% Plot mean waveforms
include_zero=true; % tells get_mean_waveforms.m that idx starts at 0, not 1
mean_wfs = get_mean_waveforms(waveforms,idx,include_zero);
max_val = max(max(max(mean_wfs)));
min_val = min(min(min(mean_wfs)));
span=max_val-min_val; margin=span*0.05;
mtop=max_val+margin; mbot=min_val-margin;

dx=(t(end)-t(1))/4;
xgrid=t(1)+[dx 2*dx 3*dx];
dy=(mtop-mbot)/4;
ygrid=mbot+[dy 2*dy 3*dy];

channel_colors=getappdata(handles.output,'channel_colors');
cluster_colors = getappdata(handles.output,'cluster_colors');

% Determine location for axes. This is as wide as and up to the bottom of
% the waveform_control_panel and

panel_selected=getappdata(handles.output,'panel_selected');
features_pos=get(handles.features_panel,'position');
feature_width=features_pos(3);
pos=get(handles.feature_control_panel,'position');
total_width=(pos(3)/feature_width)*.99;
if(strcmp(panel_selected,'features'))
    total_height=pos(2)+.015;
elseif(strcmp(panel_selected,'cell_stats'))
    total_height=.985;
end

wf_axes=zeros(1,num_clusters);
frame_handles=zeros(1,num_clusters);
wf_buttons = zeros(1,num_clusters);

plot_width=total_width;
plot_height=total_height/num_clusters;

% plot waveforms
for i = 1:num_clusters
    wf_axes(i) = axes('Parent',handles.features_panel,'units','normalized'); %#ok<LAXES>
    set(wf_axes(i),'position',[0 total_height-i*plot_height plot_width plot_height]);
    hold all;
    if(i==1)
        set(wf_axes(i),'xtick',[],'ytick',[],'xcolor',palette.text,'ycolor',palette.text,'color',palette.background);
    else
        set(wf_axes(i),'xtick',[],'ytick',[],'xcolor',palette.text,'ycolor',palette.text,'color',palette.ax_background);
    end
    
    set(wf_axes(i),'colororder',channel_colors); hold(wf_axes(i),'all');
    plot(t,squeeze(mean_wfs(i,:,:))','linewidth',2,'linesmooth','on','Parent',wf_axes(i));
    axis(wf_axes(i),[t(1) t(end) mbot mtop]);
    
    minorgrid(xgrid,ygrid,wf_axes(i));
    
    set(wf_axes(i),'box','on','xticklabel',[],'yticklabel',[],'ButtonDownFcn',{@select_axes,handles});
    
    % determine square patch location
    a = axis(wf_axes(i));
    width=a(2)-a(1); margin_x = width*.01;
    height=a(4)-a(3); margin_y = height*.01;
    square_width = width * .2;
    square_height = height * .2;
    
    % add patch
    x=[a(2)-margin_x-square_width a(2)-margin_x-square_width a(2)-margin_x a(2)-margin_x];
    y=[a(3)+margin_y a(3)+margin_y+square_height a(3)+margin_y+square_height a(3)+margin_y];
    patch(x,y,cluster_colors(i,:));
    
    % place text immediately before patch
    text(a(2)-2*margin_x-square_width,a(3)+margin_y+square_height/2,sprintf('%g',length(find(idx==i-1))),...
        'fontsize',12,'fontweight','bold','fontname','Cambria','verticalalignment','middle',...
        'horizontalalignment','right','backgroundcolor',palette.ax_background,'margin',1,'color',cluster_colors(i,:),'Parent',wf_axes(i));
    
    % set frame
    a = axis(wf_axes(i));
    span_x=a(2)-a(1); span_y=a(4)-a(3);
    margin_x=span_x*.001; margin_y=span_y*.002;
    x=[a(1) a(1) a(2)-margin_x a(2)-margin_x a(1)];
    y=[a(3)+margin_y a(4) a(4) a(3)+margin_y a(3)+margin_y];
    frame_handles(i)=plot(x,y,'linewidth',10,'color',cluster_colors(i,:),'visible','off','Parent',wf_axes(i));
    
    % turn off hittest for all children
    c=get(wf_axes(i),'Children');
    set(c,'hittest','off');
end
selected_axes=getappdata(handles.output,'selected_axes');
set(frame_handles(selected_axes),'visible','on');
setappdata(handles.output,'wf_axes',wf_axes);
setappdata(handles.output,'frame_handles',frame_handles);
setappdata(handles.output,'wf_buttons',wf_buttons);

% Reset labels
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);
drawnow;

% --- Executes on button press in recluster_button.
function recluster_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to recluster_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% --- Determines which features have been selected by the user
function feature_mask = get_feature_mask(handles)
% Determine which features to cluster along
PCs=get(handles.num_PCs_slider,'Value');
if(~iscell(PCs)),PCs={PCs}; end
PCs=cell2mat(PCs);

max_PCs=get(handles.num_PCs_slider,'Max');
if(~iscell(max_PCs)),max_PCs={max_PCs}; end
max_PCs=max([max_PCs{:}]);
num_channels=getappdata(handles.output,'num_channels');

% slope valley energy peak amplitude
feat_len=max_PCs+5;

feature_mask=false(feat_len,num_channels);

peaks=get(handles.peak_checkbox,'Value'); if(iscell(peaks)),peaks=cell2mat(peaks); end
valleys=get(handles.valley_checkbox,'Value'); if(iscell(valleys)),valleys=cell2mat(valleys); end
slopes=get(handles.slope_checkbox,'Value'); if(iscell(slopes)),slopes=cell2mat(slopes); end
energies=get(handles.energy_checkbox,'Value'); if(iscell(energies)),energies=cell2mat(energies); end
amplitudes=get(handles.amplitude_checkbox,'Value'); if(iscell(amplitudes)),amplitudes=cell2mat(amplitudes); end

for i=1:num_channels
    ind=1;
    num_PCs=PCs(i);
    feature_mask(1:num_PCs,i)=true; ind=ind+max_PCs;
    if(peaks(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(valleys(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(slopes(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(energies(i)), feature_mask(ind,i)=true; end
    ind=ind+1;
    if(amplitudes(i)), feature_mask(ind)=true; end
end

selected_channels=getappdata(handles.output,'selected_channels');
if(~isempty(selected_channels) && any(selected_channels))
    feature_mask(:,~selected_channels)=false;
end
feature_mask=reshape(feature_mask,1,[]);

% -- determine axis selection
function select_axes(src,~,handles)
% get waveform axis handles to compare
wf_axes = getappdata(handles.output,'wf_axes');

% get list of currently selected axes
selected_axes = getappdata(handles.output,'selected_axes');
%selected_axes=false(1,length(selected_axes));

% determine which axis was selected
axis_selected = find(wf_axes==src);

% get selection type
type=get(handles.output,'SelectionType');
% reset others to zero
if(strcmpi(type,'normal'))
    % remove frame from all selected axes, select/unselect current
    toggle_frame(unique([find(selected_axes) axis_selected]),handles);
else
    toggle_frame(axis_selected,handles);
end

selected_axes = getappdata(handles.output,'selected_axes');
axis_selected = find(selected_axes);
if(~isempty(axis_selected))
    setappdata(handles.output,'cluster_number',axis_selected-1);
end

panel_selected=getappdata(handles.output,'panel_selected');

if(strcmp(panel_selected,'features'))
    % update feature plot to put this plot on top
    % move current waveform view to the top
    wf_plots=getappdata(handles.output,'wf_plots');
    for i = 1:size(wf_plots,1)
        uistack(wf_plots(i,axis_selected),'top');
    end

    fix_featurepanel_display(handles);
elseif(strcmp(panel_selected,'cell_stats'))
    update=true;
    display_cell_stats(handles,update);
end

% Updates status label
idx=getappdata(handles.output,'idx');
total_wfs=false(1,length(idx));
for i = 1:length(axis_selected)
    total_wfs=total_wfs | idx==axis_selected(i)-1;
end
total_wfs=sum(total_wfs);
set(handles.status_label,'String',sprintf('%g waveforms.',total_wfs));

% Draw or remove a frame around a waveform plot
% note: if ax==-1, then turn all off
% also note: cluster_number = selected_axis-1
function selected_axes = toggle_frame(ax, handles)
frame_handles=getappdata(handles.output,'frame_handles');
selected_axes=getappdata(handles.output,'selected_axes');

% instruction to turn all off
if(ax==-1)
    if(~isempty(selected_axes))
        selected_axes(:)=0;
        setappdata(handles.output,'selected_axes',selected_axes);
    end
    return;
end

% toggle frames
for i = 1:length(ax)
    if(selected_axes(ax(i)))
        selected_axes(ax(i))=false;
        set(frame_handles(ax(i)),'visible','off','hittest','off');
    else
        selected_axes(ax(i))=true;
        set(frame_handles(ax(i)),'visible','on','hittest','off');
    end
end
drawnow;
% update merge & split buttons
num_selected=length(find(selected_axes));
if(num_selected>1)
    set(handles.merge_button,'enable','on');
else
    set(handles.merge_button,'enable','off');
end
if(num_selected==0), set(handles.split_button,'enable','off');
else set(handles.split_button,'enable','on');
end
setappdata(handles.output,'selected_axes',selected_axes);

% --- Executes on button press in merge_button.
function merge_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to merge_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Merging...'); drawnow; pause(0.01);

selected_axes=getappdata(handles.output,'selected_axes');
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);
cluster_nums=find(selected_axes)-1;
s=sprintf('Clusters %g, ',cluster_nums(1));

for i = 2:length(cluster_nums)
    idx(idx==cluster_nums(i))=cluster_nums(1);
    s=sprintf('%s%g',s,cluster_nums(i));
end

% fix cluster numbers
setappdata(handles.output,'idx',idx);
fix_clusters(handles);

s=sprintf('%s merged into cluster %g.',s,cluster_nums(1));
setappdata(handles.output,'cluster_number',cluster_nums(1));

panel_selected=getappdata(handles.output,'panel_selected');
switch panel_selected
    case 'features'
        show_featurespanel(handles);
    case 'cell_stats'
        stats_button_Callback([],[],handles);
end

set(handles.status_label,'String',s);
set(handles.status_light,'BackgroundColor',palette.green);

% --- Executes on button press in split_button.
function split_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to split_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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


% -- plays a beep or notifies android
function stop_alert(phone_alert)
if(exist('phone_alert','var') && phone_alert==true)
    notify_android('Msort finished');
end
soundbeep;


% -- updates cluster numbers so there are no clusters with zero waveforms
function fix_clusters(handles)
global waveforms;
idx=getappdata(handles.output,'idx');
if(isempty(waveforms))
    idx=[]; setappdata(handles.output,'idx',idx);
    setappdata(handles.output,'timestamps',[]); return;
end
if(isempty(idx))
    idx=zeros(1,size(waveforms,2),'uint8');
    num_clusters=1;
    setappdata(handles.output,'cluster_number',0);
else
    u=unique(idx); u(u==0)=[];
    % find missing clusters
    u(u==0)=[];
    for i = 1:length(u)
        % renumber this cluster
        if(u(i)>i)
            idx(idx==u(i))=i;
        end
    end
    u=unique([0 u]); num_clusters=length(u);
end
setappdata(handles.output,'num_clusters',num_clusters);
setappdata(handles.output,'idx',idx);

% add noise cluster if it's missing
cluster_colors=get_cluster_colors(handles);
setappdata(handles.output,'cluster_colors',cluster_colors);
cluster_number=getappdata(handles.output,'cluster_number');
if(cluster_number>num_clusters-1)
    setappdata(handles.output,'cluster_number',0);
end
setappdata(handles.output,'selected_axes',false(1,num_clusters));


% --- Executes on button press in voltage_plot_button.
function voltage_plot_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to voltage_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clear_panel(handles.features_panel);
%hide_panel(handles.features_panel);
global palette;
plot_voltages(handles);
setappdata(handles.output,'panel_selected','voltage');
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);


% --- Executes on button press in export_button.
% Note: this is the "save" button
function export_button_Callback(~, ~, handles)
% hObject    handle to export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waveforms palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Saving .mat file...'); drawnow;
mat_file=getappdata(handles.output,'mat_file');
group_num=getappdata(handles.output,'group_num');
if(isempty(mat_file))
    basename=getappdata(handles.output,'basename');
    pathname=getappdata(handles.output,'pathname');
    mat_file=sprintf('%s%s.mat',pathname,basename);
end

% save additional data to mat file - rename to group-specific variables
varout_names={...
    sprintf('original_timestamps_%g',group_num), ...
    sprintf('original_waveforms_%g',group_num),...
    sprintf('timestamps_%g',group_num)...
    sprintf('t_%g',group_num), ...
    sprintf('idx_%g',group_num), ...
    sprintf('trigger_ch_%g',group_num),...
    sprintf('num_PCs_%g',group_num), ...
    sprintf('normalized_waves_%g',group_num), ...
    sprintf('t_before_%g',group_num),...
    sprintf('t_after_%g',group_num),...
    sprintf('num_clusters_%g',group_num)};
varin_names={'original_timestamps',...
    'original_waveforms','timestamps',...
    'original_t',...
    'idx',...
    'trigger_ch',...
    'num_PCs',...
    'normalized_waves',...
    't_before',...
    't_after',...
    'num_clusters'};

% Save waveform amplitudes
t=getappdata(handles.output,'t');
if(~isempty(t))
    zeroloc=find_index(t,0);
    amplitudes=squeeze(waveforms(:,:,zeroloc)); %#ok<NASGU>
    eval(sprintf('amplitudes_%g = amplitudes;',group_num));
    save(mat_file,'-append','-v7.3',sprintf('amplitudes_%g',group_num));
end

% save variables one by one
set(handles.status_label,'String','Saving data to .mat file...'); drawnow;
for i = 1:length(varin_names)
    var=getappdata(handles.output,varin_names{i}); %#ok<NASGU>
    %if(~isempty(var))
        eval(sprintf('%s = var;',varout_names{i}));
        clear var;
        save(mat_file,'-append','-v7.3',varout_names{i});
        eval(sprintf('clear %s;',varout_names{i}));
    %end
end


% export struct containing new thresholds
exp=getappdata(handles.output,'exp');
xml_file=getappdata(handles.output,'xml_file');
struct2xml(exp,xml_file);
setappdata(handles.output,'exp',exp);

% update status label
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);

% --- Executes on button press in param_button.
function param_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to param_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calculate_features(handles);


% --- Executes on button press in stats_button.
% This function displays ISI and other statistics
function stats_button_Callback(~, ~, handles)
% hObject    handle to stats_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.output,'panel_selected','cell_stats');
set(handles.feature_control_panel,'visible','off');
clear_panel(handles.features_panel); drawnow;
display_cell_stats(handles);

% -- Displays cell statistics.  For now: ISI and rate function
function display_cell_stats(handles,update)
global Ts palette;

if(~exist('update','var') || isempty(update)), update=false; end

idx=getappdata(handles.output,'idx');

% do nothing if we haven't clustered yet
if(isempty(idx) || ~any(idx))
    set(handles.status_light,'BackgroundColor',[255 178 0]/255);
    set(handles.status_label,'String','No clusters yet!');
    return;
end

selected_axes=getappdata(handles.output,'selected_axes');
selected_axes(1)=[]; % remove noise axis
if(~any(selected_axes)),selected_axes=~selected_axes; end
% note: cluster 0 is noise, ignore it
u=unique(idx); u(u==0)=[];

margin_bot = 0.05; margin_top = .05;
margin_left = 0.18; margin_right = .01;
margin_middle_x=.02; margin_middle_y = .03;
plot_width=(1-2*margin_middle_x-margin_left-margin_right)/3;
plot_height = (1-margin_bot-margin_top - (length(u)-1)*margin_middle_y)/length(u);

%---plotting parameters
% ISI
%dt=1e-2; edges=0:dt:.5;
dt=1e-3; edges=0:dt:.1;
ISI_xtick=linspace(0,edges(end),5);
ISI_xticklabel=num2str(ISI_xtick');

% Rate Function
%span=Ts(end)-Ts(1);
resolution = 30; % in seconds
t_rate=unique(round(Ts(1)):resolution:round(Ts(end)));
xtick_ind=unique(round(linspace(1,length(t_rate),6)));
rate_xtick=t_rate(xtick_ind);
%rate_xtick=[round(Ts(1)) round(Ts(1)+span/4) round(Ts(1)+span/2) round(Ts(1)+3*span/4) round(Ts(end))];
rate_xticklabel=num2str(rate_xtick');
%locs=rate_xtick<Ts(1) | rate_xtick>Ts(end);
%rate_xtick(locs)=[];
%rate_xticklabel(locs)=[];



% Autocorrelation
lag=2e-2; % note: this is the total timespan of the autocorrelation function
num_bins=200; bin_width=2*lag/num_bins;
xspace=lag/2; autocorr_xtick=-lag:xspace:lag;
autocorr_xticklabel=num2str(autocorr_xtick');

% only make certain plots visible and leave; otherwise plot everything
if(update)
    ax=getappdata(handles.output,'ax_cell_stats');
else
    set(handles.status_light,'BackgroundColor',palette.red);
    set(handles.status_label,'String','Generating Cell Statistics...');

    timestamps=getappdata(handles.output,'timestamps');

    % grab spike trains
    num_clusters=length(u);
    spikes=cell(1,num_clusters);
    for i = 1:num_clusters
        spikes{i} = timestamps(idx==i);
    end

    % display plot
    plot_mean_waveforms(handles);

    ax=zeros(num_clusters,3);
    colors = getappdata(handles.output,'cluster_colors');
    if(size(colors,1)>num_clusters), colors(1,:)=[]; end
    for i = 1:num_clusters
        y = 1 - margin_top -(i-1)*margin_middle_y-i*plot_height;

        % ISI plot
        x = margin_left;
        ISI = diff(spikes{i});
        if(length(ISI)<2), n=0;
        else n=histc(ISI,edges); n=n./sum(n); n(isnan(n))=0; end
        ax(i,1) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        b=bar(edges,n,'histc'); set(b,'facecolor',colors(i,:),'edgecolor','none');
        min_n=0; max_n=max(n)*1.05;
        if(max_n==min_n), min_n=0; max_n=1; end
        axis([edges(1) edges(end) min_n max_n]);
        y_tick=linspace(min_n,max_n,5);
        set(ax(i,1),'ytick',y_tick,'yticklabel',[]);
        grid on;

        % put up text
        mean_rate = 1/mean(ISI);
        margin_x=(edges(end)-edges(1))*.01; margin_y=(max_n-min_n)*.05;
        text(edges(end)-margin_x,max_n-margin_y,sprintf('%2.1f spikes/s',mean_rate),'fontsize',10,...
            'fontweight','bold','color',colors(i,:),'horizontalalignment','right','verticalalignment','top','backgroundcolor',palette.ax_background);

        % rate plot
        x = margin_left + plot_width + margin_middle_x;
        %rate=smooth(histc(spikes{i},t),max(round(Ts(end)/100),1));
        rate = histc(spikes{i},t_rate)/resolution;
        ax(i,2) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        plot(t_rate,rate,'linewidth',2,'color',colors(i,:),'linesmooth','on');
        axis([round(Ts(1)) round(Ts(end)) 0 max(rate)*1.05]);
        ytick=(0:3)*max(rate)*1.05/4; ytick=unique(round(ytick*10)/10);
        set(ax(i,2),'ytick',ytick);
        grid on;

        % Autocorrelation plot
        s=spikes{i};
        x = margin_left + 2*plot_width+2*margin_middle_x;
        ax(i,3) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        
        % only use up to the some # of spikes spikes
        % note that we can't randomly sample
        MAX_SPIKES=3e3;
        if(length(s)>MAX_SPIKES), s=s(1:MAX_SPIKES); end
        [t,n]=crosscorr(s,s,bin_width,lag);
        %loc=find(t==0);
        % remove center ripple
        %n([loc-1 loc loc+1]) = (n(loc-2)+n(loc+2))/2;
        n(t==0)=0;
        tt=linspace(t(1),t(end),1000);
        nn=interp1(t,n,tt,'pchip');
        plot(ax(i,3),tt,nn,'linesmooth','on','color',colors(i,:),'linewidth',2);
        grid on;
        maxy=max(nn)*1.05;
        ytick=0:maxy/4:maxy; ytick=round(ytick*10)/10;
        set(ax(i,3),'ytick',ytick);
        if(max(nn))==0,nn=1; end
        axis(ax(i,3),[-lag lag 0 max(nn)*1.05]);
    end
end
set(ax,'color',palette.ax_background,'fontweight','bold','box','on','linewidth',1.5,'xcolor',palette.text,'ycolor',palette.text);
c_off=get(ax(~selected_axes,:),'children');
if(iscell(c_off)), c_off=cell2mat(c_off); end
off=[reshape(ax(~selected_axes,:),1,[]) c_off(:)'];
c_on=get(ax(selected_axes,:),'children');
if(iscell(c_on)), c_on=cell2mat(c_on); end
on=[reshape(ax(selected_axes,:),1,[]) c_on(:)'];
set(off,'visible','off');
set(on,'visible','on');

% Deterine axis positions
num_selected=sum(selected_axes);
width=plot_width;
height = (1-margin_bot-margin_top - (num_selected-1)*margin_middle_y)/num_selected;
x1=margin_left;
x2=margin_left+plot_width+margin_middle_x;
x3=margin_left+2*(plot_width+margin_middle_x);
y=(0:num_selected-1)*(height+margin_middle_y)+margin_bot;
y=y(end:-1:1);
pos1=[ones(num_selected,1)*x1 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos2=[ones(num_selected,1)*x2 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos3=[ones(num_selected,1)*x3 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos=num2cell([pos1;pos2;pos3],2);
ax_old=ax;
ax=ax(selected_axes,:);
ax_selected=reshape(ax,1,[]);
set(ax_selected,{'position'},pos,'xcolor',palette.text,'ycolor',palette.text);

% generate text labels
title(ax(1,1),'ISI','fontsize',14,'fontweight','bold','color',palette.text);
title(ax(1,2),'Firing Rate','fontsize',14,'fontweight','bold','color',palette.text);
title(ax(1,3),'Autocorrelatoin','fontsize',14,'fontweight','bold');
xlabel(ax(end,1),'ISI (s)','fontsize',14,'fontweight','bold','color',palette.text);
xlabel(ax(end,2),'Time (s)','fontsize',14,'fontweight','bold','color',palette.text);
xlabel(ax(end,1),'Time (s)','fontsize',14,'fontweight','bold','color',palette.text);
set(ax(:,1),'xtick',ISI_xtick,'xticklabel',ISI_xticklabel);
set(ax(:,2),'xtick',rate_xtick,'xticklabel',rate_xticklabel);
set(ax(:,3),'xtick',autocorr_xtick,'xticklabel',autocorr_xticklabel);
set(ax(1:end-1,:),'xticklabel',[]);
set(ax,'xcolor',palette.text,'ycolor',palette.text);

% generate text labels
for i = 1:length(ax_selected)
    title(ax_selected(i),[]);
    xlabel(ax_selected(i),[]);
end

if(~update)
    % update selected panel
    setappdata(handles.output,'ax_cell_stats',ax_old);
    set(handles.status_label,'String','Done!');
    set(handles.status_light,'BackgroundColor',palette.green);
end

% --- displays waveforms in each of the waveform axis plots
% note: only selected channels are displayed, and channels are displayed in
% separate axes
function show_all_waveforms(handles)
global waveforms palette;
MAX_WFS=str2double(get(handles.max_wfs_field,'String')); % max # waveforms to display
gray=[.5 .5 .5];

% create empty wf_plots matrix
channels=getappdata(handles.output,'channels');
num_channels=length(channels);

% clear any plots already around
wf_plots=getappdata(handles.output,'wf_plots');
if(~isempty(wf_plots)), delete(wf_plots(ishandle(wf_plots) & wf_plots~=0)); end
errorbar_plots=getappdata(handles.output,'errorbar_plots');
if(~isempty(errorbar_plots)), delete(errorbar_plots(ishandle(errorbar_plots) & errorbar_plots~=0)); end
gridline_plots=getappdata(handles.output,'gridline_plots');
if(~isempty(gridline_plots)), delete(gridline_plots(ishandle(gridline_plots) & gridline_plots~=0)); end

% clear current plots if they exist
ax_waveforms=getappdata(handles.output,'ax_waveforms');
delete(cell2mat(get(ax_waveforms,'children')));

% determine which channels to plot, and which groups to plot
cluster_numbers=getappdata(handles.output,'cluster_number');
if(isempty(cluster_numbers)), cluster_numbers=0; end
cluster_colors=getappdata(handles.output,'cluster_colors');
channels=getappdata(handles.output,'channels');
num_clusters=getappdata(handles.output,'num_clusters');

%waveforms=getappdata(handles.output,'waveforms');
idx=getappdata(handles.output,'idx');
t=getappdata(handles.output,'t')*1e6;

% determine ticks
supermaxval=max(waveforms(:));
superminval=min(waveforms(:));
if(get(handles.fit_wf_checkbox,'Value'))
    locs=false(1,size(waveforms,2));
    for i = 1:length(cluster_numbers)
        locs=locs | idx==cluster_numbers(i);
    end
    wfs=waveforms(:,locs,:);
    maxval=max(wfs(:));
    minval=min(wfs(:));
else
    maxval=supermaxval;
    minval=superminval;
end

if(isempty(minval)), minval=0; maxval=1; end

xtick=[-200 -100 0 100 200 300 400]; % not sure if this should be hardcoded
ytick=sort(floor(superminval*1.05)*2:2:ceil(supermaxval*1.05)*2);

set(ax_waveforms,'visible','off','xlim',[t(1) t(end)],...
    'xcolor',palette.text,'ycolor',palette.text,'xtick',xtick,'ytick',ytick,'box','on','linewidth',1.5);

display_tools=getappdata(handles.output,'display_cleanup_tools');
t_err=double(reshape([t; t; nan(1,length(t))],1,[]));

include_zero=true;
mean_wfs=get_mean_waveforms(waveforms,idx,include_zero);

% generate plots
% generate new plot handles
wf_plots=zeros(num_channels,num_clusters);
errorbar_plots=zeros(num_channels,num_clusters);
% each gridline is its own plot. Might update this in future to make single
% line since it's faster and easier to handle.
gridline_plots=zeros(1,num_channels);
ch_title_plots=zeros(1,num_channels);
max_wfs=zeros(length(channels),num_clusters);
min_wfs=zeros(length(channels),num_clusters);
for i = 1:length(channels)
    for j = 0:num_clusters-1
        wfs=squeeze(waveforms(i,idx==j,:))';
        if(isvector(wfs)),wfs=wfs'; end
        num_wfs=size(wfs,2);
        if(num_wfs>0)
            max_wfs(i,j+1)=max(wfs(:));
            min_wfs(i,j+1)=min(wfs(:));
        end
        
        % Determine std deviation error bars (display cleanup tools)
        std_dev_pts = std(wfs,[],2);
        mean_wf_pts = squeeze(mean_wfs(j+1,i,:));
        errorbar_value=get(handles.errorbar_slider,'Value');
        std_dev_pts = [mean_wf_pts-errorbar_value*std_dev_pts mean_wf_pts+errorbar_value*std_dev_pts]';
        
        
        if(num_wfs)>0
            % only plot up to MAX_WFS waveforms
            if(num_wfs>MAX_WFS)
                ind=round(linspace(1,num_wfs,MAX_WFS));
                wfs=wfs(:,ind,:);
                num_wfs=MAX_WFS;
            end
            
            wfs=[wfs; nan(1,num_wfs)]; %#ok<AGROW>
            t_wfs=repmat(t',1,num_wfs);
            t_wfs=[t_wfs; nan(1,num_wfs)]; %#ok<AGROW>
            if(any(j==cluster_numbers))
                wf_plots(i,j+1)=plot(ax_waveforms(i),...
                    t_wfs(:),wfs(:),'color',cluster_colors(j+1,:),'linesmooth','on','visible','on','hittest','off');
            else
                wf_plots(i,j+1)=plot(ax_waveforms(i),...
                    t_wfs(:),wfs(:),'color',cluster_colors(j+1,:),'linesmooth','on','visible','off','hittest','off');
            end
        end
        gridline_plots(i)=minorgrid(xtick,ytick,ax_waveforms(i),false,[t(1) t(end) superminval*1.2 supermaxval*1.2]);
        
        % plot std devs
        std_dev_pts=reshape([std_dev_pts; nan(1,length(t))],1,[]);
        errorbar_plots(i,j+1)=plot(ax_waveforms(i),t_err,std_dev_pts,'color',palette.red,'visible','off','hittest','off');
        if(display_tools && j==cluster_numbers(1))
            set(errorbar_plots(i),'visible','on');
        end
    end
    
    % put channel tag
    x=double(t(1)+(t(end)-t(1))*.01);
    y=double(maxval*1.05*.99);
    
    % add channel tag
    ch_title_plots(i)=text(x,y,sprintf('Ch %g',channels(i)),'fontsize',20,'backgroundcolor',palette.ax_background,...
        'parent',ax_waveforms(i),'color',palette.text,'horizontalalignment','left',...
        'verticalalignment','top','hittest','off');
end

setappdata(handles.output,'min_wfs',min_wfs);
setappdata(handles.output,'max_wfs',max_wfs);
setappdata(handles.output,'wf_plots',wf_plots);
setappdata(handles.output,'errorbar_plots',errorbar_plots);
setappdata(handles.output,'gridline_plots',gridline_plots);
setappdata(handles.output,'ch_title_plots',ch_title_plots);
%set([gridline_plots(:)' wf_plots(:)' errorbar_plots(:)' ch_title_plots(:)'],'hittest','off');
children=get(ax_waveforms,'children');
%if(~iscell(children)), children={children}; end
if(iscell(children)), children=cell2mat(children); end
set(children,'hittest','off');


% -- Generates merge buttons in the small panel; size depends no number of clusters
function generate_merge_buttons(handles)
num_clusters=getappdata(handles.output,'num_clusters');
merge_buttons=getappdata(handles.output,'merge_buttons');
if(~isempty(merge_buttons) && any(ishandle(merge_buttons)))
    delete(merge_buttons);
end
margin_x=.03; margin_y=.03;
margin_top=.05; margin_bot=.05;
margin_left=.05; margin_right=.05;
num_rows=floor(sqrt(num_clusters));
num_cols=ceil(num_clusters/num_rows);

button_height=(1-margin_top-margin_bot-(num_rows-1)*margin_y)/num_rows;
button_width=(1-margin_left-margin_right-(num_cols-1)*margin_x)/num_cols;

buttons=zeros(1,num_clusters);
cluster_colors=getappdata(handles.output,'cluster_colors');

ind=1;
for i = 1:num_rows
    y=1-margin_top-(i-1)*margin_y-i*button_height;
    for j = 1:num_cols
        if(ind>num_clusters), ind=ind+1; continue; end
        x=margin_left+(j-1)*(margin_x+button_width);
        buttons(ind)=uicontrol('Style','Pushbutton','units','normalized','String',[],'BackgroundColor',cluster_colors(ind,:)/3,...
            'Parent',handles.merge_panel,'selectionhighlight','off','enable','off','callback',{@move_wfs,handles,ind-1},...
            'Position',[x y button_width button_height]);
        ind=ind+1;
    end
end
setappdata(handles.output,'merge_buttons',buttons);


% -- creates a new cluster, given an index
function create_new_cluster(handles)
% Retrieve currently selectected cluster
cluster_number=getappdata(handles.output,'cluster_number');
if(isempty(cluster_number)),cluster_number=0; end
setappdata(handles.output,'cluster_number',cluster_number);

% Remove waveforms from old cluster and create new waveform
idx=getappdata(handles.output,'idx');
selected_waveforms=getappdata(handles.output,'selected_waveforms');
if(isempty(selected_waveforms))
    selected_waveforms=idx==cluster_number;
end
num_clusters=getappdata(handles.output,'num_clusters');
idx(selected_waveforms)=num_clusters;
setappdata(handles.output,'num_clusters',num_clusters+1);
setappdata(handles.output,'idx',idx);
setappdata(handles.output,'idx_old',idx);
setappdata(handles.output,'selected_waveforms',[]);
setappdata(handles.output,'cluster_number',num_clusters);

fix_clusters(handles);
show_featurespanel(handles);


% -- deselects all waveforms
function deselect_wfs(handles) %#ok<DEFNU>
global wf_plots dim_channel_color main_channel;

buttons=getappdata(handles.output,'merge_buttons');
%selected_wfs=getappdata(handles.features_panel,'selected_wfs');
set(wf_plots(selected_wfs),'color',dim_channel_color(main_channel,:));
setappdata(handles.features_panel,'selected_wfs',[]);
% re-disable buttons
for i = 1:length(buttons)
    if(strcmp(get(buttons(i),'enable'),'on'))
        set(buttons(i),'backgroundcolor',get(buttons(i),'backgroundcolor')/3,'enable','off');
    end
end

% -- Select particular waveforms with the mouse
% Note: there's some error that occurs. After merging two clusters, you
% can't reselect waveforms unless you move away from Feature view and move
% back. This should be looked into.
function select_wfs(hObject,~,handles)
global ax_coord x_ratio y_ratio ax_lim locsX locsY h x y palette;
fig=handles.output;

ax=hObject;

set([ax fig 0],'units','pixels');

ax_pos=get(hObject,'position');
ax_lim=axis(hObject);
fig_pos=get(fig,'position');

% get axis coordinates
g=hObject;
x_offset=fig_pos(1); y_offset=fig_pos(2);
while(g~=fig)
    tmp=get(g,'units'); set(g,'units','pixels');
    pos=get(g,'position');
    x_offset=x_offset+pos(1);
    y_offset=y_offset+pos(2);
    set(g,'units',tmp);
    g=get(g,'parent');
end
ax_coord = [x_offset y_offset x_offset y_offset] + [0 0 ax_pos(3) ax_pos(4)];

% determine mouse position conversion
x_ratio=(ax_lim(2)-ax_lim(1))/ax_pos(3);
y_ratio=(ax_lim(4)-ax_lim(3))/ax_pos(4);

% set up initial plot
mouse_pt = get(0,'PointerLocation');
x1=(mouse_pt(1)-ax_coord(1))*x_ratio+ax_lim(1);
y1=(mouse_pt(2)-ax_coord(2))*y_ratio+ax_lim(3);
x=x1; y=y1;
locsX=[x1 x];
locsY=[y1 y];
hold on; axis manual;
h=plot(hObject,locsX,locsY,'linewidth',1.5,'color',palette.text);
set(h,'XDataSource','locsX','YDataSource','locsY');

set(fig,'windowButtonMotionFcn',@connect_line,'windowButtonUpFcn',{@complete_line,ax,get_cluster_colors(handles)});

% --- Used while user has already clicked and is moving mouse
function connect_line(~,~)
global ax_coord ax_lim x_ratio y_ratio locsX locsY h;

mouse_pt = get(0,'PointerLocation');

% calculate new point
x=(mouse_pt(1)-ax_coord(1))*x_ratio+ax_lim(1);
y=(mouse_pt(2)-ax_coord(2))*y_ratio+ax_lim(3);

locsX(2)=x;
locsY(2)=y;
set(h,'XData',locsX,'YData',locsY);


% --- When user releases mouse, determines where the line was and selects
% all waveforms. It then draws "highlighted" waveforms overtop to indicate
% which waveforms were highlighted. Note: hidden waveforms may appear

% Note: hObject is handles.output
function complete_line(hObject,~,ax,cluster_colors)
global waveforms;
global locsX locsY h ctrl;

if(~ctrl)
    handles.output=hObject;
    deselect_button_Callback([],[],handles);
end

set(getappdata(hObject,'ax_waveforms'),...
    'units','normalized');
set(hObject,'Units','Normalized','windowButtonMotionFcn',[],'windowButtonUpFcn',[]);

delete(h);

% grab waveforms
t=getappdata(hObject,'t')*1e6;
cluster_numbers=getappdata(hObject,'cluster_number');

% determine which axis was selected
ax_waveforms=getappdata(hObject,'ax_waveforms');
channel = ax_waveforms==ax;

% get time range
t1=find(t<min(locsX),1,'last'); t2=find(t>max(locsX),1,'first');
if(isempty(t1)),t1=1; end
if(isempty(t2)),t2=length(t); end
t_range=t1:t2;

% get only waveforms in these clusters
idx=getappdata(hObject,'idx');
idx=repmat(idx,length(cluster_numbers),1);
cluster_numbers=repmat(cluster_numbers',1,size(idx,2));
wf_locs=any(idx==cluster_numbers,1);

wfs=squeeze(waveforms(channel,wf_locs,t_range));
if(length(t_range)>1 && isvector(wfs)),wfs=wfs(:)'; end
num_wfs=size(wfs,1);
t_new=t(t_range); t_length=length(t_new);

% prepare line segments for intersection test
XY2 = [locsX(1) locsY(1) locsX(2) locsY(2)];

wfs2=cat(3,wfs(:,1:end-1),wfs(:,2:end));
wfs2=permute(wfs2,[3 2 1]);
wfs2=reshape(wfs2,2,[])';
t2=repmat([t_new(1:end-1)' t_new(2:end)'],num_wfs,1);

XY1=[t2(:,1) wfs2(:,1) t2(:,2) wfs2(:,2)];
intersection=lineSegmentIntersect(XY1,XY2);
intersection=reshape(intersection,t_length-1,[]);
intersection=any(intersection,1);

if(~any(intersection))
    return;
end
% get actual waveform locations
locs=find(wf_locs);
locs=locs(intersection);

if(ctrl)
    old_locs=getappdata(hObject,'selected_waveforms');
else
    old_locs=false(1,length(idx));
end
if(length(old_locs)~=size(waveforms,2)),old_locs=false(1,length(idx)); end
old_locs(locs)=true; locs=old_locs;

% draw new selected waveforms
selected_waveforms=waveforms(:,locs,:);

selected_channels=getappdata(hObject,'selected_channels');
if(~any(selected_channels)),selected_channels=~selected_channels; end

% remove old plots
selected_wf_plots=getappdata(hObject,'selected_wf_plots');
delete(selected_wf_plots(ishandle(selected_wf_plots) & selected_wf_plots~=0));

% scroll through each waveform axis, add selected waveform lines
for i=1:length(ax_waveforms);
    wfs=selected_waveforms(i,:,:);
    num_wfs=size(wfs,2);
    wfs=squeeze(wfs);
    if(num_wfs==1),wfs=wfs'; end
    wfs=[wfs';nan(1,num_wfs)];
    wfs=reshape(wfs,1,[]);
    
    t_tmp=repmat([t nan],1,num_wfs);
    
    if(selected_channels(i))
        selected_wf_plots(i)=plot(ax_waveforms(i),t_tmp,wfs,'linesmooth','on','color','k','hittest','off','visible','on');
    else
        selected_wf_plots(i)=plot(ax_waveforms(i),t_tmp,wfs,'linesmooth','on','color','k','hittest','off','visible','off');
    end
end
drawnow;

setappdata(hObject,'selected_wf_plots',selected_wf_plots);
setappdata(hObject,'selected_waveforms',locs);

cluster_numbers=getappdata(hObject,'cluster_number');
num_clusters=getappdata(hObject,'num_clusters');
locs=1:num_clusters; locs(cluster_numbers+1)=[];
if(~isempty(locs))
    merge_buttons=getappdata(hObject,'merge_buttons');
    merge_buttons=merge_buttons(locs);
    set(merge_buttons,'enable','on');
    cluster_colors=cluster_colors(locs,:);
    for i = 1:length(merge_buttons)
        set(merge_buttons(i),'BackgroundColor',cluster_colors(i,:));
    end
end

% --- Intercept keyboard commands
% note: also called by edit menu -> screenshot
function key_press(~,eventData,handles)
global ctrl window_size palette;

% toggle ctrl key on and off
if(strcmp(eventData.Key,'control'))
    ctrl=true;

% toggle hide non-selected clusters from features panel
elseif(strcmp(eventData.Key,'h'))
    clusters_hidden=getappdata(handles.output,'clusters_hidden');
    if(clusters_hidden)
        setappdata(handles.output,'clusters_hidden',false);
        fix_featurepanel_display(handles)
        set(handles.hide_clusters_label,'visible','on');
    else
        setappdata(handles.output,'clusters_hidden',true);
        fix_featurepanel_display(handles)
        set(handles.hide_clusters_label,'visible','off');
    end
    
% This is the print screen command
elseif(strcmp(eventData.Key,'p'))
	% determine what's on the screen
    panel_selected=getappdata(handles.output,'panel_selected');
    
    % create blank figure
    scrsz = get(0,'ScreenSize');
    width = scrsz(3)*.9; height = scrsz(4)*.9;
    left=(scrsz(3)-width)*.5; bot=(scrsz(4)-height)*.5;
    f=figure('outerposition',[left bot width height],'color',palette.background,'PaperPositionMode','Auto','InvertHardCopy','off','visible','off');
    switch panel_selected
        % voltage trace
        case 'voltage'
            voltage_plot=getappdata(handles.output,'voltage_plot');
            % need to adjust voltage plot a little after we copy
            c=copyobj(voltage_plot,f);
            pos=get(c,'position');
            pos(1)=.05; pos(4)=.94; set(c,'position',pos);
            
            ch_voltage_tags=getappdata(handles.output,'ch_voltage_tags');
            c=copyobj(ch_voltage_tags,f);
            pos=get(c,'position');
            if(iscell(pos)),pos=cell2mat(pos); end
            pos(:,3)=.02;
            set(c,{'position'},num2cell(pos,2));
            
            % if black and white is checked, convert all axis children to black,
            % make plot white
            if(get(handles.screenshot_bw_checkbox,'Value'))
                ch=get(f,'children');
                control_locs=findobj(f,'type','uicontrol');
                ax_locs = findobj(f,'type','axes');
                ch2=get(ax_locs,'children');

                % set colors
                set(f,'color',[1 1 1]);
                set(ax_locs,'color',[1 1 1],'ycolor',palette.text,'xcolor',palette.text);
                set(control_locs,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0]);

                % % determine which is the voltage plot, set to black
                Ydata=get(ch2,'YData');
                L=cellfun(@length,Ydata');
                longlocs=L==max(L);
                set(ch2(longlocs),'color','k');

                attributes=arrayfun(@get,ch2,'UniformOutput',false);
                has_color = cellfun(@(x) isfield(x,'BackgroundColor'),attributes);
            end
            
        % features
        case 'features'
            view_mode=logical([get(handles.waveform_view_button,'value') get(handles.split_view_button,'value') get(handles.features_view_button,'value')]);
            wf_axes=getappdata(handles.output,'wf_axes');
            
            % stretch wf_axes to fill full height
            ax_waveforms=getappdata(handles.output,'ax_waveforms');
            ax_features=getappdata(handles.output,'ax_features');
            
            pos=get(ax_features,'position'); margin_top=1-(pos(2)+pos(4));
            pos=get(ax_waveforms,'position');
            if(iscell(pos)), pos=cell2mat(pos); margin_bot=min(pos(:,2));
            else margin_bot=pos(1);
            end
            total_height=1-(margin_top+margin_bot);
            num_ax=length(wf_axes);
            plot_height=total_height/num_ax;
            y=(0:num_ax-1)*plot_height;
            pos=get(wf_axes(1),'position');
            pos=repmat(pos,num_ax,1);
            pos(:,2)=margin_bot+y';
            pos(:,4)=plot_height;
            pos=pos(end:-1:1,:);
            c=copyobj(wf_axes,f);
            set(c,{'position'},num2cell(pos,2));

            % waveforms
            if(view_mode(1) || view_mode(2))
                copyobj(ax_waveforms,f);
            end
            
            % features
            if(view_mode(2) || view_mode(3))
                copyobj(ax_features,f);
            end
            
            % shrink text
            text_locs = findobj(f,'type','text');
            fontsize=get(text_locs,'fontsize');
            fontsize=[fontsize{:}];
            set(text_locs',{'fontsize'},num2cell(fontsize'*.75));
            
            
            % if black and white is checked, convert all axis children to black,
            % make plot white
            if(get(handles.screenshot_bw_checkbox,'Value'))
                ch=get(f,'children');
                control_locs=findobj(f,'type','uicontrol');
                ax_locs = findobj(f,'type','axes');
                ch2=get(ax_locs,'children');

                
                % set colors
                set(f,'color',[1 1 1]);
                set(ax_locs,'color',[1 1 1],'ycolor',[0 0 0],'xcolor',[0 0 0],'linewidth',1);
                set(control_locs,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0]);
                
                set(text_locs,'BackgroundColor',[1 1 1],'color',[0 0 0]);
            end
            
        
        % cell states
        case 'cell_stats'
            ax_cell_stats=getappdata(handles.output,'ax_cell_stats');
            wf_axes=getappdata(handles.output,'wf_axes');
            copyobj(ax_cell_stats,f);
            c2=copyobj(wf_axes,f);
            
            cellstat_pos=cell2mat(get(ax_cell_stats,'position'));
            left=min(cellstat_pos(:,1)); bottom=min(cellstat_pos(:,2)); top=max(cellstat_pos(:,2)+cellstat_pos(:,4));
            
            num_cells=length(c2);
            total_height=top-bottom; plot_height=total_height/num_cells;
            bottom_wf_ax=[.01 bottom left-.02 plot_height];
            pos=repmat(bottom_wf_ax,num_cells,1);
            pos(:,2)=pos(:,2)+((num_cells-1:-1:0)*plot_height)';
            
            set(c2,{'position'},num2cell(pos,2));
    end
    
    % export figure
    listing=ls('Msort_*.png');
    list_length=size(listing,1);
    filename=sprintf('Msort_%02.0f.png',list_length+1);
    set(f,'visible','on');
    print(f,'-dpng','-r300',filename,'-painters');
    close(f);
    set(handles.status_label,'String',sprintf('Exported %s/%s.',strrep(pwd,'\','/'),filename));
    
% zoom in on voltage plot
elseif(strcmp(eventData.Key,'equal') && strcmp(eventData.Modifier,'shift'))
    Fs=getappdata(handles.output,'Fs');
    min_window_size=round(.01*Fs); % .01 seconds max
    old_window_size=window_size;
    window_size=max(round(old_window_size*.8),min_window_size);
    voltage_slider=getappdata(handles.output,'voltage_slider');
    
    % determine new SliderStep
    slider_max=get(voltage_slider,'Max');
    slider_min=get(voltage_slider,'Min');
    slider_span=slider_max-slider_min;
    stepSize=[1 4]*window_size/(slider_span*5);
    set(voltage_slider,'SliderStep',stepSize);
    
    % re-center view
    pos=get(voltage_slider,'Value');
    middle=pos+old_window_size/2-window_size/2;
    set(voltage_slider,'Value',middle);
    
    voltage_slider_Callback(getappdata(handles.output,'voltage_slider'),[],handles);
% zoom out on voltage plot
elseif(any(strcmp(eventData.Key,'hyphen')) && any(strcmp(eventData.Modifier,'shift')));
    Fs=getappdata(handles.output,'Fs');
    max_window_size=round(5*Fs); % 5 seconds max
    old_window_size=window_size;
    window_size=min(round(old_window_size*1.25),max_window_size);
    voltage_slider=getappdata(handles.output,'voltage_slider');
    
    % determine new SliderStep
    slider_max=get(voltage_slider,'Max');
    slider_min=get(voltage_slider,'Min');
    slider_span=slider_max-slider_min;
    stepSize=[1 4]*window_size/(slider_span*5);
    set(voltage_slider,'SliderStep',stepSize);
    
    % re-center view
    pos=get(voltage_slider,'Value');
    middle=pos+old_window_size/2-window_size/2;
    set(voltage_slider,'Value',middle);
    
    voltage_slider_Callback(getappdata(handles.output,'voltage_slider'),[],handles);
end

function ctrl_release(~,eventData)
global ctrl;
if(strcmp(eventData.Key,'control'))
    ctrl=false;
end


% --- callback for the waveform selection buttons
function move_wfs(~,~,handles,dest_cluster)

% grab waveforms that are actively selected
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);
selected_waveforms=getappdata(handles.output,'selected_waveforms');
if(isempty(selected_waveforms))
    cluster_number=getappdata(handles.output,'cluster_number');
    % grab all waveforms that contain current clusters
    selected_waveforms=false(1,length(idx));
    for i = 1:length(cluster_number)
        selected_waveforms = selected_waveforms | idx==cluster_number(i);
    end
end
idx(selected_waveforms)=dest_cluster;
setappdata(handles.output,'idx',idx);
setappdata(handles.output,'selected_waveforms',[]);

fix_clusters(handles);
show_featurespanel(handles);


% --- Executes on button press in scale_checkbox.
function scale_checkbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scale_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of scale_checkbox
global waveforms;
if(get(hObject,'Value')==get(hObject,'Max'))
    setappdata(handles.output,'normalized_waves',true);
    if(get(handles.interp_waveforms_checkbox,'Value'))
        %setappdata(handles.output,'waveforms',getappdata(handles.output,'normalized_waveforms_interp'));
        waveforms=getappdata(handles.output,'normalized_waveforms_interp');
    else
        %setappdata(handles.output,'waveforms',getappdata(handles.output,'original_normalized_waveforms'));
        waveforms=getappdata(handles.output,'original_normalized_waveforms');
    end
else
    setappdata(handles.output,'normalized_waves',false);
    if(get(handles.interp_waveforms_checkbox,'Value'))
        %setappdata(handles.output,'waveforms',getappdata(handles.output,'waveforms_interp'));
        waveforms=getappdata(handles.output,'waveforms_interp');
    else
        %setappdata(handles.output,'waveforms',getappdata(handles.output,'original_waveforms'));
        waveforms=getappdata(handles.output,'original_waveforms');
    end
end

waveforms_norm=getappdata(handles.output,'normalized_waveforms');
if(isempty(waveforms_norm))
    %waveforms=getappdata(handles.output,'waveforms');
    set(handles.status_light,'BackgroundColor',palette.red);
    set(handles.status_label,'String','Normalizing waveforms...');
    waveforms_norm=normalize_waveforms(waveforms);
    setappdata(handles.output,'normalized_waveforms',waveforms_norm);
    set(handles.status_label,'String','Done!');
    set(handles.status_light,'BackgroundColor',palette.green);
end

calculate_features(handles);
show_featurespanel(handles);


% --- Executes on button press in ch_tags
function ch_tags_Callback(~, ~, handles)
% hObject    handle to ch_tags (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ch_tag_1
% grab selected channels
update_selected_channels(handles);

% --- Updates waveform plots and selected handles
function update_selected_channels(handles)
selected_channels=logical(cell2mat(get(handles.ch_tags,'Value'))');
setappdata(handles.output,'selected_channels',selected_channels);
fix_feature_axis_locations(handles);
fix_featurepanel_display(handles);

% --- Executes on button press in convert_button.
function convert_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to convert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Converting .plx to .mat...'); drawnow;
plx_file=getappdata(handles.output,'plx_file');
plx2mat(plx_file);
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green); drawnow;


% --- Executes during object creation, after setting all properties.
function fig_msort_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to fig_msort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
tmp=get(hObject,'units');
set(hObject,'units','pixels');
scrsz = get(0,'ScreenSize');
width = scrsz(3)*.9; height = scrsz(4)*.9;
left=(scrsz(3)-width)*.5; bot=(scrsz(4)-height)*.5;
set(hObject,'position',[left bot width height]);
set(hObject,'units',tmp);
opengl hardware;


% --- Executes on button press in restore_wfs_button.
function restore_wfs_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to restore_wfs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=getappdata(handles.output,'idx');
if(isempty(idx)), return; end
idx(:)=0;
setappdata(handles.output,'idx',idx);
fix_clusters(handles);
show_featurespanel(handles);


% --- Executes on button press in merge_group_button.
% Note: this is the "Export" button
function merge_group_button_Callback(~, ~, handles) %#ok<DEFNU>
global palette;
% hObject    handle to merge_group_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Grabbing spike times from groups...'); drawnow;

num_groups=getappdata(handles.output,'num_groups');
pathname=getappdata(handles.output,'pathname');
basename=getappdata(handles.output,'basename');
mat_file=getappdata(handles.output,'mat_file');
plx_file=getappdata(handles.output,'plx_file');
spike_file=sprintf('%s%s_spikes.mat',pathname,basename);


if(~exist(mat_file,'file'))
    set(handles.status_label,'String','.mat file not found.');
    set(handles.status_light,'BackgroundColor',palette.orange);
    return;
end

% for thresholds and channel listing
exp=xml2struct(getappdata(handles.output,'xml_file'));
exp = exp.experiment;

% cycle through groups, grab spikes
spikes = cell(1,0);
channels = cell(1,num_groups);
thresholds = cell(1,num_groups);
mean_waveforms = cell(1,0);
amplitudes = cell(1,0);
groups=zeros(1,0);
t=cell(1,num_groups); %#ok<NASGU>
for i = 1:num_groups
    % grab variables
    var_idx=sprintf('idx_%g',i);
    var_timestamps=sprintf('timestamps_%g',i);
    var_waveforms=sprintf('original_waveforms_%g',i);
    var_t=sprintf('t_%g',i);
    var_a=sprintf('amplitudes_%g',i);
    load(mat_file,var_idx,var_timestamps,var_waveforms,var_t,var_a);
    if(~exist(var_idx,'var') || ~exist(var_timestamps,'var'))
        continue;
    end
    eval(sprintf('idx=%s;',var_idx));
    eval(sprintf('timestamps=%s;',var_timestamps));
    eval(sprintf('waveforms=%s;',var_waveforms));
    eval(sprintf('t{%g}=%s;',i,var_t));
    eval(sprintf('amp = %s;',var_a));
    clear(var_idx,var_timestamps,var_waveforms,var_t);
    
    mean_wfs=get_mean_waveforms(waveforms,idx,false);
    
    u=unique(idx);
    if(~any(u)),continue; end
    s=cell(u(end),1);
    m=cell(u(end),1);
    a=cell(1,u(end));
    
    for j = 1:u(end)
        loc=idx==j;
        s{j}=timestamps(loc);
        m{j}=squeeze(mean_wfs(j,:,:));
        a{1,j}=amp(:,loc); %#ok<NODEF>
    end
    spikes = [spikes; s]; %#ok<AGROW>
    amplitudes = [amplitudes a]; %#ok<AGROW>
    mean_waveforms = [mean_waveforms; m]; %#ok<AGROW>
    groups(end+1:end+u(end))=i;
    
    % grab thresholds
    ch = zeros(1,length(exp.shank{i}.channel));
    enabled = false(1,length(ch));
    thresh = zeros(2,length(ch));
    for j = 1:length(ch)
        if(strcmpi(exp.shank{i}.channel{j}.enabled.Text,'true'))
            enabled(j)=true;
            thresh(1,j) = str2double(exp.shank{i}.channel{j}.threshold.Text);
            
            if(isfield(exp.shank{i}.channel{j},'threshold2'))
                thresh(2,j) = str2double(exp.shank{i}.channel{j}.threshold2.Text);
            end
            
            ch(j)=str2double(exp.shank{i}.channel{j}.Attributes.num);
        end
    end
    channels{i}=ch(enabled);
    thresholds{i}=thresh(:,enabled);
end

% Grab event data
% Determine location of events from plx file
set(handles.status_label,'String','Checking for event data...');
events=getappdata(handles.output,'events');
if(isempty(events))
    vars=whos('-file',mat_file);
    vars={vars.name};
    if(~any(strcmpi(vars,'events')))
        % find event list with longest length
        [~, ~, evcounts] = plx_info(plx_file,true);
        if(~any(evcounts))
            events=[];
            no_events=true;
        else
            no_events=false;
            [~,evloc]=max(evcounts);
            [~, evchans] = plx_event_chanmap(plx_file);
            event_channel=evchans(evloc);

            set(handles.status_label,'String','Grabbing events...');
            [~, events] = plx_event_ts(plx_file, event_channel);
            events=events';
        end
    end
else
    no_events=false;
end

% pull gaps out of events if they exist
gaps=getappdata(handles.output,'gaps');
if(~isempty(gaps))
    num_gaps = size(gaps,2);
    for i = 1:num_gaps
        gap_start=gaps(1,i); gap_stop=gaps(2,i); span=gap_stop-gap_start;
        if(span>0)
            locs=events>gap_stop;
            events(locs)=events(locs)-span;
            gaps(:,i+1:end)=gaps(:,i+1:end)-span;
        end
    end
end

% set first event to time 0, put spikes at beginning
%first_event = events(1);
%events=events-first_event; %#ok<NASGU>
%spikes=cellfun(@(x) x(x>first_event)-first_event,spikes,'UniformOutput',false); %#ok<NASGU>

% grab electrode type
exp=getappdata(handles.output,'exp');
electrode_type=exp.experiment.electrode.Text; %#ok<NASGU>

% new: build enabled/disabled list
enabled=cell(1,num_groups);
for i = 1:num_groups
    ch=exp.experiment.shank{i}.channel;
    enabled{i}=false(1,length(ch));
    for j = 1:length(ch)
        val=ch{j}.enabled.Text;
        if(strcmp(val,'true')),enabled{i}(j)=true;
        else enabled{i}(j)=false;
        end
    end    
end

% export data
save(spike_file,'spikes','events','groups','mean_waveforms','electrode_type','t','enabled','amplitudes','thresholds','channels');

if(no_events)
    set(handles.status_lable,'String',sprintf('No events found. Saved spikes to %s.',spike_file));
    set(handles.status_light,'BackgroundColor',palette.orange);
else
    set(handles.status_label,'String',sprintf('Saved spikes to %s.',spike_file));
    set(handles.status_light,'BackgroundColor',palette.green);
end



% --- Executes on button press in display_cleanup_button.
function display_cleanup_button_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to display_cleanup_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display_cleanup_button

% button is pressed
if(get(hObject,'Value')==get(hObject,'Max'))
    set(handles.errorbar_field,'enable','on');
    set(handles.errorbar_slider,'enable','on');
    set(handles.std_dev_label,'enable','on');
    set(handles.cleanup_button,'enable','on');
    
    setappdata(handles.output,'display_cleanup_tools',true);
    set(handles.display_cleanup_button,'String','Disable');
    panel_selected=getappdata(handles.output,'panel_selected');
    plot_mode=getappdata(handles.output,'feature_plot_mode');
    if(strcmp(panel_selected,'features') && (strcmp(plot_mode,'waveforms') || strcmp(plot_mode,'split')))
        fix_featurepanel_display(handles);
        % display cleanup tools (only on visible channels and clusters)
        %errorbar_plots=getappdata(handles.output,'errorbar_plots');
        %set(errorbar_plots,'visible','on');
    else
        set(getappdata(handles.output,'sphere_plots'),'visible','on');
    end
    
    % button is depressed
else
    set(handles.errorbar_field,'enable','off');
    set(handles.errorbar_slider,'enable','off');
    set(handles.std_dev_label,'enable','off');
    set(handles.cleanup_button,'enable','off');
    
    % display cleanup tools
    setappdata(handles.output,'display_cleanup_tools',false);
    set(handles.display_cleanup_button,'String','Enable');
    panel_selected=getappdata(handles.output,'panel_selected');
    plot_mode=getappdata(handles.output,'feature_plot_mode');
    if(strcmp(panel_selected,'features') && (strcmp(plot_mode,'waveforms') || strcmp(plot_mode,'split')))
        % display cleanup tools
        %errorbar_plot=getappdata(handles.output,'errorbar_plot');
        %set(errorbar_plot,'visible','off');
        fix_featurepanel_display(handles);
    else
        set(getappdata(handles.output,'sphere_plots'),'visible','off');
    end
end



function errorbar_field_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to errorbar_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of errorbar_field as text
%        str2double(get(hObject,'String')) returns contents of errorbar_field as a double
value=str2double(get(hObject,'String'));
minval=get(handles.errorbar_slider,'Min'); maxval=get(handles.errorbar_slider,'Max');
if(value<minval), value=minval;
elseif(value>maxval), value=maxval;
end
set(handles.errorbar_slider,'Value',value);
value=sprintf('%g',value);
set(hObject,'String',value);
%old_value=getappdata(handles.output,'errorbar_value');
setappdata(handles.output,'errorbar_value',value);

% update waveform panel if it's visible
panel=getappdata(handles.output,'panel_selected');
if(strcmp(panel,'features'))
    %change_errorbar(handles,value,old_value);
    update_spheres(handles);
end

% --- Executes during object creation, after setting all properties.
function errorbar_field_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to errorbar_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function errorbar_slider_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to errorbar_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=get(hObject,'Value');
maxval=get(hObject,'Max'); minval=get(hObject,'Min');
if(value>maxval)
    value=maxval;
    set(hObject,'Value',value);
elseif(value<minval)
    value=minval;
    set(hObject,'Value',value);
end
set(handles.errorbar_field,'String',sprintf('%g',value));
setappdata(handles.output,'errorbar_value',value);

% update waveform panel if it's visible
panel_selected=getappdata(handles.output,'panel_selected');
%plot_mode=getappdata(handles.output,'feature_plot_mode');
if(strcmp(panel_selected,'features'))
    %change_errorbar(handles,value);
    update_spheres(handles);
end

% --- Executes during object creation, after setting all properties.
function errorbar_slider_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to errorbar_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% set slider step
minval=get(hObject,'Min'); maxval=get(hObject,'Max');
dslide=0.25;
step=dslide/(maxval-minval);
set(hObject,'SliderStep',[step step]);


% --- changes the size of the ellipsoids around feature space, and
% re-determines the waveforms in each hull
function update_spheres(handles)
global features;
sphere_plots=getappdata(handles.output,'sphere_plots');
if(isempty(sphere_plots) || ~any(ishandle(sphere_plots)))
    return;
end
wfs_in_sphere=getappdata(handles.output,'wfs_in_sphere');
radius=get(handles.errorbar_slider,'Value');

% Determine which features we have selected
feature1=get(handles.feature1_menu,'Value');
feature2=get(handles.feature2_menu,'Value');
feature3=get(handles.feature3_menu,'Value');

if(any([feature1 feature2 feature3]==size(features,2)+1))
    features_temp=[features getappdata(handles.output,'timestamps')'];
else
    features_temp=features;
end

idx=getappdata(handles.output,'idx');
num_clusters=getappdata(handles.output,'num_clusters');
% 2-dimensional case
if(feature3==1)
    for i = 2:num_clusters
        ind=idx==i-1;
        if(~any(ind)),continue; end
        
        feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2);
        feats=[feat1 feat2];
        
        % calculate ellipse
        [Ex,Ey]=ellipsoid_cov(mean(feats),radius,cov(feats),30);
        %sphere_plots(i)=patch(Ex,Ey,cluster_colors(i,:),'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'parent',ax_features,'visible','off');
        set(sphere_plots(i),'XData',Ex','YData',Ey');
        wfs_in_sphere(ind,i)=inhull(feats,get(sphere_plots(i),'Vertices'));
    end
else
    feature3=feature3-1;
    for i = 2:num_clusters
        ind=idx==i-1;
        if(~any(ind)),continue; end
        feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2); feat3=features_temp(ind,feature3);
        feats=[feat1 feat2 feat3];
        
        % calculate ellipsoid
        if(size(feats,1)>2 && feature1 ~= feature2 && feature1 ~= feature3 && feature2 ~= feature3)
            [Ex,Ey,Ez]=ellipsoid_cov(mean(feats),radius,cov(feats),30);
            set(sphere_plots(i),'XData',Ex,'YData',Ey,'ZData',Ez);
            %sphere_plots(i)=surf(Ex,Ey,Ez,'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'facelighting','gouraud','parent',ax_features,'visible','off');
            [~,v]=surf2patch(Ex,Ey,Ez);
            wfs_in_sphere(ind,i)=inhull(feats,v);
        end
    end
end
setappdata(handles.output,'sphere_plots',sphere_plots);
setappdata(handles.output,'wfs_in_sphere',wfs_in_sphere);



% --- Changes the errorbars to a new value
% Note: this function may become deprecated as we've adopted the
% ellipsoid-method instead.
function change_errorbar(handles,errorbar_value)
global waveforms main_channel;
errorbar_plot=getappdata(handles.output,'errorbar_plot');
if(isempty(errorbar_plot)), return; end

% Note: the 'selected_channels' argument is optional, and indicates
% multiple channels should be plotted
t=getappdata(handles.output,'t');
if(isempty(waveforms)), return; end
idx=getappdata(handles.output,'idx');
cluster_number=getappdata(handles.output,'cluster_number');
wfs=waveforms(:,idx==cluster_number,:);

% if we're at an empty cluster
if(isempty(wfs)), return; end

% determine main plot channel
mean_wfs = squeeze(mean(wfs,2));
[~,main_channel]=min(min(mean_wfs,[],2));
channels=1:length(getappdata(handles.output,'channels'));
selected_channels=channels(getappdata(handles.output,'selected_channels'));
ch=selected_channels;

if(isempty(selected_channels))
    ch=main_channel;
end
if(isempty(find(ch==main_channel,1)))
    main_channel=ch(1);
end

% Determine 3 std deviation error bars (display cleanup tools)
wfs=squeeze(wfs(main_channel,:,:));
mean_wf_pts=mean(wfs);
std_dev_pts=std(wfs);
std_dev_pts = [mean_wf_pts-errorbar_value*std_dev_pts; mean_wf_pts+errorbar_value*std_dev_pts; nan(1,length(t))];

% plot std devs
std_dev_pts=reshape(std_dev_pts,1,[]);
set(errorbar_plot,'YData',std_dev_pts);


% --- Executes on button press in cleanup_button.
% Note: this has been updated to only utilize the cluster plots
function cleanup_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to cleanup_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selected_groups=getappdata(handles.output,'selected_axes');
wfs_in_sphere=getappdata(handles.output,'wfs_in_sphere');
wfs_in_sphere=wfs_in_sphere(:,selected_groups);
idx=getappdata(handles.output,'idx');

% we need to grab only those spikes in this cluster
selected_groups=find(selected_groups);
ind=false(1,length(idx));
for i = selected_groups
    ind=ind | idx==i-1;
end
bad_wfs=~any(wfs_in_sphere,2) & ind';

if(any(bad_wfs))
    idx(bad_wfs)=0;
    setappdata(handles.output,'idx',idx);
    fix_clusters(handles);
    plot_features(handles);
    plot_mean_waveforms(handles);
    show_all_waveforms(handles);
end



% --- Executes on selection change in cluster_menu.
function cluster_menu_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cluster_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cluster_menu


% --- Executes during object creation, after setting all properties.
function cluster_menu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when group_cmd_panel is resized.
function group_cmd_panel_ResizeFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to group_cmd_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in feature1_menu.
function feature1_menu_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to feature1_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns feature1_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature1_menu
selected_features=getappdata(handles.output,'selected_features');
selected_features(1)=get(hObject,'Value');
setappdata(handles.output,'selected_features',selected_features);
plot_features(handles);
fix_featurepanel_display(handles);

% --- Executes during object creation, after setting all properties.
function feature1_menu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to feature1_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature2_menu.
function feature2_menu_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to feature2_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns feature2_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature2_menu
selected_features=getappdata(handles.output,'selected_features');
selected_features(2)=get(hObject,'Value');
setappdata(handles.output,'selected_features',selected_features);
plot_features(handles);
fix_featurepanel_display(handles);

% --- Executes during object creation, after setting all properties.
function feature2_menu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to feature2_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in feature3_menu.
function feature3_menu_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to feature3_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns feature3_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from feature3_menu
selected_features=getappdata(handles.output,'selected_features');
selected_features(3)=get(hObject,'Value');
setappdata(handles.output,'selected_features',selected_features);
plot_features(handles);
fix_featurepanel_display(handles);

% --- Executes during object creation, after setting all properties.
function feature3_menu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to feature3_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in parameters_button.
function parameters_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to parameters_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(get(handles.parameters_panel,'visible'),'on'))
    set(handles.parameters_panel,'visible','off');
else
    set(handles.parameters_panel,'visible','on');
end


% --- Executes during object creation, after setting all properties.
function new_cluster_button_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to new_cluster_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function junk_button_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to junk_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when selected object is changed in view_panel.
function view_panel_SelectionChangeFcn(~, eventdata, handles) %#ok<DEFNU>
% hObject    handle to the selected object in view_panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
% switched to waveform view

switch(eventdata.NewValue)
    case handles.waveform_view_button
        setappdata(handles.output,'feature_plot_mode','waveforms');
        if(eventdata.OldValue==handles.split_view_button)
            fix_feature_axis_locations(handles);
        end
    case handles.split_view_button
        setappdata(handles.output,'feature_plot_mode','split');
        if(eventdata.OldValue==handles.waveform_view_button || eventdata.OldValue==handles.features_view_button)
            fix_feature_axis_locations(handles);
        end
    case handles.features_view_button
        setappdata(handles.output,'feature_plot_mode','features');
        if(eventdata.OldValue==handles.split_view_button)
            fix_feature_axis_locations(handles);
        end
end

fix_featurepanel_display(handles);

% --- Executes during object creation, after setting all properties.
function features_panel_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to features_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function deselect_button_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to deselect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in select_button.
function select_features(~, ~, handles)
% hObject    handle to select_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global features palette;
ax_features=getappdata(handles.output,'ax_features');


% Determine which features have been plotted
feature1=get(handles.feature1_menu,'Value');
feature2=get(handles.feature2_menu,'Value');
feature3=get(handles.feature3_menu,'Value')-1;

timestamps=getappdata(handles.output,'timestamps');
features_temp=[features timestamps'];

new_ax=axis(ax_features);
clusters_hidden=getappdata(handles.output,'clusters_hidden');
if(length(new_ax)==6)
    features_temp=features_temp(:,[feature1 feature2 feature3]);
    % replot features in 2D
    T=view(ax_features);
    % generate new axis limits
    new_ax=T*[reshape(axis(ax_features)',2,[]) ones(2,1)]';
    new_ax=reshape(new_ax(1:2,:)',1,[]);
    XYZ=T*[features_temp ones(size(features_temp,1),1)]';
    
    % we might have to do some data/axis flipping
    % flip X values down the middle
    if(new_ax(1)>new_ax(2))
        mean_x=mean(XYZ(1,:));
        d=mean_x-XYZ(1,:);
        XYZ(1,:)=mean_x+d;
        tmp=new_ax(1); new_ax(1)=new_ax(2); new_ax(2)=tmp;
    end
    
    % flip Y axis
    if(new_ax(3)>new_ax(4))
        mean_y=mean(XYZ(2,:));
        d=mean_y-XYZ(2,:);
        XYZ(2,:)=mean_y+d;
        tmp=new_ax(3); new_ax(3)=new_ax(4); new_ax(4)=tmp;
    end
    
    
    % put new axis temporarily over old axis
    delete(get(ax_features,'children'));
    p=plot(ax_features,XYZ(1,:),XYZ(2,:),'.','markersize',1,'linesmooth','on','color',palette.text);
    view(ax_features,0,90);
    axis(ax_features,new_ax);
    set(ax_features,'color',palette.background,'xtick',[],'ytick',[]);
    set(handles.output,'currentaxes',ax_features);
    
    % allow user to draw circle
    h=imfreehand(ax_features);
    locs=h.getPosition; delete(h);
    selected = inpolygon(XYZ(1,:)',XYZ(2,:)',locs(:,1),locs(:,2));
    
% 2 dimensional plot
else    
    if(clusters_hidden)
        idx=getappdata(handles.output,'idx');
        ind=false(1,length(idx));
        cluster_number=getappdata(handles.output,'cluster_number');
        for i = 1:length(cluster_number)
            ind = ind | idx==cluster_number(i);
        end
        features_temp=features_temp(ind,[feature1 feature2]);
    else
        features_temp=features_temp(:,[feature1 feature2]);
    end
    %set(handles.output,'currentaxes',ax_features);
    
    ax=axis(ax_features);
    delete(get(ax_features,'children'));
    p=plot(ax_features,features_temp(:,1),features_temp(:,2),'.','markersize',1,'linesmooth','on','color',palette.text);
    axis(ax_features,ax);
    %set(handles.output,'currentaxes',ax_features);
    
    % allow user to draw circle around features
    h=imfreehand(ax_features);
    locs=h.getPosition; delete(h);
    selected=inpolygon(features_temp(:,1),features_temp(:,2),locs(:,1),locs(:,2));
end

% highlight selected points
colors=[palette.text; palette.red];
X=get(p,'Xdata'); Y=get(p,'Ydata');
color_ind=ones(1,length(X));
color_ind(selected)=2;

delete(p);
colormap(ax_features,colors);
scatter(ax_features,X,Y,1,color_ind');

if(clusters_hidden)
    s=1:length(idx); s=s(ind); s(~selected)=[];
    selected=false(1,length(idx)); selected(s)=true;
end

setappdata(handles.output,'selected_waveforms',selected);

% enable merge buttons
merge_buttons=getappdata(handles.output,'merge_buttons');
set(merge_buttons,'enable','on');
cluster_colors=getappdata(handles.output,'cluster_colors');
for i = 1:length(merge_buttons)
    set(merge_buttons(i),'backgroundcolor',cluster_colors(i,:));
end

% --- Executes on button press in deselect_button.
function deselect_button_Callback(~, ~, handles)
% hObject    handle to deselect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.output,'selected_waveforms',[]);
selected_wf_plots=getappdata(handles.output,'selected_wf_plots');
delete(selected_wf_plots(ishandle(selected_wf_plots) & selected_wf_plots~=0));
selected_wf_plots(:)=0;
setappdata(handles.output,'selected_wf_plots',selected_wf_plots);

colors=get_cluster_colors(handles);
buttons=getappdata(handles.output,'merge_buttons');
for i = 1:length(buttons)
    if(strcmp(get(buttons(i),'enable'),'on'))
        set(buttons(i),'backgroundcolor',colors(i,:)/3,'enable','off');
    end
end

% --- Executes on button press in new_cluster_button.
function new_cluster_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to new_cluster_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
create_new_cluster(handles);


% --- Executes on button press in select_wfs_button.
function select_wfs_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to select_wfs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
select_wfs(handles);


% --- Executes on button press in junk_button.
function junk_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to junk_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
move_wfs([],[],handles,0);


% --- Executes on button press in sound_checkbox.
function sound_checkbox_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to sound_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sound_checkbox


% --- Executes on button press in stretch_checkbox.
function stretch_checkbox_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to stretch_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stretch_checkbox
panel_selected=getappdata(handles.output,'panel_selected');
plot_mode=getappdata(handles.output,'feature_plot_mode');
if(strcmp(panel_selected,'features') && (strcmp(plot_mode,'waveforms') || strcmp(plot_mode,'split')))
    show_all_waveforms(handles);
end


% --- Executes on button press in template_button.
% Note: this function uses clusters from a separate recording and uses
% template matching. All waveforms are matched to current templates. It is
% up to the user to create new clusters, should they arise.
function template_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to template_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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


% --- Executes on button press in apply_templates_button.
function apply_templates_button_Callback(~, ~, handles)
% hObject    handle to apply_templates_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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


% --- Executes on button press in cluster_options_button.
function cluster_options_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to cluster_options_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(get(handles.cluster_panel,'visible'),'on'))
    set(handles.cluster_panel,'visible','off');
else
    set(handles.cluster_panel,'visible','on');
end


% --- Executes on button press in link_channels_checkbox.
function link_channels_checkbox_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to link_channels_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of link_channels_checkbox
if(get(hObject,'Value')==get(hObject,'Max'))
    PCs_fields=get(handles.num_PCs_slider);
    peak_checkbox=get(handles.peak_checkbox);
    valley_checkbox=get(handles.valley_checkbox);
    slope_checkbox=get(handles.slope_checkbox);
    energy_checkbox=get(handles.energy_checkbox);
    amplitude_checkbox=get(handles.amplitude_checkbox);
    
    PCs=[PCs_fields.Value]; set(handles.num_PCs_slider,'Value',max(PCs)); set(handles.num_PCs_field,'String',sprintf('%2.0f',max(PCs)));
    peaks=[peak_checkbox.Value]; set(handles.peak_checkbox,'Value',max(peaks));
    valley=[valley_checkbox.Value]; set(handles.valley_checkbox,'Value',max(valley));
    slope=[slope_checkbox.Value]; set(handles.slope_checkbox,'Value',max(slope));
    energy=[energy_checkbox.Value]; set(handles.energy_checkbox,'Value',max(energy));
    amplitude=[amplitude_checkbox.Value]; set(handles.amplitude_checkbox,'Value',max(amplitude));
end

% --- Executes on button press in close_clusterpanel_button.
function close_clusterpanel_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to close_clusterpanel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.cluster_panel,'visible','off');


% --- Executes on selection change in group_menu.
function group_menu_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to group_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns group_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group_menu
loc=get(hObject,'value'); loc=loc-1;
if(loc>0)
    load_group(loc,handles);
end

% --- Executes during object creation, after setting all properties.
function group_menu_CreateFcn(hObject, eventdata, handles) %#ok<INUSD,DEFNU>
% hObject    handle to group_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interp_waveforms_checkbox.
function interp_waveforms_checkbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to interp_waveforms_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waveforms;

% Hint: get(hObject,'Value') returns toggle state of interp_waveforms_checkbox
if(get(hObject,'Value'))
    %setappdata(handles.output,'waveforms',getappdata(handles.output,'waveforms_interp'));
    waveforms=getappdata(handles.output,'waveforms_interp');
    setappdata(handles.output,'timestamps',getappdata(handles.output,'timestamps_interp'));
    setappdata(handles.output,'t',getappdata(handles.output,'t_interp'));
    setappdata(handles.output,'normalized_waveforms',getappdata(handles.output,'normalized_waveforms_interp'));
else
    %setappdata(handles.output,'waveforms',getappdata(handles.output,'original_waveforms'));
    waveforms=getappdata(handles.output,'original_waveforms');
    setappdata(handles.output,'timestamps',getappdata(handles.output,'original_timestamps'));
    setappdata(handles.output,'t',getappdata(handles.output,'original_t'));
    setappdata(handles.output,'normalized_waveforms',getappdata(handles.output,'original_normalized_waveforms'));
end

% if we're currently viewing waveforms, update
panel_selected=getappdata(handles.output,'panel_selected');
if(strcmp(panel_selected,'features'))
    plot_mean_waveforms(handles);
    show_all_waveforms(handles);
    fix_featurepanel_display(handles);
end
if(strcmp(panel_selected,'cell_stats'))
    plot_mean_waveforms(handles);
end


% --- Executes on button press in undo_button.
function undo_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to undo_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=getappdata(handles.output,'idx_old');
setappdata(handles.output,'idx',idx);

u=unique([0 idx]); num_clusters=length(u);
setappdata(handles.output,'num_clusters',num_clusters);
cluster_colors=get_cluster_colors(handles);
setappdata(handles.output,'cluster_colors',cluster_colors);

panel_selected=getappdata(handles.output,'panel_selected');
setappdata(handles.output,'cluster_number',0);
if(strcmp(panel_selected,'features'))
    plot_mean_waveforms(handles);
    
    % update features
    show_featurespanel(handles);
elseif(strcmp(panel_selected,'cell_stats'))
    stats_button_Callback([],[],handles);
end

% returns cluster colors
function cluster_colors = get_cluster_colors(handles)
global palette;
num_clusters=getappdata(handles.output,'num_clusters');
cluster_colors=zeros(num_clusters,3);

tmp=palette.cluster_colors;
cluster_colors(1,:)=[.5 .5 .5];
cluster_colors(2:end,:)=tmp(1:num_clusters-1,:);

% returns channel colors
function channel_colors = get_channel_colors(num_channels)
global palette;

cmap=[palette.blue;palette.red;palette.orange;palette.green; palette.teal];
channel_colors=get_colors(num_channels,cmap);
setrand(1);
channel_colors=channel_colors(randperm(num_channels),:);


% --- Executes on button press in clear_misaligned_button.
function clear_misaligned_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to clear_misaligned_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

% Re-distributes waveform axis locations (occurs when channels are selected
% and deselected)
function fix_feature_axis_locations(handles)
global palette;

channels=getappdata(handles.output,'channels');
num_channels=length(channels);

ax_waveforms=getappdata(handles.output,'ax_waveforms');

% Create new axes if they've been deleted
if(isempty(ax_waveforms) || ~any(ishandle(ax_waveforms)))
    ax_waveforms=zeros(1,num_channels);
    for i = 1:num_channels
        ax_waveforms(i)=axes('parent',handles.features_panel,'color',palette.ax_background,'visible','off',...
            'xcolor',palette.text,'ycolor',palette.text,'xtick',[],'ytick',[],'linewidth',1.5,'ButtonDownFcn',{@select_wfs,handles}'); %#ok<LAXES>
        hold(ax_waveforms(i),'all');
    end
end
setappdata(handles.output,'ax_waveforms',ax_waveforms);

% determine how many axes to make
selected_channels=getappdata(handles.output,'selected_channels');
if(~any(selected_channels)), selected_channels=~selected_channels; end
locs=find(selected_channels);
num_to_display=sum(selected_channels);

plot_mode=getappdata(handles.output,'feature_plot_mode');

% determine how many plots to make, set up basic location parameters
if(num_to_display<3), num_cols=1; else num_cols=2; end
num_rows=ceil(num_to_display/num_cols);
margin_left=0.18; margin_right=0;
margin_bot=.04;
margin_top=0.015;

% in split mode, waveform axis goes down
if(strcmp(plot_mode,'split'))
    margin_middle=.05;
    remaining=1-margin_top-margin_bot+margin_middle;
    remaining=remaining/2;
    margin_top=margin_top+remaining;
else margin_top=0.015;
end
plot_width=(1-margin_left-margin_right)/num_cols;
plot_height=(1-margin_bot-margin_top)/num_rows;

% position axes
ind=1;
for i = 1:num_rows
    y = 1-margin_top-i*plot_height;
    for j = 1:num_cols
        if(ind>length(locs)),continue; end
        x=margin_left+(j-1)*plot_width;
        set(ax_waveforms(locs(ind)),'position',[x y plot_width plot_height]);
        ind=ind+1;
    end
end

% in split mode, features axis goes up
if(strcmp(plot_mode,'split'))
    margin_top=.015;
    remaining=1-margin_top-margin_bot-margin_middle;
    remaining=remaining/2;
    margin_bot=1-margin_top-remaining;
else margin_top=.015; margin_bot=.04;
end

% ax_features is empty, need to create
ax_features=getappdata(handles.output,'ax_features');
if(isempty(ax_features) || ~ishandle(ax_features))
    ax_features=axes('parent',handles.features_panel,'color',palette.ax_background,'visible','off',...
        'xcolor',palette.text,'ycolor',palette.text,'xtick',[],'ytick',[],'linewidth',1.5,'box','on',...
        'position',[margin_left margin_bot (1-margin_left-margin_right) (1-margin_bot-margin_top)]);
    hold(ax_features,'all');
    grid(ax_features,'on');
    %set(ax_features,'ButtonDownFcn',{@select_features,handles});
    setappdata(handles.output,'ax_features',ax_features);
else
    set(ax_features,'position',[margin_left margin_bot (1-margin_left-margin_right) (1-margin_bot-margin_top)]);
end



% --- Correctly displays what is supposed to be displayed on the feature panel.
% Note: a lot of the features set here may end up being redundantly set.
function fix_featurepanel_display(handles)
margin_left=.225; margin_right=0;
margin_bot=.04; margin_top=0;

% reduces the number of "set" commands
visibles=[]; invisibles=[];

% create ax_waveforms if missing
ax_waveforms=getappdata(handles.output,'ax_waveforms');
channels=getappdata(handles.output,'channels');
num_channels=length(channels);
if(isempty(ax_waveforms) || ~any(ishandle(ax_waveforms)))
    ax_waveforms=zeros(1,num_channels);
    for i = 1:num_channels
        ax_waveforms(i)=axes('parent',handles.features_panel,'color',palette.ax_background,'visible','off',...
            'xcolor',palette.text,'ycolor',palette.text,'xtick',[],'ytick',[],'linewidth',1.5); %#ok<LAXES>
        hold(ax_waveforms(i),'all');
    end
end
setappdata(handles.output,'ax_waveforms',ax_waveforms);

% create ax_features if missing
ax_features=getappdata(handles.output,'ax_features');
if(isempty(ax_features) || ~ishandle(ax_features))
    ax_features=axes('position',[margin_left margin_bot 1-margin_left-margin_right 1-margin_top-margin_bot],...
        'parent',handles.features_panel,'color',palette.ax_background,'visible','off'); hold(ax_features,'all');
    grid(ax_features,'on');
    setappdata(handles.output,'ax_features',ax_features);
end

% determine axis selection type
view_mode=logical([get(handles.waveform_view_button,'value') get(handles.split_view_button,'value') get(handles.features_view_button,'value')]);

% determine which channels are selected
selected_channels=getappdata(handles.output,'selected_channels');
if(~any(selected_channels)), selected_channels=~selected_channels; end

% determine which groups are selected
selected_groups=getappdata(handles.output,'selected_axes');

% hide everything in channels that are not selected
hidden_axes=ax_waveforms(~selected_channels);
children=get(hidden_axes,'children');
if(iscell(children)), children=cell2mat(children); end
invisibles=[invisibles hidden_axes(:)' children(:)'];

% fix wf axis limits
min_wfs=getappdata(handles.output,'min_wfs');
max_wfs=getappdata(handles.output,'max_wfs');
checked=get(handles.fit_wf_checkbox,'Value');

% fit waveforms to selected channels
if(checked && any(selected_groups))
    min_wf=min_wfs(selected_channels,selected_groups);
    min_wf=min(min_wf(:));
    max_wf=max_wfs(selected_channels,selected_groups);
    max_wf=max(max_wf(:));
else
    min_wf=min(min_wfs(:));
    max_wf=max(max_wfs(:));
end
if(min_wf==0 || isempty(min_wf)),min_wf=0; max_wf=1; end
set(ax_waveforms,'ylim',[min_wf max_wf]);

% adjust channel labels as well
t=getappdata(handles.output,'t')*1e6;
if(isempty(t)), return; end
ch_title_plots=getappdata(handles.output,'ch_title_plots');
spanx=t(end)-t(1); spany=max_wf-min_wf;
set(ch_title_plots,'position',[t(1)+spanx*.01 min_wf+spany*.98]);

% Stuff to make visible:
%   -waveform plots (only selected groups, on selected channel axes)
%   -gridlines (only selected groups, on selected channel axes)
if(view_mode(1) || view_mode(2))
    if(view_mode(1))
        % turn off features
        children=get(ax_features,'children');
        if(iscell(children)),children=cell2mat(children); end
        invisibles=[invisibles ax_features children(:)'];
    end
    
    % waveform plot stuff
    visibles=[visibles ax_waveforms(selected_channels)];
    
    wf_plots=getappdata(handles.output,'wf_plots');
    
    if(~isempty(wf_plots) && any(ishandle(wf_plots(:))))
        gridline_plots=getappdata(handles.output,'gridline_plots');
        errorbar_plots=getappdata(handles.output,'errorbar_plots');
        
        w1=wf_plots(~selected_channels,:);
        w2=wf_plots(:,~selected_groups);
        w3=gridline_plots(~selected_channels);
        w4=ch_title_plots(~selected_channels);
        invisibles=[invisibles w1(:)' w2(:)' w3 w4(:)'];
        
        w1=wf_plots(selected_channels,selected_groups);
        w2=gridline_plots(selected_channels);
        w3=ch_title_plots(selected_channels);
        visibles=[visibles w1(:)' w2 w3(:)'];
        
        
        
        if(getappdata(handles.output,'display_cleanup_tools'))
            w1=errorbar_plots(selected_channels,selected_groups);
            visibles=[visibles w1(:)'];
            
            w1=errorbar_plots(~selected_channels,:);
            w2=errorbar_plots(:,~selected_groups);
            
            invisibles=[invisibles w1(:)' w2(:)'];
        else
            invisibles=[invisibles errorbar_plots(:)'];
        end
    end
end

%   feature plots (highlight selected groups)
if(view_mode(2) || view_mode(3))
    if(view_mode(3))
        children=get(ax_waveforms,'children');
        if(iscell(children)),children=cell2mat(children); end
        invisibles=[invisibles ax_waveforms children(:)'];
    end
    
    children=get(ax_features,'children');
    if(iscell(children)),children=cell2mat(children); end
    visibles=[visibles ax_features children(:)'];
    
    sphere_plots=getappdata(handles.output,'sphere_plots');
    if(~getappdata(handles.output,'display_cleanup_tools'))
        invisibles=[invisibles sphere_plots(:)'];
    end
    
    highlight_selected_cluster(handles);
end

set(visibles,'visible','on');
set(invisibles,'visible','off');

if(view_mode(2) || view_mode(3))
    hide_clusters(handles);
end

% fix rotate3d
if(view_mode(1))
    set(handles.output,'windowbuttonmotionfcn',[]);
else
    % make sure we're in 3D view
    twodee=get(handles.feature3_menu,'Value')==1;
    if(~twodee)
        set(handles.output,'windowbuttonmotionfcn',{@enable_rotate,get(handles.features_panel,'position'),ax_features,getappdata(handles.output,'rotate_handle')});
    end
end


% --- Executes on button press in fit_wf_checkbox.
function fit_wf_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to fit_wf_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fit_wf_checkbox
panel_selected=getappdata(handles.output,'panel_selected');
if(strcmp(panel_selected,'features'))
    ax_waveforms=getappdata(handles.output,'ax_waveforms');
    if(~isempty(ax_waveforms) && all(ishandle(ax_waveforms)))
        % determine which channels are selected
        selected_channels=getappdata(handles.output,'selected_channels');
        if(~any(selected_channels)), selected_channels=~selected_channels; end

        % determine which groups are selected
        selected_groups=getappdata(handles.output,'selected_axes');
        
        % fix wf axis limits
        min_wfs=getappdata(handles.output,'min_wfs');
        max_wfs=getappdata(handles.output,'max_wfs');
        checked=get(handles.fit_wf_checkbox,'Value');

        % fit waveforms to selected channels
        if(checked)
            min_wf=min_wfs(selected_channels,selected_groups);
            min_wf=min(min_wf(:));
            max_wf=max_wfs(selected_channels,selected_groups);
            max_wf=max(max_wf(:));
        else
            min_wf=min(min_wfs(:));
            max_wf=max(max_wfs(:));
        end
        set(ax_waveforms,'ylim',[min_wf max_wf]);
        
        % adjust channel labels as well
        t=getappdata(handles.output,'t')*1e6;
        ch_title_plots=getappdata(handles.output,'ch_title_plots');
        spanx=t(end)-t(1); spany=max_wf-min_wf;
        set(ch_title_plots,'position',[t(1)+spanx*.01 min_wf+spany*.98]);
    end
end


% --------------------------------------------------------------------
function edit_menu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edit_menu_options_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.options_panel,'visible','on');
set(handles.output,'WindowButtonMotionFcn',[]);

% --- Executes on button press in options_savebutton.
function options_savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to options_savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.options_panel,'visible','off');

% check to see if our filedir changed
if(getappdata(handles.output,'defaultdirchanged'))
    Msort_path=which('Msort.m');
    p=regexp(Msort_path,strrep(sprintf('(?<path>.*%c).*.m',filesep),'\','\\'),'names');
    param_file=sprintf('%sMsort.cfg',p.path);
    params=load_params(param_file);
    params.default_file_dir=get(handles.default_data_field,'String');
    keep_header = true; % retain the header in the .cfg file
    save_params(params,param_file,keep_header);
    setappdata(handles.output,'defaultdirchanged',false);
end

% update t_before and t_after
t_before=str2double(get(handles.t_before_field,'String'));
t_after=str2double(get(handles.t_after_field,'String'));
if(~isfinite(t_before) || t_before<1)
    t_before=300;
    set(handles.t_before,'String',sprintf('%g',t_before));
end
if(~isfinite(t_after) || t_after<1)
    t_after=500;
    set(handles.t_after,'String',sprintf('%g',t_after));
end
setappdata(handles.output,'t_before',t_before);
setappdata(handles.output,'t_after',t_after);

if(getappdata(handles.output,'slope_edited'))
    calculate_features(handles);
    setappdata(handles.output,'slope_edited',false);
end

% check to see if we need to replot voltage data
if(strcmp(getappdata(handles.output,'panel_selected'),'voltage'))
    plot_voltages(handles);
    fix_featurepanel_display(handles);
end

if(strcmp(getappdata(handles.output,'panel_selected'),'features'))
    plot_features(handles);
    fix_featurepanel_display(handles);
end




% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in square_voltage_checkbox.
function square_voltage_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to square_voltage_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of square_voltage_checkbox


% --- Executes on button press in wfs_only_checkbox.
function wfs_only_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to wfs_only_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wfs_only_checkbox



function t_before_field_Callback(hObject, eventdata, handles)
% hObject    handle to t_before_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_before_field as text
%        str2double(get(hObject,'String')) returns contents of t_before_field as a double


% --- Executes during object creation, after setting all properties.
function t_before_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_before_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_after_field_Callback(hObject, eventdata, handles)
% hObject    handle to t_after_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_after_field as text
%        str2double(get(hObject,'String')) returns contents of t_after_field as a double


% --- Executes during object creation, after setting all properties.
function t_after_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_after_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function voltage_range_field_Callback(hObject, eventdata, handles)
% hObject    handle to voltage_range_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voltage_range_field as text
%        str2double(get(hObject,'String')) returns contents of voltage_range_field as a double
val=round(2*str2double(get(hObject,'String')))/2;
minVal=get(handles.voltage_range_slider,'Min');
maxVal=get(handles.voltage_range_slider,'Max');
if(val<minVal),val=minVal; end
if(val>maxVal),val=maxVal; end
set(hObject,'String',val);
set(handles.voltage_range_slider,'Value',val);

% --- Executes during object creation, after setting all properties.
function voltage_range_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voltage_range_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function voltage_range_slider_Callback(hObject, eventdata, handles)
% hObject    handle to voltage_range_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=round(2*get(hObject,'Value'))/2;
set(hObject,'Value',val);
set(handles.voltage_range_field,'String',sprintf('%g',val));



% --- Executes during object creation, after setting all properties.
function voltage_range_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voltage_range_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in hist_optimize_checkbox.
function hist_optimize_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to hist_optimize_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hist_optimize_checkbox



function max_wfs_field_Callback(hObject, eventdata, handles)
% hObject    handle to max_wfs_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_wfs_field as text
%        str2double(get(hObject,'String')) returns contents of max_wfs_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val))
    val=5000;
    set(hObject,'String',val);
elseif(val<0)
    val=100;
    set(hObject,'String',val);
end


% --- Executes during object creation, after setting all properties.
function max_wfs_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_wfs_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function edit_menu_rethreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_rethreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global palette;
run_threshold_spikes(handles);
plot_voltages(handles);
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);


% --- Executes when selected object is changed in options_display_panel.
function options_display_panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in options_display_panel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if(eventdata.NewValue==handles.opengl_radiobutton)
    set(handles.output,'renderer','opengl');
else
    set(handles.output,'renderer','zbuffer');
end


% --- Thresholds spikes
function [timestamps, trigger_ch] = M_threshold_spikes(handles, threshold, t_before, t_after, abs_data)
%--------------------------------------------------------------------------
% threshold_spikes.m - Given set of continuous voltage traces, thresholds
% according to threshold.  The algorithm is based on that found in
% find_peaks.m.
%
% Usage: [timestamps, trigger_ch] = threshold_spikes(data, Ts, threshold, t_before, t_after, abs_data);
%
% Input:  data            * CxD data matrix consisting of C channels and D
%                           data points
%         Ts              * 1xD vector of time stamps associated with data
%         threshold       * spike discriminator threshold value in standard
%                           deviations (2xC)
%         t_before        * pre-spike duration for spike extraction
%         t_after         * post-spike duration for spike extraction
%         squared_data*   * (optional) if this is included, then thresholds
%                           are applied to the square data, but waveforms
%                           extracted from the regular data
%
% Output: timestamps      * CxW matrix of timestamps associated with spikes
%         trigger_ch      * list channels associated with each waveform
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
global data_v Ts palette;

num_channels = size(data_v,2);

if(~exist('threshold','var')), threshold=4; end
if(~exist('t_before','var')), t_before = 300e-6; end
if(~exist('t_after','var')), t_after=500e-6; end

dt=mean(diff(Ts));
freq=round(1/dt);

% Update: We've since removed this and added it in later
% determine deadzone span for a single channel
same_ch_deadzone=str2double(get(handles.same_ch_deadzone_field,'String'))*1e-6;
diff_ch_deadzone=str2double(get(handles.diff_ch_deadzone_field,'String'))*1e-6;
same_ch_deadzone_span=round(same_ch_deadzone/dt);
diff_ch_deadzone_span=round(diff_ch_deadzone/dt);

% run through each channel individually
spikelocs=cell(1,num_channels);
peaks = cell(1,num_channels);
if(nargout>1)
    h=waitbar(0,sprintf('Thresholding %g/%g...',0,num_channels),'color',palette.background);
    set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
end
for i = 1:num_channels
    peaks1=[]; peaks2=[];
    if(nargout>1),waitbar(i/num_channels,h,sprintf('Thresholding %g/%g...',i,num_channels));end
    % top threshold first
    if(threshold(1,i)~=0)
        if(exist('abs_data','var'))
            [peaks1,spikelocs1]=findpeaks2(abs_data(:,i)','minpeakheight',threshold(1,i));
        else
            [peaks1,spikelocs1]=findpeaks2(data_v(:,i)','minpeakheight',threshold(1,i));
        end
    else spikelocs1=[];
    end
    
    % bottom threshold
    if(threshold(2,i)~=0)
        [peaks2,spikelocs2]=findpeaks2(-data_v(:,i)','minpeakheight',-threshold(2,i));
    else spikelocs2=[];
    end
    
    % if we have a string of spikes that are within the deadzone, choose
    % the spike with the greatest amplitude
    [spikelocs{i},ind]=unique([spikelocs1 spikelocs2]);
    if(length(spikelocs{i})>2)
        peaks{i}=[peaks1 peaks2]; peaks{i}=abs(peaks{i}(ind));
        badlocs=diff(spikelocs{i})<same_ch_deadzone_span;
        d=diff([false badlocs]);
        starts=find(d==1); stops=find(d==-1);
        if(length(stops)<length(starts)),stops=[stops starts(end)+1]; end %#ok<AGROW>
        ranges=stops-starts;
        maxpeaks=zeros(1,length(starts));
        for j = 1:length(starts)
            ind=starts(j):starts(j)+ranges(j);
            [~,maxpeaks(j)]=max(peaks{i}(ind));
        end
        maxpeaks=maxpeaks+starts-1;

        b = [false badlocs] | [badlocs false]; b(maxpeaks)=false;
        spikelocs{i}(b)=[];
    end
end
if(nargout>1),waitbar(1,h,'Merging results...'); end

% if we have a string of spikes that are within the cross-channel deadzone,
% choose the spike with the greatest amplitude
spike_nums=cellfun(@length,spikelocs);
tc=cell(1,num_channels);
for i = 1:num_channels
    tc{i}=ones(1,spike_nums(i))*i;
end
s=[spikelocs{:}]; p=[peaks{:}]; trigger_ch=[tc{:}];
[spikelocs,ind]=unique(s); peaks=p(ind); trigger_ch=trigger_ch(ind);
if(length(spikelocs)>2)
    badlocs=diff(spikelocs)<diff_ch_deadzone_span;
    d=diff([false badlocs]);
    starts=find(d==1); stops=find(d==-1);
    if(length(stops)<length(starts)),stops=[stops starts(end)+1]; end
    ranges=stops-starts;
    maxpeaks=zeros(1,length(starts));
    for i = 1:length(starts)
        ind=starts(i):starts(i)+ranges(i);
        [~,maxpeaks(i)]=max(peaks(ind));
    end
    maxpeaks=maxpeaks+starts-1;
    
    b = [false badlocs] | [badlocs false]; b(maxpeaks)=false;
    spikelocs(b)=[];
    trigger_ch(b)=[];
end


% determine number of points per wave
pts_before=round(t_before*freq);
pts_after=round(t_after*freq);

% remove spikes at beginning and end
locs=spikelocs<=pts_before | spikelocs>length(Ts)-pts_after;
spikelocs(locs)=[]; trigger_ch(locs)=[];

% get actual spike times
timestamps = Ts(spikelocs);

delete(h);


% --- grabs spike amplitudes ---
function amplitudes = get_spike_amplitudes(handles)
waveforms=getappdata(handles.output,'waveforms');
t=getappdata(handles.output,'t');
amplitudes=squeeze(waveforms(:,:,t==0));


% --- changes the electrode configuration. Note: this wipes waveforms (user
% is informed)
function edit_menu_electrode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_electrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
panel=uipanel('Parent',handles.output,'position',[.4 .4 .2 .2],'units','normalized');

% create and populate electrode list
[names,configurations,geometry]=load_electrode_types('electrodes.cfg');
uicontrol('Style','text','String','Select Electrode Configuration','fontsize',14,'fontweight','bold',...
    'Parent',panel,'units','normalized','position',[.01 .8 .98 .18],'Horizontalalignment','center');

current_electrode=getappdata(handles.output,'electrode');
loc=strcmp(names,current_electrode);
electrode_menu=uicontrol('Style','popupmenu','String',names,'Parent',panel,...
    'units','normalized','position',[.01 .6 .98 .18],'Value',find(loc));


% "Custom..." button
uicontrol('Style','pushbutton','parent',panel,'units','normalized',...
    'Callback',{@custom_popup,panel,handles},...
    'position',[.01 .5 .25 .13],'String','Custom...','fontsize',10);

% save button
uicontrol('style','pushbutton','parent',panel,'units','normalized',...
    'CallBack',{@change_electrode,electrode_menu,panel,configurations,geometry,handles},'position',[.6 .01 .18 .13],'String','Save','fontsize',12,...
    'fontweight','bold');

% cancel button
uicontrol('style','pushbutton','parent',panel,'units','normalized',...
    'CallBack',{@delete,panel},'position',[.8 .01 .18 .13],'String','Cancel','fontsize',12);



% --- Show custom electrode popup window
function custom_popup(~,~,panel,handles)
delete(panel);
custom_panel=uipanel('Parent',handles.output,'position',[.4 .4 .2 .2],'units','normalized');
uicontrol('Style','text','String','Custom Electrode Configuration','fontsize',14,'fontweight','bold',...
    'Parent',custom_panel,'units','normalized','position',[.01 .8 .98 .18],'Horizontalalignment','center');


% create radio button group for spikes or voltage trace
auto_buttongroup = uibuttongroup('units','normalized','position',[.01 .5 .98 .3],'parent',custom_panel);

% automatic label
uicontrol('Style','text','String','Automatically compares channels to determine configuration','fontsize',10,...
    'Parent',auto_buttongroup,'units','normalized','position',[0 .6 1 .3],'Horizontalalignment','left',...
    'fontangle','italic');

% spikes radio button
use_spikes = uicontrol('Style','radiobutton','parent',auto_buttongroup,'units','normalized',...
    'position',[.4 0 .2 .5],'Value',1,'String','Spikes','fontsize',12);

% voltage trace radio button
use_voltage = uicontrol('Style','radiobutton','parent',auto_buttongroup,'units','normalized',...
    'position',[.7 0 .2 .5],'Value',0,'String','Voltage','fontsize',12);

% automatic button
auto_button = uicontrol('Style','pushbutton','parent',auto_buttongroup,'units','normalized',...
    'Callback',{@electrode_autoconfigure,use_spikes,custom_panel,handles},...
    'position',[0 0 .3 .5],'String','Automatic','fontsize',14,'horizontalalignment','left');

% Manual button group
manual_buttongroup = uibuttongroup('units','normalized','position',[.01 .2 .98 .31],'parent',custom_panel);

% manual label
uicontrol('Style','text','String','Manually input electrode configuration','fontsize',10,...
    'Parent',manual_buttongroup,'units','normalized','position',[0 .6 1 .3],'Horizontalalignment','left',...
    'fontangle','italic');

% manual button
manual_button = uicontrol('Style','pushbutton','parent',manual_buttongroup,'units','normalized',...
    'Callback',{@electrode_manualconfigure,custom_panel,handles},...
    'position',[0 0 .3 .5],'String','Manual','fontsize',14,'horizontalalignment','left');

% cancel button
uicontrol('style','pushbutton','parent',custom_panel,'units','normalized',...
    'CallBack',{@delete,custom_panel},'position',[.8 .01 .18 .13],'String','Cancel','fontsize',12);


% --- Automatically determines electrode configuration, using either spike
%     coincidence or voltage trace correlation
function electrode_autoconfigure(~,~,use_spikes,custom_panel,handles)
global data_v Ts waveforms;
use_spikes=logical(get(use_spikes,'Value'));

% make sure user realizes that this will wipe data
choice = questdlg('Automatic detection of Electrode Groups may take a while, and will wipe the current sorting data. Are you sure you want to continue?',...
      'Change Electrode Configuration','Yes','No','No');

switch(choice)
    case 'Yes'
        % free up some memory for this
        return_button_Callback([],[],handles);
        clearvars -global -except data_v Ts;
        data_v = [];
        
        xml_file = getappdata(handles.output,'xml_file');
        mat_file = getappdata(handles.output,'mat_file');
        
        % grab channels from all variables named data_*
        names = whos('-file',mat_file,'-regexp','data_(\d)*$');        
        names = {names(:).name};
        result = regexp(names,'data_(?<ch>\d*)','names');
        result = [result{:}];
        [channels,ind] = sort(str2double({result.ch}));
        names=names(ind);
        
        N = length(channels);
        
        % this will be our correlation matrix
        correlations = zeros(N);
        
        % determine spike coincidence for each channel, using a 5-std-dev-threshold
        thresh = 5;
        if(use_spikes)
            load(mat_file,'Ts');
            spikes = cell(N,1);
            h=waitbar(0,'Processing channels...','color',palette.background);
            set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
            for i = 1:N
                waitbar(i/N,h);
                load(mat_file,names{i});
                eval(sprintf('data_v = %s''; clear %s;',names{i},names{i}));
                spikes{i} = M_threshold_spikes(handles,[thresh;-thresh],0,0);
            end
            delete(h);
            
            % determine the typical coincidence rate for Poisson spikes
            % with similar firing rates
            duration=max(Ts(end),1e5);
            firing_rate=mean(cellfun(@length,spikes))/(Ts(end)-Ts(1));
            fake_spikes=num2cell(sort(rand(2,duration*firing_rate)*duration,2),2);
            
            c = analyze_coincidence(fake_spikes); c=c(1,2);
            correlations = analyze_coincidence(spikes);
            
            % *** TEMPORARY FOR LGN RECORDING PURPOSES ONLY ***
            
            % generate proper channel list
            % note: temporary for LGN
            ind=cell2mat({[11 12 6 5],[13 14 4 3],[15 16 2 1],[25 26 24 23],[27 28 22 21],[29 39 20 19],[31 32 18 17],[10 9 7 8]});
            deletelocs=false(1,length(ind));
            for i = 1:length(ind)
                if(~any(channels==ind(i)))
                    deletelocs(i)=true;
                end
            end
            ind(deletelocs)=[];
            
            ind2=zeros(1,length(ind));
            for i = 1:length(ind2)
                ind2(i)=find(channels==ind(i));
            end
            
            fig_std;
            imagesc(correlations(ind2,ind2));
            set(gca,'xtick',1:length(ind2),'xticklabel',num2str(channels(ind2)'));
            set(gca,'ytick',1:length(ind2),'yticklabel',num2str(channels(ind2)'));
            x=1;
            
            
        % Calculate correlation pairs using voltage traces
        else
            h=waitbar(0,'Processing channels...','color',palette.background);
            set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
            total_tests = N*(N-1)/2;
            ind=1;
            for i = 1:N-1
                load(mat_file,names{i});
                for j = i+1:N
                    waitbar(ind/total_tests,h);
                    load(mat_file,names{j});
                    eval(sprintf('correlations(i,j) = corr(%s'',%s'');',names{i},names{j}));
                    ind=ind+1;
                end
            end
            delete(h);
        
        end
        
        % *** TEMPORARY FOR LGN PURPOSES ONLY ***
        ind=cell2mat({[11 12 6 5],[13 14 4 3],[15 16 2 1],[25 26 24 23],[27 28 22 21],[29 39 20 19],[31 32 18 17],[10 9 7 8]});
        deletelocs=false(1,length(ind));
        for i = 1:length(ind)
            if(~any(channels==ind(i)))
                deletelocs(i)=true;
            end
        end
        ind(deletelocs)=[];
        
        ind2=zeros(1,length(ind));
        for i = 1:length(ind2)
            ind2(i)=find(channels==ind(i));
        end
        
        fig_std;
        imagesc(correlations(ind2,ind2));
        set(gca,'xtick',1:length(ind2),'xticklabel',num2str(channels(ind2)'));
        set(gca,'ytick',1:length(ind2),'yticklabel',num2str(channels(ind2)'));
        %
        set(handles.status_label,'String','Bloop!');
        delete(custom_panel);
    case 'No'
    otherwise
end

  
  
  

% --- Manual input electrode configuration
function electrode_manualconfigure(~,~,custom_panel,handles)
% Note: this function has not yet been implemented

x=1;


% --- Switches the electrode configuration
function change_electrode(~,~,electrode_menu,panel,configurations,geometry,handles)
global waveforms features data_v;
old_type=getappdata(handles.output,'electrode');
type_list=get(electrode_menu,'String');
old_loc=find(strcmp(type_list,old_type));
new_loc=get(electrode_menu,'Value');
new_type=type_list{new_loc};

if(old_loc==new_loc)
    delete(panel);
else
    current_configuration=configurations{old_loc};
    new_configuration=configurations{new_loc};
    
    old_ch_len=cellfun(@length,current_configuration);
    old_ch=horzcat(current_configuration{:});
    old_shanks=cell(1,length(current_configuration));
    for i = 1:length(old_ch_len)
        old_shanks{i}=ones(1,old_ch_len(i))*i;
    end
    old_shanks=horzcat(old_shanks{:});
    
    new_ch_len=cellfun(@length,new_configuration);
    new_ch=horzcat(new_configuration{:});
    new_shanks=cell(1,length(current_configuration));
    for i = 1:length(new_ch_len)
        new_shanks{i}=ones(1,new_ch_len(i))*i;
    end
    
    %if(~isequal(sort(old_ch),sort(new_ch)))
    %    msgbox('Error: channel configuration is not compatible with this data set.');
    %    return;
    %end
            
    % make sure user realizes that this will wipe data
    choice = questdlg('Changing electrode configuration will wipe all waveform and cluster data. Are you sure you want to continue?',...
        'Change Electrode Configuration','Yes','No','No');

    switch choice
        
        % need to switch channels
        case 'Yes'
            exp=getappdata(handles.output,'exp');
            shank=cell(1,length(new_configuration));
            
            % update configuration, retain threshold when possible
            x=1;
            num_shanks = length(exp.experiment.shank);
            old_info = cell(1,0);
            channels = [];
            ind=1;
            for i = 1:num_shanks
                old_info = [old_info exp.experiment.shank{i}.channel];
                for j = 1:length(exp.experiment.shank{i}.channel)
                    channels(ind) = str2double(exp.experiment.shank{i}.channel{j}.Attributes.num);
                    ind=ind+1;
                end
            end
            
            % grab all shank and channel information
            exp.experiment.shank=cell(1,length(new_configuration));
            for i = 1:length(new_configuration)
                exp.experiment.shank{i}.Attributes.num=i;
                exp.experiment.shank{i}.channel=cell(1,length(new_configuration{i}));
                for j = 1:length(new_configuration{i})
                    ch=new_configuration{i}(j);
                    ind=find(channels==ch);
                    if(isempty(ind))
                        exp.experiment.shank{i}.channel{j}.gain.Text=old_info{1}.gain.Text; % note: this might not be true, but probably is
                        exp.experiment.shank{i}.channel{j}.threshold.Text='5';
                        exp.experiment.shank{i}.channel{j}.threshold2.Text='0';
                        exp.experiment.shank{i}.channel{j}.enabled.Text='true';
                        exp.experiment.shank{i}.channel{j}.Attributes.num=ch;
                    end
                    ind=find_index(channels,new_configuration{i}(j));
                    exp.experiment.shank{i}.channel(j)=old_info(ind);
                end
            end

            exp.experiment.electrode=new_type;
            
            xml_file=getappdata(handles.output,'xml_file');
            struct2xml(exp,xml_file);
            setappdata(handles.output,'exp',exp);
            delete(panel);
            
            waveforms=[]; features=[]; data_v=[];
            setappdata(handles.output,'idx',[]);
            setappdata(handles.output,'timestamps',[]);
            setappdata(handles.output,'t',[]);
            setappdata(handles.output,'electrode_geometry',geometry{new_loc});
            
            export_button_Callback([],[],handles);
            loadfile_item_Callback([],[],handles,sprintf('%s.plx',getappdata(handles.output,'basename')),getappdata(handles.output,'pathname'));
    end
end


% --- Executes on button press in ISI_move_button.
function ISI_move_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to ISI_move_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global palette;
set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Moving spikes with ISI < 1ms to noise cluster...'); drawnow;
ISI_THRESH=1e-3;

current_clusters=getappdata(handles.output,'cluster_number');
idx=getappdata(handles.output,'idx');
setappdata(handles.output,'idx_old',idx);
timestamps=getappdata(handles.output,'timestamps');
m=max(idx)+1;
for i = 1:length(current_clusters)
    locs=find(idx==current_clusters(i));
    ts=timestamps(locs);
    badlocs1=[diff(ts)<ISI_THRESH false];
    badlocs2=[false diff(ts)<ISI_THRESH];
    idx(locs(badlocs1))=m;
    idx(locs(badlocs2))=m+1;
end
setappdata(handles.output,'idx',idx);
fix_clusters(handles);

% update screen
panel_selected=getappdata(handles.output,'panel_selected');
if(strcmp(panel_selected,'features'))
    show_featurespanel(handles);
elseif(strcmp(panel_selected,'cell_stats'))
    display_cell_stats(handles);
end
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);


% --- Executes on button press in move_misaligned_button.
function move_misaligned_button_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
% hObject    handle to move_misaligned_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global waveforms;
%waveforms=getappdata(handles.output,'waveforms');

% note: This isn't the actual channel #, just the index
selected_channels=getappdata(handles.output,'selected_channels');
if(~any(selected_channels)),selected_channels=~selected_channels; end
current_channel=find(selected_channels);

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

% take a few points to the left and right of t
margin=2;
ind=(max(center_loc-margin,1):min(center_loc+margin,length(t)));

waveforms_temp=waveforms(current_channel,cluster_locs,ind);

center_loc=find_index(t(ind),0);
minlocs=zeros(size(waveforms_temp,2),length(current_channel));
maxlocs=zeros(size(waveforms_temp,2),length(current_channel));

% find bad peaks
m=max(idx)+1;
for i = 1:length(current_channel)
    [~,minlocs(:,i)]=min(waveforms_temp(i,:,:),[],3);
    [~,maxlocs(:,i)]=max(waveforms_temp(i,:,:),[],3);
end
badlocs=all(minlocs~=center_loc,2) & all(maxlocs~=center_loc,2);
idx(cluster_locs(badlocs))=m;
setappdata(handles.output,'idx',idx);

% update display
fix_clusters(handles);
show_featurespanel(handles);
set(handles.status_label,'String',sprintf('Moved %g waveforms.\n',sum(badlocs)));



function slope_start_field_Callback(hObject, eventdata, handles)
% hObject    handle to slope_start_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slope_start_field as text
%        str2double(get(hObject,'String')) returns contents of slope_start_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val)), val=50; end
val=round(val);
t_before=str2double(get(handles.t_before_field,'String'));
t_after=str2double(get(handles.t_after_field,'String'));
if(val<-t_before), val=-t_before;
elseif(val>t_after), val=t_after;
end
set(hObject,'string',sprintf('%g',val));
setappdata(handles.output,'slope_edited',true);

% --- Executes during object creation, after setting all properties.
function slope_start_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slope_start_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slope_stop_field_Callback(hObject, eventdata, handles)
% hObject    handle to slope_stop_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slope_stop_field as text
%        str2double(get(hObject,'String')) returns contents of slope_stop_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val)), val=150; end
val=round(val);
t_before=str2double(get(handles.t_before_field,'String'));
t_after=str2double(get(handles.t_after_field,'String'));
if(val<-t_before), val=-t_before;
elseif(val>t_after), val=t_after;
end
set(hObject,'string',sprintf('%g',val));
setappdata(handles.output,'slope_edited',true);

% --- Executes during object creation, after setting all properties.
function slope_stop_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slope_stop_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- extracts features from waveforms
function M_extract_features(t, params, handles, PCA_type)
global features waveforms;
%--------------------------------------------------------------------------
% extract_features.m - Given a set of waveforms, extracts the requested
% features and returns them in a matrix.
%
% Usage: features = extract_features(waveforms);
%
% Input:  waveforms                * CxWxT matrix of waveforms (C channels,
%                                    W waveforms, T time points)
%         t                        * time series t
%         params (optional)        * parameter file (optional)
% Output: features                 * CxWxF matrix of feature projections.
%                                    Features are concatenated into the last
%                                    dimension.
%
% List of possible features:
%     - PCA
%     - initial slope
%     - peak-trough amplitude
%     - peak-trough time interval
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
[num_channels, num_wf, ~] = size(waveforms);
% Deal with parameter file
if(~exist('params','var') || isempty(params))
    params = struct('PCA',3,'slope',0,'pt_amplitude',0,'pt_interval',0);
else
    if(ischar(params))
        params = load_params(params);
    end
end
feature_names=fieldnames(params);
locs = structfun(@(x) x==0,params);
params = rmfield(params,feature_names(locs));
feature_names(locs)=[];
num_features = length(feature_names);
if(any(strcmp('PCA',feature_names)))
    num_PCs=params.PCA;
    feature_length=num_features-1+params.PCA;
else
    num_PCs=0;
    feature_length=num_features;
end

features = zeros(num_channels,num_wf, feature_length);
idx=getappdata(handles.output,'idx');
switch PCA_type
    case 'all'
        selection = true(1,num_wf);
    case 'active'
        selection = idx~=0;
    case 'clusters'
        cluster_number=getappdata(handles.output,'cluster_number');
        u=unique(cluster_number);
        selection=false(1,num_wf);
        for i = 1:length(u)
            selection(idx==u(i))=true;
        end
end

% Extract each type of feature one by one
if(any(strcmp(feature_names,'valleys')) || any(strcmp(feature_names,'peaks')) || any(strcmp(feature_names,'slope')))
    % check 100us away from the valley
    Fs=abs(round(1/(t(2)-t(1))));
    
    x_start=str2double(get(handles.slope_start_field,'String'))*1e-6;
    x_stop=str2double(get(handles.slope_stop_field,'String'))*1e-6;
    dx=x_stop-x_start;
    
    pts_start=round(Fs*x_start); pts_stop=round(Fs*x_stop);
    
    zeroloc=find_index(t,0);    
    valleys=squeeze(waveforms(:,:,zeroloc));
    slope_start=squeeze(waveforms(:,:,zeroloc+pts_start));
    slope_stop=squeeze(waveforms(:,:,zeroloc+pts_stop));
    slopes=(slope_stop-slope_start)./dx;
    
    max_peakspan=300e-6;
    pts_peak=round(Fs*max_peakspan);
    peaks=max(waveforms(:,:,zeroloc:zeroloc+pts_peak),[],3);
    
    energy=sum(waveforms.^2,3);

    nanlocs=isnan(slopes);
    slopes(nanlocs)=0;
end

for c = 1:num_channels
    feature_index=1;
    if(any(strcmp(feature_names,'amplitude')))
        wfs=squeeze(waveforms(c,:,:));
    end
    for f = 1:num_features
        switch(feature_names{f})	% select out the feature
            % run Principal Component Analysis
            case 'PCA'
                % if only active waveforms, we have to run PCA on a subset
                % and then project waveforms onto the PCs
                
                if(num_wf>1)
                    wfs=squeeze(waveforms(c,:,:));
                    % note: this overrides using "active waveforms"
                    if(strcmp(PCA_type,'moving'))
                        timestamps=getappdata(handles.output,'timestamps');
                        window_size=round(str2double(get(handles.PCA_window_field,'String')));
                        [~,scores]=princomp_time(wfs,timestamps,window_size,'pchip');
                    else
                        [coeff,~]=pca(wfs(selection,:));
                        wfs_nomean=wfs-repmat(mean(wfs,1),size(wfs,1),1);
                        scores=wfs_nomean*coeff;
                    end
                else
                    scores=zeros(num_wf,params.PCA);
                end
                
                features(c,:,feature_index+(0:num_PCs-1))=scores(:,1:num_PCs);
                feature_index=feature_index+num_PCs;
            % Maximum point
            case 'peaks'
                features(c,:,feature_index)=peaks(c,:);
                feature_index=feature_index+1;
            % Minimum point
            case 'valleys'
                features(c,:,feature_index)=valleys(c,:);
                feature_index=feature_index+1;
            % Slope from peak to trough
            case 'slope'
                features(c,:,feature_index)=slopes(c,:);
                feature_index=feature_index+1;
            % integral of squared waveform
            case 'energy'
                features(c,:,feature_index)=energy(c,:);
                feature_index=feature_index+1;
            % maximum - minimum span
            case 'amplitude'
                w=squeeze(waveforms(c,:,:));
                features(c,:,feature_index)=max(wfs,[],2)-min(wfs,[],2);
        end
    end
end


% --------------------------------------------------------------------
function PCA_menu_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function PCA_menu_calculate_PCs_active_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu_calculate_PCs_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calculate_features(handles,'active');
plot_features(handles);

% --------------------------------------------------------------------
function PCA_menu_calculate_PCs_all_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu_calculate_PCs_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calculate_features(handles,'all');
plot_features(handles);

% --------------------------------------------------------------------
function PCA_menu_calculate_PCs_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu_calculate_PCs_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calculate_features(handles,'clusters');
plot_features(handles);


% --------------------------------------------------------------------
function PCA_menu_calculate_PCs_moving_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu_calculate_PCs_moving (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
calculate_features(handles,'moving');
plot_features(handles);


% --------------------------------------------------------------------
function PCA_menu_remove_outliers_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_menu_remove_outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
remove_outliers(handles);



function same_ch_deadzone_field_Callback(hObject, eventdata, handles)
% hObject    handle to same_ch_deadzone_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of same_ch_deadzone_field as text
%        str2double(get(hObject,'String')) returns contents of same_ch_deadzone_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val < 0), val=1000; end
set(hObject,'String',sprintf('%g',round(val)));


% --- Executes during object creation, after setting all properties.
function same_ch_deadzone_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to same_ch_deadzone_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function diff_ch_deadzone_field_Callback(hObject, eventdata, handles)
% hObject    handle to diff_ch_deadzone_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diff_ch_deadzone_field as text
%        str2double(get(hObject,'String')) returns contents of diff_ch_deadzone_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val < 0), val=50; end
set(hObject,'String',sprintf('%g',round(val)));

% --- Executes during object creation, after setting all properties.
function diff_ch_deadzone_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diff_ch_deadzone_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fastsort_field_Callback(hObject, eventdata, handles)
% hObject    handle to fastsort_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fastsort_field as text
%        str2double(get(hObject,'String')) returns contents of fastsort_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val<0), val=1000; end
set(hObject,'String',sprintf('%g',val));

% --- Executes during object creation, after setting all properties.
function fastsort_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fastsort_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fastsort_checkbox.
function fastsort_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to fastsort_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fastsort_checkbox
if(get(hObject,'Value'))
    set(handles.fastsort_field,'Enable','on');
else
    set(handles.fastsort_field,'Enable','off');
end


% --- Executes on button press in channel_sort_button.
function channel_sort_button_Callback(hObject, eventdata, handles)
% hObject    handle to channel_sort_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global palette;
trigger_ch=getappdata(handles.output,'trigger_ch');
if(~isempty(trigger_ch))
    setappdata(handles.output,'idx_old',getappdata(handles.output,'idx'));
    setappdata(handles.output,'idx',trigger_ch);
    fix_clusters(handles);
    setappdata(handles.output,'cluster_number',1);
    show_featurespanel(handles);
else
    set(handles.status_label,'String','Trigger channels not set. Re-threshold to generate.');
    set(handles.status_light,'backgroundcolor',palette.orange);
end





function outliers_field_Callback(hObject, eventdata, handles)
% hObject    handle to outliers_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outliers_field as text
%        str2double(get(hObject,'String')) returns contents of outliers_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val<0), val=10; end
set(hObject,'String',sprintf('%g',val));

% --- Executes during object creation, after setting all properties.
function outliers_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outliers_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
% realign all waveforms in the selected cluster to the currently selected
% channel. If multiple channels are selected, pick the first one, since
% this shouldn't happen.
function edit_menu_realign_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_realign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% nope: this function is currently defunct.
return;

cluster_number=getappdata(handles.output,'cluster_number');
selected_channel=find(getappdata(handles.output,'selected_channels'),1,'first');
idx=getappdata(handles.output,'idx');
timestamps=getappdata(handles.output,'timestamps');
t=getappdata(handles.output,'t');

locs=false(1,length(idx));
for i = 1:length(cluster_number)
    locs(idx==selected_channel)=true;
end

trigger_ch=ones(1,sum(locs))*selected_channel;
% look up to 50us in each direction
Fs=getappdata(handles.output,'Fs');
search_width = round(50e-6*Fs);
if(get(handles.interp_waveforms_checkbox,'Value'))
    search_width=search_width*4;
end
[waveforms2,timestamps2,t2]=align_waveforms(waveforms(:,locs,:),timestamps(locs),t,search_width,trigger_ch);

zero_pt=find_index(t,0);
zero_pt_new=find_index(t2,0);
margin_start=zero_pt-zero_pt_new;
margin_end=(length(t)-zero_pt)-(length(t2)-zero_pt_new);

waveforms=waveforms(:,:,margin_start:end-margin_end+1);
waveforms(:,locs,:)=waveforms2;
timestamps(locs)=timestamps2;
t=t2;
setappdata(handles.output,'timestamps',timestamps);
setappdata(handles.output,'t',t);



function interp_factor_field_Callback(hObject, eventdata, handles)
% hObject    handle to interp_factor_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interp_factor_field as text
%        str2double(get(hObject,'String')) returns contents of interp_factor_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val<1 || mod(val,1)~=0), val = 4; end
set(hObject,'String',sprintf('%g',val));

% --- Executes during object creation, after setting all properties.
function interp_factor_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interp_factor_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PCA_window_field_Callback(hObject, eventdata, handles)
% hObject    handle to PCA_window_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PCA_window_field as text
%        str2double(get(hObject,'String')) returns contents of PCA_window_field as a double
val=str2double(get(hObject,'String'));
if(~isfinite(val) || val<0), val=60; end
set(hObject,'String',sprintf('%g',round(val)));

% --- Executes during object creation, after setting all properties.
function PCA_window_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PCA_window_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function fig_msort_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to fig_msort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fclose('all');


% --- Executes on button press in delete_junk_button.
function delete_junk_button_Callback(hObject, eventdata, handles)
% hObject    handle to delete_junk_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global features waveforms;
idx=getappdata(handles.output,'idx');
bad_locs=idx==0;

% update features
features(bad_locs,:)=[];

% update interpolated data
w=getappdata(handles.output,'waveforms_interp');
w(:,bad_locs,:)=[];
setappdata(handles.output,'waveforms_interp',w);
w=getappdata(handles.output,'normalized_waveforms_interp');
w(:,bad_locs,:)=[];
setappdata(handles.output,'normalized_waveforms_interp',w);
t=getappdata(handles.output,'timestamps_interp');
t(bad_locs)=[];
setappdata(handles.output,'timestamps_interp',t);

% update regular data
w=getappdata(handles.output,'original_waveforms');
w(:,bad_locs,:)=[];
setappdata(handles.output,'original_waveforms',w); clear w;
w=getappdata(handles.output,'original_normalized_waveforms');
w(:,bad_locs,:)=[];
setappdata(handles.output,'original_normalized_waveforms',w); clear w;
t=getappdata(handles.output,'original_timestamps');
t(bad_locs)=[];
setappdata(handles.output,'original_timestamps',t); clear t;

% update current data
waveforms(:,bad_locs,:)=[];
w=getappdata(handles.output,'normalized_waveforms');
w(:,bad_locs,:)=[];
setappdata(handles.output,'normalized_waveforms',w);
t=getappdata(handles.output,'timestamps');
t(bad_locs)=[];
setappdata(handles.output,'timestamps',t);

% update trigger channels
trigger_ch=getappdata(handles.output,'trigger_ch');
trigger_ch(bad_locs)=[];
setappdata(handles.output,'trigger_ch',trigger_ch);

% update idx
%idx=getappdata(handles.output,'idx');
idx_old=getappdata(handles.output,'idx_old');
if(isempty(idx_old)), idx_old=zeros(1,length(idx),'uint8'); end
idx(bad_locs)=[];
idx_old(bad_locs)=[];

setappdata(handles.output,'idx',idx);
setappdata(handles.output,'idx_old',idx_old);

features_button_Callback([],[],handles);
set(handles.status_label,'String',sprintf('%g waveforms removed.',length(bad_locs)));



function default_data_field_Callback(hObject, eventdata, handles)
% hObject    handle to default_data_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of default_data_field as text
%        str2double(get(hObject,'String')) returns contents of default_data_field as a double


% --- Executes during object creation, after setting all properties.
function default_data_field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to default_data_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in default_data_button.
function default_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to default_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir;
if(~isempty(folder_name))
    set(handles.default_data_field,'String',folder_name);
    setappdata(handles.output,'defaultdirchanged',true);
end


% --------------------------------------------------------------------
function edit_menu_screenshot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_screenshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eventData.Key = 'p';
key_press([],eventData,handles);


% --- Executes on button press in screenshot_bw_checkbox.
function screenshot_bw_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to screenshot_bw_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of screenshot_bw_checkbox


% --- Executes on button press in waveforms_left_button.
function waveforms_left_button_Callback(hObject, eventdata, handles)
% hObject    handle to waveforms_left_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slide_waveforms(handles,-1);

% --- Executes on button press in waveforms_right_button.
function waveforms_right_button_Callback(hObject, eventdata, handles)
% hObject    handle to waveforms_right_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
slide_waveforms(handles,1);


% --- Moves waveforms left or right (left = -1, right = 1)
function slide_waveforms(handles,direction)
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

% --- Grabs waveforms from voltage truce using selected timestamps
% note: if "offsets" is provided, then a small series of offsets is added
% to some of the waveforms. 
function capture_waveforms(handles)
global data_v Ts waveforms;

interp_factor=round(str2double(get(handles.interp_factor_field,'String')));

timestamps=getappdata(handles.output,'timestamps');
if(isempty(timestamps))
    waveforms=[];
    t=[]; %#ok<NASGU>
    return;
end
trigger_ch=getappdata(handles.output,'trigger_ch');
t_before=getappdata(handles.output,'t_before')*1e-6;
t_after=getappdata(handles.output,'t_after')*1e-6;
Fs=getappdata(handles.output,'Fs');
pts_before=round(Fs*t_before)+1;
pts_after=round(Fs*t_after)+1;

% convert timestamps to Ts locations
timestamp_locs = int32(round((timestamps-Ts(1))*Fs)+1);
badlocs=timestamp_locs<pts_before | timestamp_locs>length(Ts)-pts_after;
timestamp_locs(badlocs)=[]; timestamps(badlocs)=[]; trigger_ch(badlocs)=[];
N=length(timestamp_locs);
C = size(data_v,2);
L=pts_before+pts_after+1;
adder=repmat(int32(-pts_before:pts_after)',1,N);

index = repmat(timestamp_locs,L,1)+adder;
index = index(:);
waveforms = permute(reshape(data_v(index,:),[L N C]),[3 2 1]);

% interpolate waveforms
t=(-pts_before:pts_after)/Fs;
[w,t_interp]=interp_waveforms(waveforms,t,interp_factor);
    
% note: doing this may push two waveforms together.
[waveforms_interp,timestamps_interp,t_interp]=align_waveforms(w,timestamps,t_interp,interp_factor-1,trigger_ch);
badlocs=t_interp<-t_before | t_interp>t_after;
t_interp(badlocs)=[]; waveforms_interp(:,:,badlocs)=[];
normalized_waveforms=normalize_waveforms(waveforms);
normalized_waveforms_interp=normalize_waveforms(waveforms_interp);
waveforms=waveforms(:,:,2:end-1); t=t(2:end-1);
normalized_waveforms=normalized_waveforms(:,:,2:end-1);

% update application data
setappdata(handles.output,'original_waveforms',waveforms);
setappdata(handles.output,'original_timestamps',timestamps);
setappdata(handles.output,'original_normalized_waveforms',normalized_waveforms);
setappdata(handles.output,'original_t',t);
setappdata(handles.output,'waveforms_interp',waveforms_interp);
setappdata(handles.output,'timestamps_interp',timestamps_interp);
setappdata(handles.output,'t_interp',t_interp);
setappdata(handles.output,'normalized_waveforms_interp',normalized_waveforms_interp);

% set current data based on check box
if(get(handles.interp_waveforms_checkbox,'Value'))
    waveforms=waveforms_interp;
    normalized_waveforms=normalized_waveforms_interp;
    t=t_interp;
    timestamps=timestamps_interp;
end
setappdata(handles.output,'waveforms',waveforms);
setappdata(handles.output,'normalized_waveforms',normalized_waveforms);
setappdata(handles.output,'timestamps',timestamps);
setappdata(handles.output,'t',t);

% --------------------------------------------------------------------
function edit_menu_show_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_show_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eventData.Key = 'h';
key_press([],eventData,handles);


% --------------------------------------------------------------------
function edit_menu_subtract_common_mean_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu_subtract_common_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data_v palette;

set(handles.status_light,'BackgroundColor',palette.orange);
set(handles.status_label,'String','Subtracting mean...');
drawnow;

data_mean=mean(data_v,2);
data_v=data_v-repmat(data_mean,1,size(data_v,2));
data_mean=mean(data_v,1);
data_v=data_v-repmat(data_mean,size(data_v,1),1);
s=std(data_v,1);
data_v=data_v./repmat(s,size(data_v,1),1);

set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);
drawnow;