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
if(get(handles.link_channels_checkbox,'Value'))
    set(row,'Value',get(hObject,'Value'));
end


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

function loadfile_item_Callback(~, ~, handles,filename,pathname)
%LOADFILE_ITEM_CALLBACK   Load .plx or .mat file
%   LOADFILE_ITEM_CALLBACK Executes upon the user pressing the Load file button. This loads all metadata from a .mat
%   file. If the associated .mat and .xml files are not available, auto-generates .xml and .mat files and populates
%   appropriate fields.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
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

 
  
  

% --- Manual input electrode configuration
function electrode_manualconfigure(~,~,custom_panel,handles)
% Note: this function has not yet been implemented



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


% --- Grabs waveforms from voltage truce using selected timestamps
% note: if "offsets" is provided, then a small series of offsets is added
% to some of the waveforms. 
function capture_waveforms(handles)


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