function varargout = voltage_trace(varargin)
% VOLTAGE_TRACE MATLAB code for voltage_trace.fig
%      VOLTAGE_TRACE, by itself, creates a new VOLTAGE_TRACE or raises the existing
%      singleton*.
%
%      H = VOLTAGE_TRACE returns the handle to a new VOLTAGE_TRACE or the handle to
%      the existing singleton*.
%
%      VOLTAGE_TRACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VOLTAGE_TRACE.M with the given input arguments.
%
%      VOLTAGE_TRACE('Property','Value',...) creates a new VOLTAGE_TRACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before voltage_trace_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to voltage_trace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help voltage_trace

% Last Modified by GUIDE v2.5 18-Jan-2013 14:27:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @voltage_trace_OpeningFcn, ...
                   'gui_OutputFcn',  @voltage_trace_OutputFcn, ...
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


% --- Executes just before voltage_trace is made visible.
function voltage_trace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to voltage_trace (see VARARGIN)

% Choose default command line output for voltage_trace
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes voltage_trace wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = voltage_trace_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadfile_button.
function loadfile_button_Callback(hObject, eventdata, handles,reload)
% hObject    handle to loadfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% add Plexon folders to path
Plexon_folder = sprintf('%s\\Plexon',pwd);
mexPlex_folder = sprintf('%s\\Plexon\\mexPlex',pwd);
addpath(Plexon_folder,mexPlex_folder);

if(~exist('reload','var') || ~reload)
    [filename,pathname] = uigetfile('*.plx','Select the .plx file.');
    
    % parse filename
    expr = regexp(filename,'(?<basename>.*)\.plx','names');
    basename = expr.basename;
    setappdata(handles.output,'basename',basename);
    setappdata(handles.output,'pathname',pathname);
else
    pathname=getappdata(handles.output,'pathname');
    basename=getappdata(handles.output,'basename');
end

%filename='depth=23642_fast_spike.plx';
%pathname='E:\Data\2013-01-07-Monkey\Pen 1\';


% look for PLX file
plx_file = sprintf('%s%s.plx',pathname,basename);
setappdata(handles.output,'plx_file',plx_file);

% look for XML file
xml_file = sprintf('%s%s.xml',pathname,basename);
if(~exist(xml_file,'file'))
    % generate XML file
    %set(handles.status_label,'String','Generating associated .xml file...');
    generate_standard_Plexon_xml(plx_file);
end
setappdata(handles.output,'xml_file',xml_file);

% look for MAT file
mat_exists=false;
mat_file=sprintf('%s%s.mat',pathname,basename);
if(~exist(mat_file,'file'))
    save(mat_file,'');
else
    mat_exists=true;
end
setappdata(handles.output,'mat_file',mat_file);

% load struct into memory
exp = xml2struct(xml_file);
setappdata(handles.output,'exp',exp);

% Populate Information fields
num_shanks = length(exp.experiment.shank);
channel_shanks = cell(1,num_shanks);
shank_nums = zeros(1,num_shanks);
shank_thresholds = cell(1,num_shanks);
disabled=cell(1,num_shanks);
for i = 1:num_shanks
    % get shank numbers
    shank_nums(i) = str2double(exp.experiment.shank{i}.Attributes.num);
    num_channels = length(exp.experiment.shank{i}.channel);
    thresh = zeros(1,num_channels);
    disabled{i}=false(1,num_channels);
    
    % Note: this special case is due to the fact that xml2struct doesn't
    % place single items into cell arrays
    if(num_channels==1)
        channel_shanks{i} = str2double(exp.experiment.shank{i}.channel.Attributes.num);
        thresh = str2double(exp.experiment.shank{i}.channel.threshold.Text);
        if(strcmp(exp.experiment.shank{i}.channel.enabled.Text,'false'))
            disabled{i}=true;
        end
    else
        for c = 1:num_channels
            channel_shanks{i}(c) = str2double(exp.experiment.shank{i}.channel{c}.Attributes.num);
            thresh(c) = str2double(exp.experiment.shank{i}.channel{c}.threshold.Text);
            if(strcmp(exp.experiment.shank{i}.channel{c}.enabled.Text,'false'))
                disabled{i}(c)=true;
            end
        end
    end
    thresh(disabled{i})=[];
    channel_shanks{i}(disabled{i})=[];
    shank_thresholds{i} = thresh;
end

Fs = str2double(exp.experiment.ADRate.Text);
num_total_channels = sum(cell2mat(cellfun(@length,channel_shanks,'UniformOutput',false)));

% duration of experiment
exp_duration = str2double(exp.experiment.duration.Text);
setappdata(handles.output,'duration',exp_duration);

% save data to application
setappdata(handles.output,'num_shanks',num_shanks);
setappdata(handles.output,'channel_shanks',channel_shanks);
setappdata(handles.output,'shank_nums',shank_nums);
setappdata(handles.output,'shank_thresholds',shank_thresholds);
setappdata(handles.output,'Fs',Fs);
setappdata(handles.output,'exp_duration',exp_duration);
setappdata(handles.output,'num_total_channels',num_total_channels);

% update Information fields
set(handles.filename_label,'String',basename);
set(handles.info_shanks_label,'String',num2str(num_shanks));
set(handles.info_channels_label,'String',num2str(num_total_channels));
set(handles.info_duration_label,'String',time2str(exp_duration));
set(handles.info_Fs_label,'String',sprintf('%g kHz',Fs/1e3));

% Create shank buttons near bottom of page
axis(handles.ax,'manual'); grid on; hold all;
ax_pos=get(handles.ax,'position');
ax_width=ax_pos(3);

% center buttons underneath axis
margin_x=0.03; margin_y=.03; margin_left=ax_pos(1);
button_width=(1-margin_x*(num_shanks+1))/num_shanks;
remaining_height=ax_pos(2)-.03; % leave .01 top and bottom
button_height=1-margin_y;

% generate buttons
shankbutton_panel=uibuttongroup('visible','on','bordertype','none',...
    'units','normalized',...
    'position',[margin_left .01 ax_width remaining_height]);
shank_buttons=zeros(1,num_shanks);
for i = 1:num_shanks
    x=i*margin_x+(i-1)*button_width;
    shank_buttons(i)=uicontrol(shankbutton_panel,'Style','togglebutton',...
        'units','normalized','String',sprintf('Shank %g',i),...
        'fontsize',14,'visible','off','enable','off',...
        'Position',[x margin_y button_width button_height]);
end
set(shank_buttons,'value',0);
set(shankbutton_panel,'SelectionChangeFcn',{@toggle_button_callback,handles});

if(~mat_exists)
    set(handles.convert_button,'enable','on','visible','on');
else
    set(handles.convert_button,'enable','off','visible','off');
    set(shankbutton_panel,'visible','on');
    set(shank_buttons,'visible','on','enable','on');
end

setappdata(handles.output,'shankbutton_panel',shankbutton_panel);
setappdata(handles.output,'shank_buttons',shank_buttons);

% set the x_slider bar so that we don't jump too far with each click
% (should be 1/10th of x window)
window_time=get(handles.horz_slider,'Value');
duration=getappdata(handles.output,'duration');
shift_amount=(window_time*.1)/duration;
set(handles.horz_slider,'SliderStep',[shift_amount shift_amount*10]);


% --- Executes on button press in convert_button.
function convert_button_Callback(hObject, eventdata, handles)
% hObject    handle to convert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plx_file=getappdata(handles.output,'plx_file');
plx2mat(plx_file);
set(hObject,'visible','off','enable','off');

% reload everything
reload=true;
loadfile_button_Callback(handles.output,[],handles,reload);


% --- Executes on button press of toggle buttons
function toggle_button_callback(hObject,eventdata,handles)
global window_size data voltage_plot Ts;
cla(handles.ax);
set(handles.loading_label,'visible','on'); drawnow;
ch_tags=getappdata(handles.output,'ch_tags');
if(any(ishandle(ch_tags))), delete(ch_tags); end
buttons=getappdata(handles.output,'shank_buttons');
selected_shank=find(buttons==eventdata.NewValue);

% load up shank
channel_shanks=getappdata(handles.output,'channel_shanks');
channels=channel_shanks{selected_shank}; %#ok<FNDSB>

% temporary, 1000s of data only
duration=getappdata(handles.output,'duration');
mat_file=getappdata(handles.output,'mat_file');
Fs=getappdata(handles.output,'Fs');
load(mat_file,'Ts'); num_pts=length(Ts);
Ts=Ts-Ts(1);
num_channels=length(channels);
data=zeros(num_pts,num_channels,'single');

for i = 1:length(channels)
    varname=sprintf('data_%g',channels(i));
    load(mat_file,varname);
    eval(sprintf('data(:,i)=%s;',varname));
    clear(varname);
end

% note: data has been normalized to std dev of +/-1
setappdata(handles.output,'data',data);
setappdata(handles.output,'Ts',Ts);

% adjust each row of data so that it occupies a different range
std_dev_zoom=get(handles.vert_slider,'Value');
offset=(0:num_channels-1)*std_dev_zoom*2;
offset=repmat(offset,num_pts,1);
data=data+offset;
setappdata(handles.output,'std_dev_zoom',std_dev_zoom);

% generate plot range
window_time=get(handles.horz_slider,'Value'); % put in something to adjust this later?
window_size=round(window_time*Fs);
start_ind=1; stop_ind=start_ind+window_size;

set(handles.ax,'xlim',[Ts(start_ind) Ts(stop_ind)],'ylim',[-std_dev_zoom offset(end)+std_dev_zoom])
gray=[.2 .2 .2];
%colors=distinguishable_colors(num_channels,'gray');
colors=get_colors(num_channels); colors=1-((1-colors)*.5);
set(handles.ax,'ColorOrder',colors); hold all;
voltage_plot=plot(handles.ax,Ts(start_ind:stop_ind),data(start_ind:stop_ind,:),'linesmooth','on'); grid on;
set(handles.ax,'color',gray);

% set up channel labels to the left of the plot
ax_pos=get(handles.ax,'position');
channel_height=ax_pos(4)/num_channels;
channel_middle=(0:num_channels-1)*channel_height+ax_pos(2)+channel_height/2;
tag_halfheight=.015;
channel_bot=channel_middle-tag_halfheight;
ch_tags=zeros(1,num_channels);
ch_width=.02;
for i = 1:num_channels
    ch_tags(i)=uicontrol('Style','text','String',sprintf('%g',channels(i)),...
        'fontsize',16,'fontweight','bold','horizontalalignment','right',...'
        'units','normalized',...
        'position',[0 channel_bot(i) ch_width tag_halfheight*2],...
        'parent',handles.output);
end
setappdata(handles.output,'ch_tags',ch_tags);
set(handles.plot_slider,'visible','on','enable','on');

% tie plot to slider
set(handles.plot_slider,'Min',1,'Max',num_pts-window_size,'value',1);
set([handles.horz_slider handles.vert_slider],'enable','on');

set(handles.loading_label,'visible','off'); drawnow;


% --- Executes on slider movement.
function plot_slider_Callback(hObject, eventdata, handles)
% hObject    handle to plot_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global voltage_plot window_size data Ts;
newval=get(hObject,'Value');
start_ind=round(newval); stop_ind=round(newval+window_size);
Xdata=Ts(start_ind:stop_ind); Ydata=data(start_ind:stop_ind,:);
num_channels=size(data,2);
for i = 1:num_channels
    set(voltage_plot(i),'XData',Xdata,'YData',Ydata(:,i));
end
xtick=round(linspace(Xdata(1),Xdata(end),10)*10)/10;
if(length(xtick)>10)
    skip=round(length(xtick)/10);
    xtick=xtick(1:skip:end);
end
if(any(diff(xtick)<=0)),xtick(diff(xtick)<=0)=[]; end
set(handles.ax,'xlim',[Xdata(1) Xdata(end)],'xtick',xtick);

% --- Executes during object creation, after setting all properties.
function plot_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function horz_slider_Callback(hObject, eventdata, handles)
% hObject    handle to horz_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global window_size;
window_time=get(handles.horz_slider,'Value');
Fs=getappdata(handles.output,'Fs');
window_size=round(window_time*Fs);
plot_slider_Callback(handles.plot_slider,[],handles);

% set the x_slider bar so that we don't jump too far with each click
% (should be 1/10th of x window)
window_time=get(handles.horz_slider,'Value');
duration=getappdata(handles.output,'duration');
shift_amount=(window_time*.1)/duration;
set(handles.plot_slider,'SliderStep',[shift_amount shift_amount*10]);

set(handles.window_label,'String',sprintf('%2.1fs',window_time));


% --- Executes during object creation, after setting all properties.
function horz_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to horz_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function vert_slider_Callback(hObject, eventdata, handles)
% hObject    handle to vert_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global data;
old_zoom=getappdata(handles.output,'std_dev_zoom');
new_zoom=get(hObject,'Value');

[num_pts,num_channels]=size(data);
offset=repmat(0:num_channels-1,num_pts,1);
offset=offset*2*(new_zoom-old_zoom);
data=data+offset;
plot_slider_Callback(handles.plot_slider,[],handles);
setappdata(handles.output,'std_dev_zoom',new_zoom);
set(handles.ax,'ylim',[-new_zoom num_channels*2*new_zoom-new_zoom]);

set(handles.std_dev_label,'String',sprintf('%2.0f std dev',new_zoom));


% --- Executes during object creation, after setting all properties.
function vert_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vert_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
