function varargout = Mmerge(varargin)
% MMERGE MATLAB code for Mmerge.fig
%      MMERGE, by itself, creates a new MMERGE or raises the existing
%      singleton*.
%
%      H = MMERGE returns the handle to a new MMERGE or the handle to
%      the existing singleton*.
%
%      MMERGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MMERGE.M with the given input arguments.
%
%      MMERGE('Property','Value',...) creates a new MMERGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mmerge_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mmerge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Mmerge

% Last Modified by GUIDE v2.5 19-Nov-2012 14:28:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mmerge_OpeningFcn, ...
                   'gui_OutputFcn',  @Mmerge_OutputFcn, ...
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


% --- Executes just before Mmerge is made visible.
function Mmerge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mmerge (see VARARGIN)

% Choose default command line output for Mmerge
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Mmerge wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Mmerge_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function merge_title_Callback(hObject, eventdata, handles)
% hObject    handle to merge_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of merge_title as text
%        str2double(get(hObject,'String')) returns contents of merge_title as a double


% --- Executes during object creation, after setting all properties.
function merge_title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to merge_title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadfile_button.
function loadfile_button_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global plot_width plot_height ax_centers ax_pos margin_y cell_positions a;

[filename,pathname] = uigetfile('*.plx','Select the .plx file.');
%filename='Pen2_depth=1910um_right_eye_natural_information_stimulus_run1_spl_t001.plx';
%pathname='C:\Data\2012-02-03-Monkey\';

if(~filename)% %#ok<BDSCI,BDLGI>
    set(handles.status_label,'String','Please load a .plx file');
    set(handles.status_light,'BackgroundColor','green');
    return
end

expr = regexp(filename,'(?<basename>.*)\.plx','names');
basename = expr.basename;
expr=regexp(basename,'(?<number>\d*$)','names');
number=expr.number;
num_digits=length(number);
basename=basename(1:end-num_digits);

% Grab all files in this sequence
ls=dir(pathname);
names={ls.name};
locs=regexp(names,sprintf('^%s(\\d)+(\\.mat$)',basename));
locs=~cellfun(@isempty,locs);
names=names(locs);
num_files=length(names);

% TEMPORARY
%num_files=1;

suffix=num2str(1:num_files,sprintf('%%0%g.0f',num_digits));
suffix=reshape(suffix',num_digits,[])';
mat_files=num2cell([repmat(pathname,num_files,1) repmat(basename,num_files,1) suffix repmat('.mat',num_files,1)],2);

% save application data
setappdata(handles.output,'basename',basename);
setappdata(handles.output,'files',basename);
setappdata(handles.output,'pathname',pathname);
setappdata(handles.output,'mat_files',mat_files);

% hide stuff
set([handles.text_instructions handles.text_example handles.text_example_list handles.loadfile_button],'Visible','off');
set(handles.loadfile_button,'Enable','off');

% Load all waveforms, calculate mean waveform templates.
set(handles.status_light,'backgroundcolor','red');
set(handles.status_label,'String','Loading shank data. Please wait...');
drawnow;

% we need # shanks, # cells per shank per file
num_shanks=zeros(1,num_files);
for i = 1:num_files
    % generate file of waveforms
    files=whos('-file',mat_files{i});
    wf_list={files.name};
    wf_list=wf_list(~cellfun(@isempty,regexp(wf_list,'^waveforms')));
    
    num_shanks(i)=length(wf_list);
end
templates=cell(num_files,num_shanks);
num_cells=zeros(num_files,max(num_shanks));
for i = 1:num_files
    % generate file of waveforms
    files=whos('-file',mat_files{i});
    wf_list={files.name};
    wf_list=wf_list(~cellfun(@isempty,regexp(wf_list,'^waveforms')));
    
    % load index
    idx_list={files.name};
    idx_list=idx_list(~cellfun(@isempty,regexp(idx_list,'^idx')));
        
    % run through each shank
    for j = 1:length(wf_list)

        load(mat_files{i},wf_list{j},idx_list{j});
        
        % determine number of cells
        eval(sprintf('u=unique(%s);',idx_list{j}));
        
        eval(sprintf('waveforms=get_mean_waveforms(%s,%s,false);',wf_list{j},idx_list{j}));

        % normalize waveforms
        max_waveform=max(waveforms(:)); min_waveform=min(waveforms(:));
        span=max_waveform-min_waveform;
        if(span~=0)
            waveforms=(waveforms-min_waveform)/span;
        end

        templates{i,j}=waveforms;
        num_cells(i,j)=size(waveforms,1);
    end    
end
setappdata(handles.output,'num_cells',num_cells);
setappdata(handles.output,'num_files',num_files);
setappdata(handles.output,'templates',templates);

% create buttons up top
num_rows=max(num_cells); num_cols=num_files;
margin_top=0.07; margin_left=0.05;
margin_x=.002; margin_y=.02;

button_margin_x=.01;
num_buttons=max(num_shanks);
button_width=(1-(num_buttons+1)*button_margin_x)/num_buttons;

shank_buttongroup=uibuttongroup('units','normalized','Title','Shank','BackgroundColor',[0 0 0],'Parent',handles.main_panel,...
    'selectionhighlight','off','position',[0 1-margin_top+.01 1 margin_top-.01],'SelectionChangeFcn',{@select_shank,handles});
shank_buttons=zeros(1,max(num_shanks));
for i = 1:max(num_shanks)
    % determine x vals
    x=button_margin_x*i+(i-1)*button_width;
    shank_buttons(i)=uicontrol('Style','Togglebutton','units','normalized','String',sprintf('Shank %g',i),...
            'Parent',shank_buttongroup,'selectionhighlight','off','Position',[x .1 button_width .9]);
end
setappdata(handles.output,'shank_buttongroup',shank_buttongroup);
setappdata(handles.output,'shank_buttons',shank_buttons);

saved_cellpositions=cell(1,num_shanks);
max_cells=max(num_cells);
for i = 1:num_shanks
    saved_cellpositions{i}=repmat((1:max_cells(i))',1,num_files);
end
setappdata(handles.output,'saved_cellpositions',saved_cellpositions);
setappdata(handles.output,'shank_num',1);
x=[];
x.OldValue=1; x.NewValue=shank_buttons(1);
cell_positions=saved_cellpositions{1};
a=[];
select_shank([],x,handles);

set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor','green');

% --- Selects an axis upon clicking
% ax: selected axis
function select_axes(hObject,evendata,handles,ax)
global clickpos initial_ax_pos plot_height current_row a current_col cell_positions;
setappdata(handles.output,'selected_axes',ax);
set(handles.output,'WindowButtonMotionFcn',{@move_axis,handles,ax});
uistack(ax,'top');
[current_row,current_col]=find(ax==a);

% determine where the user initially clicked
pos=get(ax,'position');
initial_ax_pos=pos(2);
[~,y]=getcurpt(ax);
y=pos(2)+y*plot_height;
clickpos=y;

% set pointer
set(handles.output,'pointer','fleur');


% --- Deselects an axis upon clicking
% ax: selected axis
function deselect_axes(hObject,evendata,handles)
global current_row current_col ax_pos a;
setappdata(handles.output,'selected_axes',[]);
set(handles.output,'WindowButtonMotionFcn',[]);

% snap to grid
ax=a(current_row,current_col);
set(ax,'position',ax_pos(current_row,current_col,:));

% set pointer
set(handles.output,'pointer','arrow');

% --- Moves an axis after clicking and dragging
% ax: selected axis
function move_axis(hObject,eventdata,handles,ax)
global clickpos initial_ax_pos plot_width plot_height ax_centers a current_row current_col ax_pos margin_y cell_positions;
% calculate actual current mouse position
pos=get(ax,'position');
[~,y]=getcurpt(ax);

mouseY=pos(2)+y*plot_height;

% determine how far mouse has moved
y_offset=mouseY-clickpos;

% determine new axis position
newpos=initial_ax_pos+y_offset;
if(newpos<0),newpos=0; end
if(newpos>1-plot_height),newpos=1-plot_height; end

% determine closest center
center=newpos+plot_height/2;
dist=center-ax_centers(:,current_col,2);
[~,closest_row]=min(abs(dist));

% need to move the other axis
if(closest_row~=current_row)
    % move axis to adjacent location
    closest_ax=a(closest_row,current_col);
    pos=get(closest_ax,'position');
    if(closest_row>current_row)
        % move axis up one
        change=1;
    elseif(closest_row<current_row)
        % move axis down one
        change=-1;
    end
    
    set(closest_ax,'position',pos+[0 plot_height+margin_y 0 0]*change);
    
    % swap axis positions
    tmp=ax_pos(closest_row,current_col);
    ax_pos(closest_row,current_col)=ax_pos(current_row,current_col);
    ax_pos(current_row,current_col)=tmp;
    
    tmp=a(closest_row,current_col);
    a(closest_row,current_col)=a(current_row,current_col);
    a(current_row,current_col)=tmp;

    tmp=cell_positions(closest_row,current_col);
    cell_positions(closest_row,current_col)=cell_positions(current_row,current_col);
    cell_positions(current_row,current_col)=tmp;
    
    current_row=current_row+change;
end

set(ax,'position',[pos(1) newpos plot_width plot_height]); drawnow;

    
% --- Switch shanks when user presses shank button
function select_shank(hObject,eventdata,handles)
global cell_positions a ax_pos ax_centers plot_height plot_width;
num_cells=getappdata(handles.output,'num_cells');
templates=getappdata(handles.output,'templates');

shank_buttons=getappdata(handles.output,'shank_buttons');

% determine which shank number we are
old_shanknum=getappdata(handles.output,'shank_num');

% save cellposition
saved_cellpositions=getappdata(handles.output,'saved_cellpositions');
saved_cellpositions{old_shanknum}=cell_positions;
setappdata(handles.output,'saved_cellpositions',saved_cellpositions);

shank_num = find(eventdata.NewValue==shank_buttons);
setappdata(handles.output,'shank_num',shank_num);
cell_positions=saved_cellpositions{shank_num};

% determine how many cells we have
%num_cells=max(num_cells(:,shank_num));
num_cells=num_cells(:,shank_num);
num_files=getappdata(handles.output,'num_files');

% create plot axes (shank1 for now)
num_rows=max(num_cells); num_cols=num_files;
margin_top=0.07; margin_left=0.05;
margin_x=.002; margin_y=.02;
plot_width=(1-(num_cols-1)*margin_x-margin_left)/num_cols; plot_height=(1-(num_rows-1)*margin_y-margin_top)/num_rows;

% determine axis centers and borders
ax_centers=zeros(num_rows,num_cols,2);
ax_centers(:,:,1)=repmat(linspace(plot_width/2+margin_left,1-plot_width/2,num_cols),num_rows,1);
ax_centers(:,:,2)=repmat((linspace(1-margin_top-plot_height/2,plot_height/2,num_rows))',1,num_cols);

ax_pos=zeros(num_rows,num_cols,4);
ax_pos(:,:,1)=ax_centers(:,:,1)-plot_width/2;
ax_pos(:,:,2)=ax_centers(:,:,2)-plot_height/2;
ax_pos(:,:,3)=plot_width;
ax_pos(:,:,4)=plot_height;

% we need to use the old cell positions if it's been modified already
if(isempty(cell_positions)), 
    cell_positions=repmat((1:num_cells)',1,num_files);
end

% generate axes
if(~isempty(a))
    delete(a);
end
a=zeros(num_rows,num_cols);

for i = 1:num_rows % which cell
    for j = 1:num_cols % which file
        a(i,j)=axes('position',squeeze(ax_pos(i,j,:)),'parent',handles.main_panel,'ytick',[],'xtick',[]); %#ok<LAXES>
        
        % plot templates
        if(i<=num_cells(j))
            wfs=templates{j,shank_num};
            plot(a(i,j),squeeze(wfs(i,:,:))','linesmooth','on','linewidth',2,'hittest','off');
            set(a(i,j),'buttonDownFcn',{@select_axes,handles,a(i,j)});
        end
        set(a(i,j),'color','k','ycolor',[.5 .5 .5],'xcolor',[.5 .5 .5]);
        set(a(i,j),'xticklabel',[],'yticklabel',[]);
        grid(a(i,j),'on');
        set(a(i,j),'ylim',[0 1]);
    end
end

% re-order cells
for i = 1:num_cols
    
    % scroll through rows
    for j = 1:size(a,1)
        set(a(j,i),'position',ax_pos(cell_positions(j,i),i,:));
    end
    
    a(:,i)=a(cell_positions(:,i),i);
end

set(handles.output,'WindowButtonUpFcn',{@deselect_axes,handles});


% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.status_light,'backgroundcolor','red');
set(handles.status_label,'String','Grabbing spikes...');

saved_cellpositions=getappdata(handles.output,'saved_cellpositions');
mat_files=getappdata(handles.output,'mat_files');
num_cells=getappdata(handles.output,'num_cells');
num_files=getappdata(handles.output,'num_files');
max_cells=max(num_cells);

% cycle through, get spikes and concatenate them
num_shanks=length(getappdata(handles.output,'shank_buttons'));
spikes=cell(num_shanks,num_files,max(max_cells));

for i = 1:num_shanks
    for j = 1:num_files
        % load idx and timestamps
        load(mat_files{j},sprintf('idx_%g',i),sprintf('timestamps_%g',i));
        for k = 1:max_cells(i)
            eval(sprintf('spikes{i,j,k}=timestamps_%g(idx_%g==k);',i,i));
        end
    end
end

% now we have to reorder
set(handles.status_label,'String','Merging datasets...');
for i = 1:num_shanks
    cell_positions=saved_cellpositions{i};
    for j = 1:num_files
        pos=cell_positions(:,j);
        spikes(i,j,1:length(pos))=spikes(i,j,pos);
    end
end

% concatenate
spikes=permute(spikes,[2 1 3]);
spikes=reshape(spikes,num_files,[])';

num_cells=size(spikes,1);
for i = 1:num_cells
    spikes{i,1}=cat(2,spikes{i,:});
end
spikes=spikes(:,1);
locs=~cellfun(@isempty,spikes);
spikes=spikes(locs);
shank_nums=repmat(1:num_shanks,1,max(max_cells));
shank_nums=shank_nums(locs)'; %#ok<NASGU>
spikes=cellfun(@sort,spikes,'UniformOutput',false); %#ok<NASGU>

% save output file
basename=getappdata(handles.output,'basename');
pathname=getappdata(handles.output,'pathname');
filename=sprintf('%s%s.mat',pathname,basename);

if(exist(filename,'file'))
    save(filename,'-append','spikes','shank_nums');
else
    save(filename,'spikes','shank_nums');
end

set(handles.status_label,'String',sprintf('Saved to %s.',filename));
set(handles.status_light,'backgroundcolor','green');