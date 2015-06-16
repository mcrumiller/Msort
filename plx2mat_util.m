function varargout = plx2mat_util(varargin)
% PLX2MAT_UTIL MATLAB code for plx2mat_util.fig
%      PLX2MAT_UTIL, by itself, creates a new PLX2MAT_UTIL or raises the existing
%      singleton*.
%
%      H = PLX2MAT_UTIL returns the handle to a new PLX2MAT_UTIL or the handle to
%      the existing singleton*.
%
%      PLX2MAT_UTIL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLX2MAT_UTIL.M with the given input arguments.
%
%      PLX2MAT_UTIL('Property','Value',...) creates a new PLX2MAT_UTIL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plx2mat_util_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plx2mat_util_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plx2mat_util

% Last Modified by GUIDE v2.5 24-Jan-2013 18:35:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @plx2mat_util_OpeningFcn, ...
    'gui_OutputFcn',  @plx2mat_util_OutputFcn, ...
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


% --- Executes just before plx2mat_util is made visible.
function plx2mat_util_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plx2mat_util (see VARARGIN)

% Choose default command line output for plx2mat_util
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes plx2mat_util wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plx2mat_util_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in file_listbox.
function file_listbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to file_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns file_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from file_listbox
loc=get(hObject,'Value');
num_selected=length(loc);
full_filenames=getappdata(handles.output,'full_filenames');

% determine if we have any
code=getappdata(handles.output,'code');
num_unprocessed=sum(code(loc)==0 | code(loc)==2);
if(num_unprocessed>0)
    set(handles.convert_button,'enable','on');
else
    set(handles.convert_button,'enable','off');
end

if(num_selected>1)
    set(handles.directory_label,'String',sprintf('%g files selected (%g unprocessed) ',length(loc),num_unprocessed));
else
    set(handles.directory_label,'String',sprintf('%s ',full_filenames{loc}));
end

% --- Executes during object creation, after setting all properties.
function file_listbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to file_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_directory_button.
function load_directory_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to load_directory_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name=uigetdir(pwd);
if(~folder_name)
    set(handles.status_label,'String','Please select a directory');
    set(handles.status_light,'BackgroundColor','green');
    return
end
folder_name=strrep(folder_name,'\','/');

% grab all files from the folder and subfolders, populate list
[filenames,pathnames]=get_files(folder_name,'.plx');
full_filenames=strcat(pathnames(:),'/',filenames(:));
setappdata(handles.output,'full_filenames',full_filenames);

% keep track of what we've done for each file
code=zeros(length(full_filenames),1);

% run a scan of each file; determine whether a plx file exists for each
% set appropriate code:
%      0 = unprocessed (xml & mat file found)
%      1 = processed (xml & mat file found)
num_files=length(full_filenames);
deepread=get(handles.deepread_tag,'Value');

if(deepread)
    h=waitbar(0,sprintf('Checking (%g/%g)...',0,num_files),'name','Checking .mat files for errors...');
end
for i = 1:num_files
    abort=false;
    current_filename=full_filenames{i};
    temp_filename=current_filename(1:end-4);
    xml_version=sprintf('%s.xml',temp_filename);
    mat_version=sprintf('%s.mat',temp_filename);
    
    if(deepread)
        waitbar(i/num_files,h,sprintf('Checking (%g/%g)...',i,num_files));
    end
    
    % check if both exist
    if(exist(xml_version,'file') && exist(mat_version,'file'))
        % do a deep read on .mat file
        if(deepread)
            
            % get all channels from XML file
            s=xml2struct(xml_version);
            num_shanks=length(s.experiment.shank);
            channels=[];
            for j = 1:num_shanks
                if(abort), break; end
                num_channels=length(s.experiment.shank{j}.channel);
                for c = 1:num_channels
                    ch=str2double(s.experiment.shank{j}.channel{c}.Attributes.num);
                    if(~isfield(s.experiment.shank{j}.channel{c},'enabled'))
                        abort=true; break;
                    end
                    enabled=s.experiment.shank{j}.channel{c}.enabled.Text;
                    if(strcmp(enabled,'true'))
                        channels=[channels ch]; %#ok<AGROW>
                    end
                end
            end
            if(abort)
                % delete xml and .mat file
                delete(xml_version,mat_version);
                code(i)=0;
                continue;
            end
            
            % make sure we can write to the .mat file
            writeable=true;
            try
                test_write=[]; %#ok<NASGU>
                save(mat_version,'-v7.3','-append','test_write');
            catch %#ok<CTCH>
                writeable=false;
            end
            
            if(~writeable)
                delete(xml_version,mat_version);
                code(i)=0;
                continue;
            end
            
            % make sure we have all data_# variables
            vars=whos('-file',mat_version);
            vars={vars.name}; expr=regexp(vars,'data_(?<channel>\d*$)','names');
            locs=cellfun(@(x) ~isempty(x),expr); vars_data=expr(locs);
            variable_list=cellfun(@(x) str2double(x.channel),vars_data);
            
            
            
            if(isequal(sort(channels),sort(variable_list)))
                code(i)=1;
            else
                code(i)=2;
            end
        else
            code(i)=1;
        end
    else
        code(i)=0;
    end
end
if(deepread), delete(h); end

setappdata(handles.output,'code',code);

update_listbox(handles);

% --- Returns all files from within a folder, and calls itself recursively
% on sub-folders. The complete filename is returned.
function [filenames,pathnames] = get_files(directory, extension)
results=dir(directory);
names={results(:).name};

% filter out this and parent directory
parents=strcmp(names,'.') | strcmp(names,'..');
results(parents)=[]; names(parents)=[];

if(isempty(results)), filenames=[]; pathnames=[]; return; end

dirs=logical(cell2mat({results(:).isdir}));

% filter out filenames that don't end with the proper extension
if(extension(1)=='.'),expr=sprintf('%s$',extension);
else expr=sprintf('\\%s$',extension);
end
badlocs=cellfun(@isempty,regexpi(names,expr));
badlocs=badlocs & (~dirs);
names(badlocs)=[]; dirs(badlocs)=[];

% list actual files
filenames=names(~dirs);
num_files=length(filenames);
pathnames=repmat({directory},num_files,1);

% scroll through remaining directories
dirlocs=find(dirs);
for i = 1:length(dirlocs)
    [filenames2,pathnames2]=get_files(sprintf('%s/%s',directory,names{dirlocs(i)}),extension);
    filenames = [filenames(:); filenames2(:)];
    pathnames = [pathnames(:); pathnames2(:)];
end

% --- Executes on button press in convert_button.
function convert_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to convert_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global electrode_type;
loc=get(handles.file_listbox,'Value');
code=getappdata(handles.output,'code');
full_filenames=getappdata(handles.output,'full_filenames');
listnames=get(handles.file_listbox,'String');

% determine how many we have to process
currentcode=code(loc);
toprocess=loc(currentcode==0 | currentcode==2);
num_toprocess=length(toprocess);

if(num_toprocess==0), return; end

% load electrode types from file
[electrode_names,electrode_configuration]=load_electrode_types('electrodes.cfg');

% create electrode selection dialog
scrsz = get(0,'ScreenSize'); width = scrsz(3)*.2; height = scrsz(4)*.1;
left=(scrsz(3)-width)*.5; bot=(scrsz(4)-height)*.5;
dial=dialog('position',[left bot width height]);
label_width=.8;
margin=(1-label_width)/2;
uicontrol('Style','text','String','Please select electrode type',...
    'fontsize',16,'fontweight','bold','horizontalalignment','left',...
    'Units','normalized','position',[margin .5 label_width .4],...
    'Parent',dial);
menu_width=.8;
margin=(1-menu_width)/2;
menu=uicontrol('Style', 'popup',...
    'String', electrode_names, 'fontsize',12,...
    'Units','normalized','Position',[margin 0 menu_width .6],'Parent',dial);
button_width=.2;
margin=(1-button_width)/2;
uicontrol('Style','pushbutton','String','OK','fontsize',16,'Units','normalized',...
    'Position',[margin .05 button_width .25],'Parent',dial,...
    'Callback',{@set_electrode_type,menu});
waitfor(dial);

% set ones to process in bold
listnames(toprocess) = strcat('<html><font color=''red''><b>',listnames(toprocess),'</b></font></html>');
set(handles.file_listbox,'String',listnames);
set(handles.status_light,'BackgroundColor','red');

failed=[];
for i = 1:num_toprocess
    % bold current file
    current_loc=toprocess(i);
    current_file=full_filenames{current_loc};
    set(handles.file_listbox,'String',listnames);
    set(handles.file_listbox,'Value',current_loc);
    
    % plx2 mat the current file
    set(handles.status_label,'String',sprintf('Converting %s...',current_file));
    try
        plx2mat(current_file,electrode_type);
    catch %#ok<CTCH>
        listnames{current_loc}=strcat('<html><font color=''red''>',listnames{current_loc},'</font></html>');
        set(handles.file_listbox,'String',listnames);
        failed=[failed current_loc]; %#ok<AGROW>
        continue;
    end
    set(handles.directory_label,'String',sprintf('%g/%g processed ',i,num_toprocess));
    code(toprocess)=1;
    setappdata(handles.output,'code',code);
    
    listnames{current_loc}=full_filenames{current_loc};
    set(handles.file_listbox,'String',listnames);
end
set(handles.file_listbox,'Value',[]);

% success
if(isempty(failed))
    set(handles.status_label,'String','Done!');
    set(handles.status_light,'BackgroundColor','green');
    % set to orange
else
    set(handles.status_label,'String',sprintf('%g files failed conversion (highlighted in bold).',length(failed)));
    set(handles.status_light,'BackgroundColor',[255 200 0]/255);
end
set(handles.directory_label,'String','No files selected ');


% --- Update listbox view
function update_listbox(handles)
code=getappdata(handles.output,'code');
full_filenames=getappdata(handles.output,'full_filenames');
num_files=length(full_filenames);
prefix=cell(num_files,1); suffix=cell(num_files,1);

prefix(code==0)={'<html><font color=''red''>* '};
prefix(code==1)={''};
prefix(code==2)={'<html><font color=''#FFA000''>* '}; % orange


suffix(code==true)={''};
suffix(code==0 | code==2)={'</font></html>'};

display_full_filenames=strcat(prefix,full_filenames,suffix);

% set listbox strings
set(handles.file_listbox,'String',display_full_filenames);
set(handles.status_label,'String',sprintf('%g files found.',num_files));


function set_electrode_type(~,~,menu_obj)
global electrode_type;
s=get(menu_obj,'String');
electrode_type=s{get(menu_obj,'Value')};
%setappdata(0,'electrode_type',result);
delete(get(menu_obj,'parent'));


% --- Executes on button press in deepread_tag.
function deepread_tag_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to deepread_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deepread_tag
