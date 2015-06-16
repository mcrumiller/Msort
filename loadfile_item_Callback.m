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