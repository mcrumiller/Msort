function plx2mat(plx_file,electrode_type,selected_shanks)
%--------------------------------------------------------------------------
% plx2mat2.m - Gets as input a plx file containing continuous channel
% traces. An associated .xml file is generated if it doesn't exist, along
% with a .mat file that contains all of the filtered and unfiltered spike
% data.
%
% If provided, the user can specify which specific shanks should be
% converted to .mat; this saves time if one only wants to convert a subset
% of the channels to .mat.
%
% Usage: plx2mat2(plx_file,selected_shanks);
%
% Input:  plx_file                  * .plx file containing continuous
%                                      recording
%         selected_shanks           * vector of shank numbers
% Output:
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
% get subdirectory and base name of file
expr = regexp(plx_file,'(?<dirname>(.*[\\/])*)(?<basename>.*)(?<ext>\.plx|\.xml)+','names');
if(isempty(expr)), error('Invalid filename.'); end
dirname = expr.dirname;
if(isempty(dirname))
    dirname = pwd;
    if(dirname(end)~='/' || dirname(end)~='\')
        dirname=sprintf('%s\\',dirname);
    end
else dirname=expr.dirname;
end;
basename = expr.basename;

% get file names
plx_file = sprintf('%s%s.plx',dirname,basename);
xml_file = sprintf('%s%s.xml',dirname,basename);
mat_file = sprintf('%s%s.mat',dirname,basename);

if(~exist(mat_file,'file'))
    % create empty mat file
    save(mat_file,'-v7.3','');
end

% Load data from parameter file
if(~exist(plx_file,'file')), error('.plx file does not exist.'); end
if(~exist(xml_file,'file'))
    fprintf('.xml file does not exist.  Generating standard xml...\n');
    if(exist('electrode_type','var'))
        generate_standard_Plexon_xml(plx_file,electrode_type);
    else
        generate_standard_Plexon_xml(plx_file);
    end
end

% grab shank and channel information from the .plx file
s = xml2struct(xml_file);

num_shanks = length(s.experiment.shank);

% Determine which channels have already been processed in the .mat file
file_vars=whos('-file',mat_file);
vars={file_vars.name}; expr=regexp(vars,'data_(?<channel>\d*$)','names');

% Determine if we have 'Ts' or not and, if so, get its length
Ts_loc=find(strcmp(vars,'Ts'),1,'first');
contains_Ts=any(Ts_loc);
events_loc=find(strcmp(vars,'events'),1,'first');
contains_events=any(events_loc);

Ts_length=0;
if(contains_Ts), Ts_length=max(file_vars(Ts_loc).size); end
locs=cellfun(@(x) ~isempty(x),expr); vars=expr(locs);
completed_flat=sort(cellfun(@(x) str2double(x.channel),vars));

shanks=1:num_shanks;
process_shanks=true(1,num_shanks);
if(exist('selected_shanks','var') && ~isempty(selected_shanks))
    process_shanks(selected_shanks)=false;
    process_shanks=~process_shanks;
end

% load channel and shank information
channels=cell(1,num_shanks);
disabled=cell(1,num_shanks);
completed=cell(1,num_shanks);
for i = 1:num_shanks
    if(~process_shanks(i)),continue; end
    len=length(s.experiment.shank{i}.channel);
    shank=str2double(s.experiment.shank{i}.Attributes.num); shanks(i)=shank;
    channels{i}=zeros(1,len);
    disabled{i}=false(1,len);
    
    % Note: this special case is due to the fact that xml2struct doesn't
    % place single items into cell arrays
    if(len==1)
        channels{i} = str2double(s.experiment.shank{i}.channel.Attributes.num);
        if(strcmp(s.experiment.shank{i}.channel.enabled.Text,'false'))
            disabled{i}=true;
        end
        if(any(completed_flat==channels{i})), completed{i}=true;
        else completed{i}=false;
        end
    else
        for j = 1:len
            channels{i}(j)=str2double(s.experiment.shank{i}.channel{j}.Attributes.num);
            if(strcmp(s.experiment.shank{i}.channel{j}.enabled.Text,'false'))
                disabled{i}(j)=true;
            else
                disabled{i}(j)=false;
            end
            if(any(completed_flat==channels{i}(j))), completed{i}(j)=true;
            else completed{i}(j)=false;
            end
        end
    end
end

channels_flat=horzcat(channels{:});
disabled_flat=horzcat(disabled{:});
completed_flat=horzcat(completed{:});

num_remaining=length(channels_flat(~disabled_flat & ~completed_flat));

Fs=str2double(s.experiment.ADRate.Text);

% run through each shank, process .mat files
num_total=num_remaining; num_completed=1;
if(num_total~=0)
    
    h=waitbar(0,sprintf('Processing %g/%g channels...',0,num_remaining));
    for i = 1:num_shanks
        if(~process_shanks(i)), continue; end
        for j = 1:length(channels{i})
            current_channel=channels{i}(j);
            
            % skip if disabled
            channel_loc=find(current_channel==channels_flat);
            if(disabled_flat(channel_loc) || completed_flat(channel_loc)), continue; end
            
            waitbar(num_completed/num_total,h,sprintf('Processing Channel %g/%g',num_completed,num_total)); drawnow;
            
            % grab channel data from plx file
            [data,Ts]=get_data_from_plx_c(plx_file,current_channel);
            data=int16(data);
            
            % save timestamps if this one is larger
            if(length(Ts)>Ts_length),
                save(mat_file,'-v7.3','-append','Ts');
                Ts_length=length(Ts);
            end
            
            % we have a bad channel, disable it and update struct/xml file
            if(~any(data))
                disabled{i}(j)=true;
                s.experiment.shank{i}.channel{j}.enabled.Text='false';
                struct2xml(s,xml_file);
                % data is fine, save to .mat file
            else
                % save unfiltered version
                unfiltered_var_name=sprintf('data_%g_unfiltered',current_channel);
                eval(sprintf('%s = data;',unfiltered_var_name));
                save(mat_file,'-append','-v7.3',unfiltered_var_name);
                clear(unfiltered_var_name);
                
                % save filtered version
                waitbar(num_completed/num_total,h,sprintf('Processing Channel %g/%g (normalizing data)',num_completed,num_total)); drawnow;
                data=single(data);
                data=filter_voltage(data,Fs);
                
                % set standard deviation to 1
                std_data=std(data);
                data=data./std_data; %#ok<NASGU>
                
                % save data variable
                waitbar(num_completed/num_total,h,sprintf('Processing Channel %g/%g (saving data)',num_completed,num_total)); drawnow;
                var_name=sprintf('data_%g',current_channel);
                std_varname=sprintf('std_%g',current_channel);
                eval(sprintf('%s=std_data;',std_varname));
                eval(sprintf('%s=data;',var_name)); clear('data');
                save(mat_file,'-v7.3','-append',var_name,std_varname);
                clear('data',var_name);
            end
            
            num_completed=num_completed+1;
        end
    end
    
end


if(~contains_events)
    % grab event information
    [~, ~, evcounts] = plx_info(plx_file,false);
    if(~any(evcounts))
        [~, ~, evcounts] = plx_info(plx_file,true);
    end
    if(~any(evcounts))
        events=[]; %#ok<NASGU>
    else
        [~,evloc]=max(evcounts);
        [~, evchans] = plx_event_chanmap(plx_file);
        event_channel=evchans(evloc);
        [~, events] = plx_event_ts(plx_file, event_channel);
        events=events'; %#ok<NASGU>
    end
    save(mat_file,'-append','events');
end

delete(h);