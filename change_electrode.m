function change_electrode(~,~,electrode_menu,panel,configurations,geometry,handles)
%CHANGE_ELECTRODE   switch electrode configuration
%   CHANGE_ELECTRODE switches the current electrode configuration from one type of electrode to another by remapping
%   channels. This involves a somewhat extensive level of reassignment of data and, as such, may take some time as the
%   voltage traces are loaded and remapped into separate portions of the .mat file.
%
%   This function may not currently work properly.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
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