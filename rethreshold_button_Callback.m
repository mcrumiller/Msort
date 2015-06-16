function rethreshold_button_Callback(~, ~, handles)
%RETHRESHOLD_BUTTON_CALLBACK   Set spike threshold levels
%   RETHRESHOLD_BUTTON_CALLBACK Sets the spike thresholding by either applying an automated threshold level by user
%   prompt, or the user manually selects a threshold with a mouse.
%
%   RUN_THRESHOLD_SPIKES() is called from this function, which is the function that actually applies the thresholds.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
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
run_threshold_spikes(handles); % actually apply thresholds
plot_voltages(handles);

% play beep
if(get(handles.sound_checkbox,'Value')), stop_alert; end
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);