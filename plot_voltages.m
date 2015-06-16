function plot_voltages(handles)
%PLOT_VOLTAGES   Generate voltage plot
%   PLOT_VOLTAGES generates a voltage plot of each channel, accompanied by a horizontal slider. It also allows for
%   temporal zooming in/out.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
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