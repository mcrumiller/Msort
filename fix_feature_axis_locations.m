function fix_feature_axis_locations(handles)
%FIX_FEATURE_AXIS_LOCATIONS   Re-distribute waveform axis locations (occurs when channels are selected and deselected)
%   [] = FIX_FEATURE_AXIS_LOCATIONS() 
%   
%   Input:  
%   Output: 
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
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