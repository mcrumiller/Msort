function plot_mean_waveforms(handles)
%PLOT_MEAN_WAVEFORMS   Generate plot of average waveform for each neuron
%   PLOT_MEAN_WAVEFORMS() Plots the mean waveform of each channel for each neuron and displays the result in the
%   small sub-axes located to the left of the main plot.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms palette;
idx=getappdata(handles.output,'idx');
num_clusters=getappdata(handles.output,'num_clusters');
if(isempty(num_clusters) || num_clusters==0), return; end

t=getappdata(handles.output,'t');

set(handles.status_light,'BackgroundColor',palette.red);
set(handles.status_label,'String','Plotting waveforms...');

% Plot mean waveforms
include_zero=true; % tells get_mean_waveforms.m that idx starts at 0, not 1
mean_wfs = get_mean_waveforms(waveforms,idx,include_zero);
max_val = max(max(max(mean_wfs)));
min_val = min(min(min(mean_wfs)));
span=max_val-min_val; margin=span*0.05;
mtop=max_val+margin; mbot=min_val-margin;

dx=(t(end)-t(1))/4;
xgrid=t(1)+[dx 2*dx 3*dx];
dy=(mtop-mbot)/4;
ygrid=mbot+[dy 2*dy 3*dy];

channel_colors=getappdata(handles.output,'channel_colors');
cluster_colors = getappdata(handles.output,'cluster_colors');

% Determine location for axes. This is as wide as and up to the bottom of
% the waveform_control_panel and

panel_selected=getappdata(handles.output,'panel_selected');
features_pos=get(handles.features_panel,'position');
feature_width=features_pos(3);
pos=get(handles.feature_control_panel,'position');
total_width=(pos(3)/feature_width)*.99;
if(strcmp(panel_selected,'features'))
    total_height=pos(2)+.015;
elseif(strcmp(panel_selected,'cell_stats'))
    total_height=.985;
end

wf_axes=zeros(1,num_clusters);
frame_handles=zeros(1,num_clusters);
wf_buttons = zeros(1,num_clusters);

plot_width=total_width;
plot_height=total_height/num_clusters;

% plot waveforms
for i = 1:num_clusters
    wf_axes(i) = axes('Parent',handles.features_panel,'units','normalized'); %#ok<LAXES>
    set(wf_axes(i),'position',[0 total_height-i*plot_height plot_width plot_height]);
    hold all;
    if(i==1)
        set(wf_axes(i),'xtick',[],'ytick',[],'xcolor',palette.text,'ycolor',palette.text,'color',palette.background);
    else
        set(wf_axes(i),'xtick',[],'ytick',[],'xcolor',palette.text,'ycolor',palette.text,'color',palette.ax_background);
    end
    
    set(wf_axes(i),'colororder',channel_colors); hold(wf_axes(i),'all');
    plot(t,squeeze(mean_wfs(i,:,:))','linewidth',2,'linesmooth','on','Parent',wf_axes(i));
    axis(wf_axes(i),[t(1) t(end) mbot mtop]);
    
    minorgrid(xgrid,ygrid,wf_axes(i));
    
    set(wf_axes(i),'box','on','xticklabel',[],'yticklabel',[],'ButtonDownFcn',{@select_axes,handles});
    
    % determine square patch location
    a = axis(wf_axes(i));
    width=a(2)-a(1); margin_x = width*.01;
    height=a(4)-a(3); margin_y = height*.01;
    square_width = width * .2;
    square_height = height * .2;
    
    % add patch
    x=[a(2)-margin_x-square_width a(2)-margin_x-square_width a(2)-margin_x a(2)-margin_x];
    y=[a(3)+margin_y a(3)+margin_y+square_height a(3)+margin_y+square_height a(3)+margin_y];
    patch(x,y,cluster_colors(i,:));
    
    % place text immediately before patch
    text(a(2)-2*margin_x-square_width,a(3)+margin_y+square_height/2,sprintf('%g',length(find(idx==i-1))),...
        'fontsize',12,'fontweight','bold','fontname','Cambria','verticalalignment','middle',...
        'horizontalalignment','right','backgroundcolor',palette.ax_background,'margin',1,'color',cluster_colors(i,:),'Parent',wf_axes(i));
    
    % set frame
    a = axis(wf_axes(i));
    span_x=a(2)-a(1); span_y=a(4)-a(3);
    margin_x=span_x*.001; margin_y=span_y*.002;
    x=[a(1) a(1) a(2)-margin_x a(2)-margin_x a(1)];
    y=[a(3)+margin_y a(4) a(4) a(3)+margin_y a(3)+margin_y];
    frame_handles(i)=plot(x,y,'linewidth',10,'color',cluster_colors(i,:),'visible','off','Parent',wf_axes(i));
    
    % turn off hittest for all children
    c=get(wf_axes(i),'Children');
    set(c,'hittest','off');
end
selected_axes=getappdata(handles.output,'selected_axes');
set(frame_handles(selected_axes),'visible','on');
setappdata(handles.output,'wf_axes',wf_axes);
setappdata(handles.output,'frame_handles',frame_handles);
setappdata(handles.output,'wf_buttons',wf_buttons);

% Reset labels
set(handles.status_label,'String','Done!');
set(handles.status_light,'BackgroundColor',palette.green);