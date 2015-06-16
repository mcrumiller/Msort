function show_all_waveforms(handles)
%SHOW_ALL_WAVEFORMS   Generate waveform plot
%   SHOW_ALL_WAVEFORMS Generates waveform plots for the selected neurons, shown on the selected channels. Additional
%   errorbars are also generated.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms palette;
MAX_WFS=str2double(get(handles.max_wfs_field,'String')); % max # waveforms to display

% create empty wf_plots matrix
channels=getappdata(handles.output,'channels');
num_channels=length(channels);

% clear any plots already around
wf_plots=getappdata(handles.output,'wf_plots');
if(~isempty(wf_plots)), delete(wf_plots(ishandle(wf_plots) & wf_plots~=0)); end
errorbar_plots=getappdata(handles.output,'errorbar_plots');
if(~isempty(errorbar_plots)), delete(errorbar_plots(ishandle(errorbar_plots) & errorbar_plots~=0)); end
gridline_plots=getappdata(handles.output,'gridline_plots');
if(~isempty(gridline_plots)), delete(gridline_plots(ishandle(gridline_plots) & gridline_plots~=0)); end

% clear current plots if they exist
ax_waveforms=getappdata(handles.output,'ax_waveforms');
delete(cell2mat(get(ax_waveforms,'children')));

% determine which channels to plot, and which groups to plot
cluster_numbers=getappdata(handles.output,'cluster_number');
if(isempty(cluster_numbers)), cluster_numbers=0; end
cluster_colors=getappdata(handles.output,'cluster_colors');
channels=getappdata(handles.output,'channels');
num_clusters=getappdata(handles.output,'num_clusters');

%waveforms=getappdata(handles.output,'waveforms');
idx=getappdata(handles.output,'idx');
t=getappdata(handles.output,'t')*1e6;

% determine ticks
supermaxval=max(waveforms(:));
superminval=min(waveforms(:));
if(get(handles.fit_wf_checkbox,'Value'))
    locs=false(1,size(waveforms,2));
    for i = 1:length(cluster_numbers)
        locs=locs | idx==cluster_numbers(i);
    end
    wfs=waveforms(:,locs,:);
    maxval=max(wfs(:));
    minval=min(wfs(:));
else
    maxval=supermaxval;
    minval=superminval;
end

if(isempty(minval)), minval=0; maxval=1; end

xtick=[-200 -100 0 100 200 300 400]; % not sure if this should be hardcoded
ytick=sort(floor(superminval*1.05)*2:2:ceil(supermaxval*1.05)*2);

set(ax_waveforms,'visible','off','xlim',[t(1) t(end)],...
    'xcolor',palette.text,'ycolor',palette.text,'xtick',xtick,'ytick',ytick,'box','on','linewidth',1.5);

display_tools=getappdata(handles.output,'display_cleanup_tools');
t_err=double(reshape([t; t; nan(1,length(t))],1,[]));

include_zero=true;
mean_wfs=get_mean_waveforms(waveforms,idx,include_zero);

% generate plots
% generate new plot handles
wf_plots=zeros(num_channels,num_clusters);
errorbar_plots=zeros(num_channels,num_clusters);
% each gridline is its own plot. Might update this in future to make single
% line since it's faster and easier to handle.
gridline_plots=zeros(1,num_channels);
ch_title_plots=zeros(1,num_channels);
max_wfs=zeros(length(channels),num_clusters);
min_wfs=zeros(length(channels),num_clusters);
for i = 1:length(channels)
    for j = 0:num_clusters-1
        wfs=squeeze(waveforms(i,idx==j,:))';
        if(isvector(wfs)),wfs=wfs'; end
        num_wfs=size(wfs,2);
        if(num_wfs>0)
            max_wfs(i,j+1)=max(wfs(:));
            min_wfs(i,j+1)=min(wfs(:));
        end
        
        % Determine std deviation error bars (display cleanup tools)
        std_dev_pts = std(wfs,[],2);
        mean_wf_pts = squeeze(mean_wfs(j+1,i,:));
        errorbar_value=get(handles.errorbar_slider,'Value');
        std_dev_pts = [mean_wf_pts-errorbar_value*std_dev_pts mean_wf_pts+errorbar_value*std_dev_pts]';
        
        
        if(num_wfs)>0
            % only plot up to MAX_WFS waveforms
            if(num_wfs>MAX_WFS)
                ind=round(linspace(1,num_wfs,MAX_WFS));
                wfs=wfs(:,ind,:);
                num_wfs=MAX_WFS;
            end
            
            wfs=[wfs; nan(1,num_wfs)]; %#ok<AGROW>
            t_wfs=repmat(t',1,num_wfs);
            t_wfs=[t_wfs; nan(1,num_wfs)]; %#ok<AGROW>
            if(any(j==cluster_numbers))
                wf_plots(i,j+1)=plot(ax_waveforms(i),...
                    t_wfs(:),wfs(:),'color',cluster_colors(j+1,:),'linesmooth','on','visible','on','hittest','off');
            else
                wf_plots(i,j+1)=plot(ax_waveforms(i),...
                    t_wfs(:),wfs(:),'color',cluster_colors(j+1,:),'linesmooth','on','visible','off','hittest','off');
            end
        end
        gridline_plots(i)=minorgrid(xtick,ytick,ax_waveforms(i),false,[t(1) t(end) superminval*1.2 supermaxval*1.2]);
        
        % plot std devs
        std_dev_pts=reshape([std_dev_pts; nan(1,length(t))],1,[]);
        errorbar_plots(i,j+1)=plot(ax_waveforms(i),t_err,std_dev_pts,'color',palette.red,'visible','off','hittest','off');
        if(display_tools && j==cluster_numbers(1))
            set(errorbar_plots(i),'visible','on');
        end
    end
    
    % put channel tag
    x=double(t(1)+(t(end)-t(1))*.01);
    y=double(maxval*1.05*.99);
    
    % add channel tag
    ch_title_plots(i)=text(x,y,sprintf('Ch %g',channels(i)),'fontsize',20,'backgroundcolor',palette.ax_background,...
        'parent',ax_waveforms(i),'color',palette.text,'horizontalalignment','left',...
        'verticalalignment','top','hittest','off');
end

setappdata(handles.output,'min_wfs',min_wfs);
setappdata(handles.output,'max_wfs',max_wfs);
setappdata(handles.output,'wf_plots',wf_plots);
setappdata(handles.output,'errorbar_plots',errorbar_plots);
setappdata(handles.output,'gridline_plots',gridline_plots);
setappdata(handles.output,'ch_title_plots',ch_title_plots);
%set([gridline_plots(:)' wf_plots(:)' errorbar_plots(:)' ch_title_plots(:)'],'hittest','off');
children=get(ax_waveforms,'children');
%if(~iscell(children)), children={children}; end
if(iscell(children)), children=cell2mat(children); end
set(children,'hittest','off');