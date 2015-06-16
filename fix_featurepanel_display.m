function fix_featurepanel_display()
%FIX_FEATUREPANEL_DISPLAY   Set features panel to appropriate settings
%   FIX_FEATUREPANEL_DISPLAY correctly displays what is supposed to be displayed on the feature panel.
%   Note: a lot of the features set here may end up being redundantly set.
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
margin_left=.225; margin_right=0;
margin_bot=.04; margin_top=0;

% reduces the number of "set" commands
visibles=[]; invisibles=[];

% create ax_waveforms if missing
ax_waveforms=getappdata(handles.output,'ax_waveforms');
channels=getappdata(handles.output,'channels');
num_channels=length(channels);
if(isempty(ax_waveforms) || ~any(ishandle(ax_waveforms)))
    ax_waveforms=zeros(1,num_channels);
    for i = 1:num_channels
        ax_waveforms(i)=axes('parent',handles.features_panel,'color',palette.ax_background,'visible','off',...
            'xcolor',palette.text,'ycolor',palette.text,'xtick',[],'ytick',[],'linewidth',1.5); %#ok<LAXES>
        hold(ax_waveforms(i),'all');
    end
end
setappdata(handles.output,'ax_waveforms',ax_waveforms);

% create ax_features if missing
ax_features=getappdata(handles.output,'ax_features');
if(isempty(ax_features) || ~ishandle(ax_features))
    ax_features=axes('position',[margin_left margin_bot 1-margin_left-margin_right 1-margin_top-margin_bot],...
        'parent',handles.features_panel,'color',palette.ax_background,'visible','off'); hold(ax_features,'all');
    grid(ax_features,'on');
    setappdata(handles.output,'ax_features',ax_features);
end

% determine axis selection type
view_mode=logical([get(handles.waveform_view_button,'value') get(handles.split_view_button,'value') get(handles.features_view_button,'value')]);

% determine which channels are selected
selected_channels=getappdata(handles.output,'selected_channels');
if(~any(selected_channels)), selected_channels=~selected_channels; end

% determine which groups are selected
selected_groups=getappdata(handles.output,'selected_axes');

% hide everything in channels that are not selected
hidden_axes=ax_waveforms(~selected_channels);
children=get(hidden_axes,'children');
if(iscell(children)), children=cell2mat(children); end
invisibles=[invisibles hidden_axes(:)' children(:)'];

% fix wf axis limits
min_wfs=getappdata(handles.output,'min_wfs');
max_wfs=getappdata(handles.output,'max_wfs');
checked=get(handles.fit_wf_checkbox,'Value');

% fit waveforms to selected channels
if(checked && any(selected_groups))
    min_wf=min_wfs(selected_channels,selected_groups);
    min_wf=min(min_wf(:));
    max_wf=max_wfs(selected_channels,selected_groups);
    max_wf=max(max_wf(:));
else
    min_wf=min(min_wfs(:));
    max_wf=max(max_wfs(:));
end
if(min_wf==0 || isempty(min_wf)),min_wf=0; max_wf=1; end
set(ax_waveforms,'ylim',[min_wf max_wf]);

% adjust channel labels as well
t=getappdata(handles.output,'t')*1e6;
if(isempty(t)), return; end
ch_title_plots=getappdata(handles.output,'ch_title_plots');
spanx=t(end)-t(1); spany=max_wf-min_wf;
set(ch_title_plots,'position',[t(1)+spanx*.01 min_wf+spany*.98]);

% Stuff to make visible:
%   -waveform plots (only selected groups, on selected channel axes)
%   -gridlines (only selected groups, on selected channel axes)
if(view_mode(1) || view_mode(2))
    if(view_mode(1))
        % turn off features
        children=get(ax_features,'children');
        if(iscell(children)),children=cell2mat(children); end
        invisibles=[invisibles ax_features children(:)'];
    end
    
    % waveform plot stuff
    visibles=[visibles ax_waveforms(selected_channels)];
    
    wf_plots=getappdata(handles.output,'wf_plots');
    
    if(~isempty(wf_plots) && any(ishandle(wf_plots(:))))
        gridline_plots=getappdata(handles.output,'gridline_plots');
        errorbar_plots=getappdata(handles.output,'errorbar_plots');
        
        w1=wf_plots(~selected_channels,:);
        w2=wf_plots(:,~selected_groups);
        w3=gridline_plots(~selected_channels);
        w4=ch_title_plots(~selected_channels);
        invisibles=[invisibles w1(:)' w2(:)' w3 w4(:)'];
        
        w1=wf_plots(selected_channels,selected_groups);
        w2=gridline_plots(selected_channels);
        w3=ch_title_plots(selected_channels);
        visibles=[visibles w1(:)' w2 w3(:)'];
        
        
        
        if(getappdata(handles.output,'display_cleanup_tools'))
            w1=errorbar_plots(selected_channels,selected_groups);
            visibles=[visibles w1(:)'];
            
            w1=errorbar_plots(~selected_channels,:);
            w2=errorbar_plots(:,~selected_groups);
            
            invisibles=[invisibles w1(:)' w2(:)'];
        else
            invisibles=[invisibles errorbar_plots(:)'];
        end
    end
end

%   feature plots (highlight selected groups)
if(view_mode(2) || view_mode(3))
    if(view_mode(3))
        children=get(ax_waveforms,'children');
        if(iscell(children)),children=cell2mat(children); end
        invisibles=[invisibles ax_waveforms children(:)'];
    end
    
    children=get(ax_features,'children');
    if(iscell(children)),children=cell2mat(children); end
    visibles=[visibles ax_features children(:)'];
    
    sphere_plots=getappdata(handles.output,'sphere_plots');
    if(~getappdata(handles.output,'display_cleanup_tools'))
        invisibles=[invisibles sphere_plots(:)'];
    end
    
    highlight_selected_cluster(handles);
end

set(visibles,'visible','on');
set(invisibles,'visible','off');

if(view_mode(2) || view_mode(3))
    hide_clusters(handles);
end

% fix rotate3d
if(view_mode(1))
    set(handles.output,'windowbuttonmotionfcn',[]);
else
    % make sure we're in 3D view
    twodee=get(handles.feature3_menu,'Value')==1;
    if(~twodee)
        set(handles.output,'windowbuttonmotionfcn',{@enable_rotate,get(handles.features_panel,'position'),ax_features,getappdata(handles.output,'rotate_handle')});
    end
end