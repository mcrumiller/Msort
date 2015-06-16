function key_press(~,eventData,handles)
%KEY_PRESS   Intercept keyboard commands
%   KEY_PRESS intercepts keyboard commands and runs the appropriate function.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global ctrl window_size palette;

% toggle ctrl key on and off
if(strcmp(eventData.Key,'control'))
    ctrl=true;

% toggle hide non-selected clusters from features panel
elseif(strcmp(eventData.Key,'h'))
    clusters_hidden=getappdata(handles.output,'clusters_hidden');
    if(clusters_hidden)
        setappdata(handles.output,'clusters_hidden',false);
        fix_featurepanel_display(handles)
        set(handles.hide_clusters_label,'visible','on');
    else
        setappdata(handles.output,'clusters_hidden',true);
        fix_featurepanel_display(handles)
        set(handles.hide_clusters_label,'visible','off');
    end
    
% This is the print screen command
elseif(strcmp(eventData.Key,'p'))
	% determine what's on the screen
    panel_selected=getappdata(handles.output,'panel_selected');
    
    % create blank figure
    scrsz = get(0,'ScreenSize');
    width = scrsz(3)*.9; height = scrsz(4)*.9;
    left=(scrsz(3)-width)*.5; bot=(scrsz(4)-height)*.5;
    f=figure('outerposition',[left bot width height],'color',palette.background,'PaperPositionMode','Auto','InvertHardCopy','off','visible','off');
    switch panel_selected
        % voltage trace
        case 'voltage'
            voltage_plot=getappdata(handles.output,'voltage_plot');
            % need to adjust voltage plot a little after we copy
            c=copyobj(voltage_plot,f);
            pos=get(c,'position');
            pos(1)=.05; pos(4)=.94; set(c,'position',pos);
            
            ch_voltage_tags=getappdata(handles.output,'ch_voltage_tags');
            c=copyobj(ch_voltage_tags,f);
            pos=get(c,'position');
            if(iscell(pos)),pos=cell2mat(pos); end
            pos(:,3)=.02;
            set(c,{'position'},num2cell(pos,2));
            
            % if black and white is checked, convert all axis children to black,
            % make plot white
            if(get(handles.screenshot_bw_checkbox,'Value'))
                ch=get(f,'children');
                control_locs=findobj(f,'type','uicontrol');
                ax_locs = findobj(f,'type','axes');
                ch2=get(ax_locs,'children');

                % set colors
                set(f,'color',[1 1 1]);
                set(ax_locs,'color',[1 1 1],'ycolor',palette.text,'xcolor',palette.text);
                set(control_locs,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0]);

                % % determine which is the voltage plot, set to black
                Ydata=get(ch2,'YData');
                L=cellfun(@length,Ydata');
                longlocs=L==max(L);
                set(ch2(longlocs),'color','k');

                attributes=arrayfun(@get,ch2,'UniformOutput',false);
                has_color = cellfun(@(x) isfield(x,'BackgroundColor'),attributes);
            end
            
        % features
        case 'features'
            view_mode=logical([get(handles.waveform_view_button,'value') get(handles.split_view_button,'value') get(handles.features_view_button,'value')]);
            wf_axes=getappdata(handles.output,'wf_axes');
            
            % stretch wf_axes to fill full height
            ax_waveforms=getappdata(handles.output,'ax_waveforms');
            ax_features=getappdata(handles.output,'ax_features');
            
            pos=get(ax_features,'position'); margin_top=1-(pos(2)+pos(4));
            pos=get(ax_waveforms,'position');
            if(iscell(pos)), pos=cell2mat(pos); margin_bot=min(pos(:,2));
            else margin_bot=pos(1);
            end
            total_height=1-(margin_top+margin_bot);
            num_ax=length(wf_axes);
            plot_height=total_height/num_ax;
            y=(0:num_ax-1)*plot_height;
            pos=get(wf_axes(1),'position');
            pos=repmat(pos,num_ax,1);
            pos(:,2)=margin_bot+y';
            pos(:,4)=plot_height;
            pos=pos(end:-1:1,:);
            c=copyobj(wf_axes,f);
            set(c,{'position'},num2cell(pos,2));

            % waveforms
            if(view_mode(1) || view_mode(2))
                copyobj(ax_waveforms,f);
            end
            
            % features
            if(view_mode(2) || view_mode(3))
                copyobj(ax_features,f);
            end
            
            % shrink text
            text_locs = findobj(f,'type','text');
            fontsize=get(text_locs,'fontsize');
            fontsize=[fontsize{:}];
            set(text_locs',{'fontsize'},num2cell(fontsize'*.75));
            
            
            % if black and white is checked, convert all axis children to black,
            % make plot white
            if(get(handles.screenshot_bw_checkbox,'Value'))
                ch=get(f,'children');
                control_locs=findobj(f,'type','uicontrol');
                ax_locs = findobj(f,'type','axes');
                ch2=get(ax_locs,'children');

                
                % set colors
                set(f,'color',[1 1 1]);
                set(ax_locs,'color',[1 1 1],'ycolor',[0 0 0],'xcolor',[0 0 0],'linewidth',1);
                set(control_locs,'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0]);
                
                set(text_locs,'BackgroundColor',[1 1 1],'color',[0 0 0]);
            end
            
        
        % cell states
        case 'cell_stats'
            ax_cell_stats=getappdata(handles.output,'ax_cell_stats');
            wf_axes=getappdata(handles.output,'wf_axes');
            copyobj(ax_cell_stats,f);
            c2=copyobj(wf_axes,f);
            
            cellstat_pos=cell2mat(get(ax_cell_stats,'position'));
            left=min(cellstat_pos(:,1)); bottom=min(cellstat_pos(:,2)); top=max(cellstat_pos(:,2)+cellstat_pos(:,4));
            
            num_cells=length(c2);
            total_height=top-bottom; plot_height=total_height/num_cells;
            bottom_wf_ax=[.01 bottom left-.02 plot_height];
            pos=repmat(bottom_wf_ax,num_cells,1);
            pos(:,2)=pos(:,2)+((num_cells-1:-1:0)*plot_height)';
            
            set(c2,{'position'},num2cell(pos,2));
    end
    
    % export figure
    listing=ls('Msort_*.png');
    list_length=size(listing,1);
    filename=sprintf('Msort_%02.0f.png',list_length+1);
    set(f,'visible','on');
    print(f,'-dpng','-r300',filename,'-painters');
    close(f);
    set(handles.status_label,'String',sprintf('Exported %s/%s.',strrep(pwd,'\','/'),filename));
    
% zoom in on voltage plot
elseif(strcmp(eventData.Key,'equal') && strcmp(eventData.Modifier,'shift'))
    Fs=getappdata(handles.output,'Fs');
    min_window_size=round(.01*Fs); % .01 seconds max
    old_window_size=window_size;
    window_size=max(round(old_window_size*.8),min_window_size);
    voltage_slider=getappdata(handles.output,'voltage_slider');
    
    % determine new SliderStep
    slider_max=get(voltage_slider,'Max');
    slider_min=get(voltage_slider,'Min');
    slider_span=slider_max-slider_min;
    stepSize=[1 4]*window_size/(slider_span*5);
    set(voltage_slider,'SliderStep',stepSize);
    
    % re-center view
    pos=get(voltage_slider,'Value');
    middle=pos+old_window_size/2-window_size/2;
    set(voltage_slider,'Value',middle);
    
    voltage_slider_Callback(getappdata(handles.output,'voltage_slider'),[],handles);
% zoom out on voltage plot
elseif(any(strcmp(eventData.Key,'hyphen')) && any(strcmp(eventData.Modifier,'shift')));
    Fs=getappdata(handles.output,'Fs');
    max_window_size=round(5*Fs); % 5 seconds max
    old_window_size=window_size;
    window_size=min(round(old_window_size*1.25),max_window_size);
    voltage_slider=getappdata(handles.output,'voltage_slider');
    
    % determine new SliderStep
    slider_max=get(voltage_slider,'Max');
    slider_min=get(voltage_slider,'Min');
    slider_span=slider_max-slider_min;
    stepSize=[1 4]*window_size/(slider_span*5);
    set(voltage_slider,'SliderStep',stepSize);
    
    % re-center view
    pos=get(voltage_slider,'Value');
    middle=pos+old_window_size/2-window_size/2;
    set(voltage_slider,'Value',middle);
    
    voltage_slider_Callback(getappdata(handles.output,'voltage_slider'),[],handles);
end

function ctrl_release(~,eventData)
global ctrl;
if(strcmp(eventData.Key,'control'))
    ctrl=false;
end