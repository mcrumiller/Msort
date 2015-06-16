function complete_line(hObject,~,ax,cluster_colors)
%COMPLETE_LINE   Complete drawing current line upon mouse recently
%   COMPLETE_LINE When user releases mouse, determines where the line was and selects all waveforms. It then draws
%   "highlighted" waveforms overtop to indicate which waveforms were highlighted. Note: hidden waveforms may appear.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms;
global locsX locsY h ctrl;

if(~ctrl)
    handles.output=hObject;
    deselect_button_Callback([],[],handles);
end

set(getappdata(hObject,'ax_waveforms'),...
    'units','normalized');
set(hObject,'Units','Normalized','windowButtonMotionFcn',[],'windowButtonUpFcn',[]);

delete(h);

% grab waveforms
t=getappdata(hObject,'t')*1e6;
cluster_numbers=getappdata(hObject,'cluster_number');

% determine which axis was selected
ax_waveforms=getappdata(hObject,'ax_waveforms');
channel = ax_waveforms==ax;

% get time range
t1=find(t<min(locsX),1,'last'); t2=find(t>max(locsX),1,'first');
if(isempty(t1)),t1=1; end
if(isempty(t2)),t2=length(t); end
t_range=t1:t2;

% get only waveforms in these clusters
idx=getappdata(hObject,'idx');
idx=repmat(idx,length(cluster_numbers),1);
cluster_numbers=repmat(cluster_numbers',1,size(idx,2));
wf_locs=any(idx==cluster_numbers,1);

wfs=squeeze(waveforms(channel,wf_locs,t_range));
if(length(t_range)>1 && isvector(wfs)),wfs=wfs(:)'; end
num_wfs=size(wfs,1);
t_new=t(t_range); t_length=length(t_new);

% prepare line segments for intersection test
XY2 = [locsX(1) locsY(1) locsX(2) locsY(2)];

wfs2=cat(3,wfs(:,1:end-1),wfs(:,2:end));
wfs2=permute(wfs2,[3 2 1]);
wfs2=reshape(wfs2,2,[])';
t2=repmat([t_new(1:end-1)' t_new(2:end)'],num_wfs,1);

XY1=[t2(:,1) wfs2(:,1) t2(:,2) wfs2(:,2)];
intersection=lineSegmentIntersect(XY1,XY2);
intersection=reshape(intersection,t_length-1,[]);
intersection=any(intersection,1);

if(~any(intersection))
    return;
end
% get actual waveform locations
locs=find(wf_locs);
locs=locs(intersection);

if(ctrl)
    old_locs=getappdata(hObject,'selected_waveforms');
else
    old_locs=false(1,length(idx));
end
if(length(old_locs)~=size(waveforms,2)),old_locs=false(1,length(idx)); end
old_locs(locs)=true; locs=old_locs;

% draw new selected waveforms
selected_waveforms=waveforms(:,locs,:);

selected_channels=getappdata(hObject,'selected_channels');
if(~any(selected_channels)),selected_channels=~selected_channels; end

% remove old plots
selected_wf_plots=getappdata(hObject,'selected_wf_plots');
delete(selected_wf_plots(ishandle(selected_wf_plots) & selected_wf_plots~=0));

% scroll through each waveform axis, add selected waveform lines
for i=1:length(ax_waveforms);
    wfs=selected_waveforms(i,:,:);
    num_wfs=size(wfs,2);
    wfs=squeeze(wfs);
    if(num_wfs==1),wfs=wfs'; end
    wfs=[wfs';nan(1,num_wfs)];
    wfs=reshape(wfs,1,[]);
    
    t_tmp=repmat([t nan],1,num_wfs);
    
    if(selected_channels(i))
        selected_wf_plots(i)=plot(ax_waveforms(i),t_tmp,wfs,'linesmooth','on','color','k','hittest','off','visible','on');
    else
        selected_wf_plots(i)=plot(ax_waveforms(i),t_tmp,wfs,'linesmooth','on','color','k','hittest','off','visible','off');
    end
end
drawnow;

setappdata(hObject,'selected_wf_plots',selected_wf_plots);
setappdata(hObject,'selected_waveforms',locs);

cluster_numbers=getappdata(hObject,'cluster_number');
num_clusters=getappdata(hObject,'num_clusters');
locs=1:num_clusters; locs(cluster_numbers+1)=[];
if(~isempty(locs))
    merge_buttons=getappdata(hObject,'merge_buttons');
    merge_buttons=merge_buttons(locs);
    set(merge_buttons,'enable','on');
    cluster_colors=cluster_colors(locs,:);
    for i = 1:length(merge_buttons)
        set(merge_buttons(i),'BackgroundColor',cluster_colors(i,:));
    end
end