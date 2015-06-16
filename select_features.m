function select_features(~, ~, handles)
%SELECT_FEATURES   Select features with mouse
%   [] = SELECT_FEATURES() 
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global features palette;
ax_features=getappdata(handles.output,'ax_features');


% Determine which features have been plotted
feature1=get(handles.feature1_menu,'Value');
feature2=get(handles.feature2_menu,'Value');
feature3=get(handles.feature3_menu,'Value')-1;

timestamps=getappdata(handles.output,'timestamps');
features_temp=[features timestamps'];

new_ax=axis(ax_features);
clusters_hidden=getappdata(handles.output,'clusters_hidden');
if(length(new_ax)==6)
    features_temp=features_temp(:,[feature1 feature2 feature3]);
    % replot features in 2D
    T=view(ax_features);
    % generate new axis limits
    new_ax=T*[reshape(axis(ax_features)',2,[]) ones(2,1)]';
    new_ax=reshape(new_ax(1:2,:)',1,[]);
    XYZ=T*[features_temp ones(size(features_temp,1),1)]';
    
    % we might have to do some data/axis flipping
    % flip X values down the middle
    if(new_ax(1)>new_ax(2))
        mean_x=mean(XYZ(1,:));
        d=mean_x-XYZ(1,:);
        XYZ(1,:)=mean_x+d;
        tmp=new_ax(1); new_ax(1)=new_ax(2); new_ax(2)=tmp;
    end
    
    % flip Y axis
    if(new_ax(3)>new_ax(4))
        mean_y=mean(XYZ(2,:));
        d=mean_y-XYZ(2,:);
        XYZ(2,:)=mean_y+d;
        tmp=new_ax(3); new_ax(3)=new_ax(4); new_ax(4)=tmp;
    end
    
    
    % put new axis temporarily over old axis
    delete(get(ax_features,'children'));
    p=plot(ax_features,XYZ(1,:),XYZ(2,:),'.','markersize',1,'linesmooth','on','color',palette.text);
    view(ax_features,0,90);
    axis(ax_features,new_ax);
    set(ax_features,'color',palette.background,'xtick',[],'ytick',[]);
    set(handles.output,'currentaxes',ax_features);
    
    % allow user to draw circle
    h=imfreehand(ax_features);
    locs=h.getPosition; delete(h);
    selected = inpolygon(XYZ(1,:)',XYZ(2,:)',locs(:,1),locs(:,2));
    
% 2 dimensional plot
else    
    if(clusters_hidden)
        idx=getappdata(handles.output,'idx');
        ind=false(1,length(idx));
        cluster_number=getappdata(handles.output,'cluster_number');
        for i = 1:length(cluster_number)
            ind = ind | idx==cluster_number(i);
        end
        features_temp=features_temp(ind,[feature1 feature2]);
    else
        features_temp=features_temp(:,[feature1 feature2]);
    end
    %set(handles.output,'currentaxes',ax_features);
    
    ax=axis(ax_features);
    delete(get(ax_features,'children'));
    p=plot(ax_features,features_temp(:,1),features_temp(:,2),'.','markersize',1,'linesmooth','on','color',palette.text);
    axis(ax_features,ax);
    %set(handles.output,'currentaxes',ax_features);
    
    % allow user to draw circle around features
    h=imfreehand(ax_features);
    locs=h.getPosition; delete(h);
    selected=inpolygon(features_temp(:,1),features_temp(:,2),locs(:,1),locs(:,2));
end

% highlight selected points
colors=[palette.text; palette.red];
X=get(p,'Xdata'); Y=get(p,'Ydata');
color_ind=ones(1,length(X));
color_ind(selected)=2;

delete(p);
colormap(ax_features,colors);
scatter(ax_features,X,Y,1,color_ind');

if(clusters_hidden)
    s=1:length(idx); s=s(ind); s(~selected)=[];
    selected=false(1,length(idx)); selected(s)=true;
end

setappdata(handles.output,'selected_waveforms',selected);

% enable merge buttons
merge_buttons=getappdata(handles.output,'merge_buttons');
set(merge_buttons,'enable','on');
cluster_colors=getappdata(handles.output,'cluster_colors');
for i = 1:length(merge_buttons)
    set(merge_buttons(i),'backgroundcolor',cluster_colors(i,:));
end