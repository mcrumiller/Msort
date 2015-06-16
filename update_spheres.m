function update_spheres(handles)
%UPDATE_SPHERES   Change size of ellipsoids around neuron clusters
%   UPDATE_SPHERES() Updates the 2- or 3-D ellipsoids surrounding clusters. The sizes are determined by the user's
%   settings in the "cleanup" box.
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
global features;
sphere_plots=getappdata(handles.output,'sphere_plots');
if(isempty(sphere_plots) || ~any(ishandle(sphere_plots)))
    return;
end
wfs_in_sphere=getappdata(handles.output,'wfs_in_sphere');
radius=get(handles.errorbar_slider,'Value');

% Determine which features we have selected
feature1=get(handles.feature1_menu,'Value');
feature2=get(handles.feature2_menu,'Value');
feature3=get(handles.feature3_menu,'Value');

if(any([feature1 feature2 feature3]==size(features,2)+1))
    features_temp=[features getappdata(handles.output,'timestamps')'];
else
    features_temp=features;
end

idx=getappdata(handles.output,'idx');
num_clusters=getappdata(handles.output,'num_clusters');
% 2-dimensional case
if(feature3==1)
    for i = 2:num_clusters
        ind=idx==i-1;
        if(~any(ind)),continue; end
        
        feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2);
        feats=[feat1 feat2];
        
        % calculate ellipse
        [Ex,Ey]=ellipsoid_cov(mean(feats),radius,cov(feats),30);
        %sphere_plots(i)=patch(Ex,Ey,cluster_colors(i,:),'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'parent',ax_features,'visible','off');
        set(sphere_plots(i),'XData',Ex','YData',Ey');
        wfs_in_sphere(ind,i)=inhull(feats,get(sphere_plots(i),'Vertices'));
    end
else
    feature3=feature3-1;
    for i = 2:num_clusters
        ind=idx==i-1;
        if(~any(ind)),continue; end
        feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2); feat3=features_temp(ind,feature3);
        feats=[feat1 feat2 feat3];
        
        % calculate ellipsoid
        if(size(feats,1)>2 && feature1 ~= feature2 && feature1 ~= feature3 && feature2 ~= feature3)
            [Ex,Ey,Ez]=ellipsoid_cov(mean(feats),radius,cov(feats),30);
            set(sphere_plots(i),'XData',Ex,'YData',Ey,'ZData',Ez);
            %sphere_plots(i)=surf(Ex,Ey,Ez,'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'facelighting','gouraud','parent',ax_features,'visible','off');
            [~,v]=surf2patch(Ex,Ey,Ez);
            wfs_in_sphere(ind,i)=inhull(feats,v);
        end
    end
end
setappdata(handles.output,'sphere_plots',sphere_plots);
setappdata(handles.output,'wfs_in_sphere',wfs_in_sphere);