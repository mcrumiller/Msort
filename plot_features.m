function plot_features(handles)
%PLOT_FEATURES   Generate feature plot
%   PLOT_FEATURES(handles) Updates the features axis with a plot dependent on current plot parameters.
%   
%   Input:  handles..................vector of GUI object handles
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global rotate features palette;
rotate = false;
if(isempty(features))
    calculate_features(handles);
end
timestamps=getappdata(handles.output,'timestamps');
features_temp=[features timestamps(:)];

set(handles.status_light,'backgroundcolor',palette.red);
set(handles.status_label,'String','Plotting features...'); drawnow;

% show features axis
ax_features=getappdata(handles.output,'ax_features');
feat_children=get(ax_features,'children'); if(iscell(feat_children)),feat_children=cell2mat(feat_children); end
delete(feat_children);
set(ax_features,'XTickMode','auto','YTickMode','auto','xlimmode','manual','ylimmode','manual','zlimmode','manual');

% Determine which features we have selected
selected_features=getappdata(handles.output,'selected_features');
if(isempty(selected_features))
    feature1=get(handles.feature1_menu,'Value');
    feature2=get(handles.feature2_menu,'Value');
    feature3=get(handles.feature3_menu,'Value');
    setappdata(handles.output,'selected_features',[feature1 feature2 feature3]);
else
    set(handles.feature1_menu,'Value',selected_features(1));
    set(handles.feature2_menu,'Value',selected_features(2));
    set(handles.feature3_menu,'Value',selected_features(3));
    feature1=selected_features(1);
    feature2=selected_features(2);
    feature3=selected_features(3);
end
feature1_list=get(handles.feature1_menu,'String');
feature2_list=get(handles.feature2_menu,'String');
feature3_list=get(handles.feature3_menu,'String');
cluster_colors=getappdata(handles.output,'cluster_colors');
num_clusters=getappdata(handles.output,'num_clusters');
feature_plots=zeros(1,num_clusters);

feature_names=getappdata(handles.output,'feature_names');

%feature_names=get(handles.feature1_menu,'String');


sphere_plots=zeros(1,num_clusters);
num_stds=get(handles.errorbar_slider,'Value');
wfs_in_sphere=false(size(features,1),num_clusters);
show_spheres=getappdata(handles.output,'display_cleanup_tools');
if(show_spheres),visible='on'; else visible='off'; end
idx=getappdata(handles.output,'idx');

if(feature2==length(feature2_list) || feature3==length(feature3_list))
    hist_mode=true;
    optimize_width=get(handles.hist_optimize_checkbox,'Value');
else
    hist_mode=false; hist_plots=[];
end
setappdata(handles.output,'hist_mode',hist_mode);
num_pts=size(features_temp,1);

% We have a 2D plot
if(feature3==1)
    % create histogram instead of scatterplot
    if(hist_mode)
        minFeat=min(features_temp(:,feature1)); maxFeat=max(features_temp(:,feature1));
        % optimize bin width
        if(~optimize_width)
            num_bins=round(20*log10(num_pts));
            edges=linspace(minFeat,maxFeat,num_bins);
            % double the X to make histogram
            x=[reshape(repmat(edges,2,1),1,[]) edges(end-1:-1:1)];
        end
        
        
        max_N=0;
        hist_plots=zeros(1,num_clusters);
        
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            feat1=features_temp(ind,feature1);
            if(optimize_width)
                edges=best_bin_width(feat1);
                edge_diff=edges(2)-edges(1);
                num_bins=ceil((maxFeat-minFeat)/edge_diff);
                edges=linspace(minFeat,maxFeat,num_bins);
                % double the X to make histogram
                x=[reshape(repmat(edges,2,1),1,[]) edges(end-1:-1:1)];
            end
            N=histc(feat1,edges); N(end)=[]; max_N=max(max_N,max(N));
            N2=reshape(repmat(N',2,1),1,[]); N2=[0 N2 0]; %#ok<AGROW>
            hist_plots(i)=patch(x,[N2 zeros(1,length(N))],cluster_colors(i,:),'linesmooth','on','parent',ax_features,'edgecolor','none','FaceAlpha','flat');
        end
        h=ylabel(ax_features,'Count','color',palette.text,'fontsize',16,'fontweight','bold');
        span=edges(end)-edges(1);
        set(ax_features,'ButtonDownFcn',[],'xlim',[edges(1)-span*.01 edges(end)+span*.01],'ylim',[0 max_N*1.05]);
        alpha(hist_plots,.6);
    else
        MaxF=[-inf -inf]; MinF=[inf inf];
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1);
            feat2=features_temp(ind,feature2);
            feats=[feat1 feat2];
            MaxF=max([MaxF;max(feats,[],1)]);
            MinF=min([MinF;min(feats,[],1)]);
            
            % calculate ellipse
            feature_plots(i)=plot(ax_features,feat1,feat2,'.','markersize',1,'linesmooth','on','color',cluster_colors(i,:));
            if(i~=1 && size(feats,1)>2)
                [Ex,Ey]=ellipsoid_cov(mean(feats),num_stds,cov(feats),30);
                sphere_plots(i)=patch(Ex,Ey,cluster_colors(i,:),'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'parent',ax_features,'visible',visible);
                wfs_in_sphere(ind,i)=inhull(feats,get(sphere_plots(i),'Vertices'));
            end
        end
        alpha(sphere_plots(ishandle(sphere_plots) & sphere_plots~=0),0.2);
        set(ax_features,'xlim',[MinF(1) MaxF(1)],'ylim',[MinF(2) MaxF(2)],'ButtonDownFcn',{@select_features,handles});
        h=ylabel(ax_features,feature_names{feature2},'color',palette.text,'fontsize',16,'fontweight','bold');
    end
    view(ax_features,2); set(ax_features,'projection','orthographic');
    xlabel(ax_features,feature_names{feature1},'color',palette.text,'fontsize',16,'fontweight','bold');
    ylim=get(ax_features,'ylim');
    pos=get(h,'position'); set(h,'position',[pos(1) (ylim(2)-ylim(1))/3 pos(3)]);
    set(handles.output,'WindowButtonMotionFcn',[]);

% we have a 3D plot
else
    feature3=feature3-1;
    if(hist_mode)
        if(feature2==length(feature2_list)),feature2=feature1; end

        minFeat1=min(features_temp(:,feature1)); maxFeat1=max(features_temp(:,feature1));
        minFeat2=min(features_temp(:,feature2)); maxFeat2=max(features_temp(:,feature2));
        del = .001; % space between bars, relative to bar size
        if(~optimize_width)
            num_binsX=round(20*log10(num_pts)); num_binsY=num_binsX;
            edgesX=linspace(minFeat1,maxFeat1,num_binsX);
            edgesY=linspace(minFeat2,maxFeat2,num_binsY);
            %dX=edgesX(2)-edgesX(1); dY=edgesY(2)-edgesY(1);
            
            % Build x-coords for the eight corners of each bar.
            dX = edgesX(2)-edgesX(1);
            xx = [edgesX edgesX(end)+dX];
            xx = [xx(1:num_binsX)+del*dX; xx(2:num_binsX+1)-del*dX];
            xx = [reshape(repmat(xx(:)',2,1),4,num_binsX); NaN(1,num_binsX)];
            xx = [repmat(xx(:),1,4) NaN(5*num_binsX,1)];
            xx = repmat(xx,1,num_binsX);
            
            % Build y-coords for the eight corners of each bar.
            dY = edgesY(2)-edgesY(1);
            yy = [edgesY edgesY(end)+dY];
            yy = [yy(1:num_binsY)+del*dY; yy(2:num_binsY+1)-del*dY];
            yy = [reshape(repmat(yy(:)',2,1),4,num_binsY); NaN(1,num_binsY)];
            yy = [repmat(yy(:),1,4) NaN(5*num_binsY,1)];
            yy = repmat(yy',num_binsY,1);
        end
        
        max_N=0;
        hist_plots=zeros(1,num_clusters);
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2);
            if(optimize_width)
                edgesX=best_bin_width(feat1); edgesY=best_bin_width(feat2);
                edge_diffX=edgesX(2)-edgesX(1);
                edge_diffY=edgesY(2)-edgesY(1);
                num_binsX=ceil((maxFeat1-minFeat1)/edge_diffX);
                num_binsY=ceil((maxFeat2-minFeat2)/edge_diffY);
                edgesX=linspace(minFeat1,maxFeat1,num_binsX);
                edgesY=linspace(minFeat2,maxFeat2,num_binsY);
                % Build x-coords for the eight corners of each bar.
                dX = edgesX(2)-edgesX(1);
                xx = [edgesX edgesX(end)+dX];
                xx = [xx(1:num_binsX)+del*dX; xx(2:num_binsX+1)-del*dX];
                xx = [reshape(repmat(xx(:)',2,1),4,num_binsX); NaN(1,num_binsX)];
                xx = [repmat(xx(:),1,4) NaN(5*num_binsX,1)];
                xx = repmat(xx,1,num_binsY);

                % Build y-coords for the eight corners of each bar.
                dY = edgesY(2)-edgesY(1);
                yy = [edgesY edgesY(end)+dY];
                yy = [yy(1:num_binsY)+del*dY; yy(2:num_binsY+1)-del*dY];
                yy = [reshape(repmat(yy(:)',2,1),4,num_binsY); NaN(1,num_binsY)];
                yy = [repmat(yy(:),1,4) NaN(5*num_binsY,1)];
                yy = repmat(yy',num_binsX,1);
            end
            n=hist3([feat1 feat2],{edgesX edgesY},'edgecolor','none','parent',ax_features);
            max_N=max(max_N,max(n(:)));
            % Build z-coords for the eight corners of each bar.
            
            zz = zeros(5*num_binsX, 5*num_binsY);
            n(n==0)=NaN;
            zz(5*(1:num_binsX)-3, 5*(1:num_binsY)-3) = n;
            zz(5*(1:num_binsX)-3, 5*(1:num_binsY)-2) = n;
            zz(5*(1:num_binsX)-2, 5*(1:num_binsY)-3) = n;
            zz(5*(1:num_binsX)-2, 5*(1:num_binsY)-2) = n;
            hist_plots(i)=surf(ax_features, xx, yy, zz,'facecolor',cluster_colors(i,:),'edgecolor','none','tag','hist3','linesmooth','off','facelighting','flat','Facealpha',0.6);
        end
        
        %alpha(hist_plots,.6);
        set(ax_features,'xlim',[edgesX(1) edgesX(end)],'ylim',[edgesY(1) edgesY(end)],'zlim',[0 max_N],'box','off');
    else
        MaxF=[-inf -inf -inf]; MinF=[inf inf inf];
        for i = 1:num_clusters
            ind=idx==i-1;
            if(~any(ind)),continue; end
            
            feat1=features_temp(ind,feature1); feat2=features_temp(ind,feature2); feat3=features_temp(ind,feature3);
            feats=[feat1 feat2 feat3];
            MaxF=max([MaxF;max(feats,[],1)]);
            MinF=min([MinF;min(feats,[],1)]);
            feature_plots(i)=plot3(ax_features,feat1,feat2,feat3,'.','markersize',1,'color',cluster_colors(i,:),'buttondownfcn',[],'hittest','off');
            
            % calculate ellipsoid
            if(size(feats,1)>2 && feature1 ~= feature2 && feature1 ~= feature3 && feature2 ~= feature3 && i~=1)
                [Ex,Ey,Ez]=ellipsoid_cov(mean(feats),num_stds,cov(feats),30);
                % sometimes we get a weird result with complex numbers with
                % zero imaginary component
                if(any(~isreal(Ex))),continue; end
                sphere_plots(i)=surf(Ex,Ey,Ez,'linesmooth','on','edgecolor','none','facecolor',cluster_colors(i,:),'facelighting','gouraud','parent',ax_features,'visible',visible);
                [~,v]=surf2patch(Ex,Ey,Ez);
                wfs_in_sphere(ind,i)=inhull(feats,v);
            end
        end
        set(ax_features,'xlim',[MinF(1) MaxF(1)],'ylim',[MinF(2) MaxF(2)],'zlim',[MinF(3) MaxF(3)],'ButtonDownFcn',[]);
        zlabel(feature3_list{feature3},'fontsize',16,'fontweight','bold','color',palette.text);
        alpha(sphere_plots(ishandle(sphere_plots) & sphere_plots~=0),0.2);
    end
    
    xlabel(feature1_list{feature1},'fontsize',16,'fontweight','bold','color',palette.text);
    ylabel(feature2_list{feature2},'fontsize',16,'fontweight','bold','color',palette.text);
    
    view(ax_features,[15 45]);
    h=light('parent',ax_features);
    camlight(h,45,45);
    set(ax_features,'projection','perspective');
    
    h=rotate3d(ax_features);
    setappdata(handles.output,'rotate_handle',h);
    set(handles.output,'WindowButtonMotionFcn',{@enable_rotate,get(handles.features_panel,'position'),ax_features,h});
end
set(ax_features,'color',palette.ax_background,'xcolor',palette.gridlines,'ycolor',palette.gridlines,'fontweight','bold','zcolor',palette.gridlines);
grid(ax_features,'on');

setappdata(handles.output,'feature_plots',feature_plots);
setappdata(handles.output,'sphere_plots',sphere_plots);
setappdata(handles.output','wfs_in_sphere',wfs_in_sphere);
setappdata(handles.output,'hist_plots',hist_plots);

% highlight selected
highlight_selected_cluster(handles);
hide_clusters(handles);

set(handles.status_label,'String','Done!');
set(handles.status_light,'backgroundcolor',palette.green);