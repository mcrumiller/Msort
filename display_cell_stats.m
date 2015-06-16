function display_cell_stats(handles,update)
%DISPLAY_CELL_STATS   Display statistics plots of cell
%   DISPLAY_CELL_STATS shows an ISI histogram, a firing rate plot, and an autocorrelation plot of the selected neuron.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global Ts palette;

if(~exist('update','var') || isempty(update)), update=false; end

idx=getappdata(handles.output,'idx');

% do nothing if we haven't clustered yet
if(isempty(idx) || ~any(idx))
    set(handles.status_light,'BackgroundColor',[255 178 0]/255);
    set(handles.status_label,'String','No clusters yet!');
    return;
end

selected_axes=getappdata(handles.output,'selected_axes');
selected_axes(1)=[]; % remove noise axis
if(~any(selected_axes)),selected_axes=~selected_axes; end
% note: cluster 0 is noise, ignore it
u=unique(idx); u(u==0)=[];

margin_bot = 0.05; margin_top = .05;
margin_left = 0.18; margin_right = .01;
margin_middle_x=.02; margin_middle_y = .03;
plot_width=(1-2*margin_middle_x-margin_left-margin_right)/3;
plot_height = (1-margin_bot-margin_top - (length(u)-1)*margin_middle_y)/length(u);

%---plotting parameters
% ISI
%dt=1e-2; edges=0:dt:.5;
dt=1e-3; edges=0:dt:.1;
ISI_xtick=linspace(0,edges(end),5);
ISI_xticklabel=num2str(ISI_xtick');

% Rate Function
%span=Ts(end)-Ts(1);
resolution = 30; % in seconds
t_rate=unique(round(Ts(1)):resolution:round(Ts(end)));
xtick_ind=unique(round(linspace(1,length(t_rate),6)));
rate_xtick=t_rate(xtick_ind);
%rate_xtick=[round(Ts(1)) round(Ts(1)+span/4) round(Ts(1)+span/2) round(Ts(1)+3*span/4) round(Ts(end))];
rate_xticklabel=num2str(rate_xtick');
%locs=rate_xtick<Ts(1) | rate_xtick>Ts(end);
%rate_xtick(locs)=[];
%rate_xticklabel(locs)=[];



% Autocorrelation
lag=2e-2; % note: this is the total timespan of the autocorrelation function
num_bins=200; bin_width=2*lag/num_bins;
xspace=lag/2; autocorr_xtick=-lag:xspace:lag;
autocorr_xticklabel=num2str(autocorr_xtick');

% only make certain plots visible and leave; otherwise plot everything
if(update)
    ax=getappdata(handles.output,'ax_cell_stats');
else
    set(handles.status_light,'BackgroundColor',palette.red);
    set(handles.status_label,'String','Generating Cell Statistics...');

    timestamps=getappdata(handles.output,'timestamps');

    % grab spike trains
    num_clusters=length(u);
    spikes=cell(1,num_clusters);
    for i = 1:num_clusters
        spikes{i} = timestamps(idx==i);
    end

    % display plot
    plot_mean_waveforms(handles);

    ax=zeros(num_clusters,3);
    colors = getappdata(handles.output,'cluster_colors');
    if(size(colors,1)>num_clusters), colors(1,:)=[]; end
    for i = 1:num_clusters
        y = 1 - margin_top -(i-1)*margin_middle_y-i*plot_height;

        % ISI plot
        x = margin_left;
        ISI = diff(spikes{i});
        if(length(ISI)<2), n=0;
        else n=histc(ISI,edges); n=n./sum(n); n(isnan(n))=0; end
        ax(i,1) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        b=bar(edges,n,'histc'); set(b,'facecolor',colors(i,:),'edgecolor','none');
        min_n=0; max_n=max(n)*1.05;
        if(max_n==min_n), min_n=0; max_n=1; end
        axis([edges(1) edges(end) min_n max_n]);
        y_tick=linspace(min_n,max_n,5);
        set(ax(i,1),'ytick',y_tick,'yticklabel',[]);
        grid on;

        % put up text
        mean_rate = 1/mean(ISI);
        margin_x=(edges(end)-edges(1))*.01; margin_y=(max_n-min_n)*.05;
        text(edges(end)-margin_x,max_n-margin_y,sprintf('%2.1f spikes/s',mean_rate),'fontsize',10,...
            'fontweight','bold','color',colors(i,:),'horizontalalignment','right','verticalalignment','top','backgroundcolor',palette.ax_background);

        % rate plot
        x = margin_left + plot_width + margin_middle_x;
        %rate=smooth(histc(spikes{i},t),max(round(Ts(end)/100),1));
        rate = histc(spikes{i},t_rate)/resolution;
        ax(i,2) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        plot(t_rate,rate,'linewidth',2,'color',colors(i,:),'linesmooth','on');
        axis([round(Ts(1)) round(Ts(end)) 0 max(rate)*1.05]);
        ytick=(0:3)*max(rate)*1.05/4; ytick=unique(round(ytick*10)/10);
        set(ax(i,2),'ytick',ytick);
        grid on;

        % Autocorrelation plot
        s=spikes{i};
        x = margin_left + 2*plot_width+2*margin_middle_x;
        ax(i,3) = axes('Parent',handles.features_panel,'Position',[x y plot_width plot_height]); %#ok<LAXES>
        
        % only use up to the some # of spikes spikes
        % note that we can't randomly sample
        MAX_SPIKES=3e3;
        if(length(s)>MAX_SPIKES), s=s(1:MAX_SPIKES); end
        [t,n]=crosscorr(s,s,bin_width,lag);
        %loc=find(t==0);
        % remove center ripple
        %n([loc-1 loc loc+1]) = (n(loc-2)+n(loc+2))/2;
        n(t==0)=0;
        tt=linspace(t(1),t(end),1000);
        nn=interp1(t,n,tt,'pchip');
        plot(ax(i,3),tt,nn,'linesmooth','on','color',colors(i,:),'linewidth',2);
        grid on;
        maxy=max(nn)*1.05;
        ytick=0:maxy/4:maxy; ytick=round(ytick*10)/10;
        set(ax(i,3),'ytick',ytick);
        if(max(nn))==0,nn=1; end
        axis(ax(i,3),[-lag lag 0 max(nn)*1.05]);
    end
end
set(ax,'color',palette.ax_background,'fontweight','bold','box','on','linewidth',1.5,'xcolor',palette.text,'ycolor',palette.text);
c_off=get(ax(~selected_axes,:),'children');
if(iscell(c_off)), c_off=cell2mat(c_off); end
off=[reshape(ax(~selected_axes,:),1,[]) c_off(:)'];
c_on=get(ax(selected_axes,:),'children');
if(iscell(c_on)), c_on=cell2mat(c_on); end
on=[reshape(ax(selected_axes,:),1,[]) c_on(:)'];
set(off,'visible','off');
set(on,'visible','on');

% Deterine axis positions
num_selected=sum(selected_axes);
width=plot_width;
height = (1-margin_bot-margin_top - (num_selected-1)*margin_middle_y)/num_selected;
x1=margin_left;
x2=margin_left+plot_width+margin_middle_x;
x3=margin_left+2*(plot_width+margin_middle_x);
y=(0:num_selected-1)*(height+margin_middle_y)+margin_bot;
y=y(end:-1:1);
pos1=[ones(num_selected,1)*x1 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos2=[ones(num_selected,1)*x2 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos3=[ones(num_selected,1)*x3 y' ones(num_selected,1)*width ones(num_selected,1)*height];
pos=num2cell([pos1;pos2;pos3],2);
ax_old=ax;
ax=ax(selected_axes,:);
ax_selected=reshape(ax,1,[]);
set(ax_selected,{'position'},pos,'xcolor',palette.text,'ycolor',palette.text);

% generate text labels
title(ax(1,1),'ISI','fontsize',14,'fontweight','bold','color',palette.text);
title(ax(1,2),'Firing Rate','fontsize',14,'fontweight','bold','color',palette.text);
title(ax(1,3),'Autocorrelatoin','fontsize',14,'fontweight','bold');
xlabel(ax(end,1),'ISI (s)','fontsize',14,'fontweight','bold','color',palette.text);
xlabel(ax(end,2),'Time (s)','fontsize',14,'fontweight','bold','color',palette.text);
xlabel(ax(end,1),'Time (s)','fontsize',14,'fontweight','bold','color',palette.text);
set(ax(:,1),'xtick',ISI_xtick,'xticklabel',ISI_xticklabel);
set(ax(:,2),'xtick',rate_xtick,'xticklabel',rate_xticklabel);
set(ax(:,3),'xtick',autocorr_xtick,'xticklabel',autocorr_xticklabel);
set(ax(1:end-1,:),'xticklabel',[]);
set(ax,'xcolor',palette.text,'ycolor',palette.text);

% generate text labels
for i = 1:length(ax_selected)
    title(ax_selected(i),[]);
    xlabel(ax_selected(i),[]);
end

if(~update)
    % update selected panel
    setappdata(handles.output,'ax_cell_stats',ax_old);
    set(handles.status_label,'String','Done!');
    set(handles.status_light,'BackgroundColor',palette.green);
end