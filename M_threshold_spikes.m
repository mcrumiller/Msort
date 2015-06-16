function [timestamps, trigger_ch] = M_threshold_spikes(handles, threshold, t_before, t_after, abs_data)
%M_THRESHOLD_SPIKES   Threshold voltage traces and extract spikes
%   [timestamps, trigger_ch] = M_THRESHOLD_SPIKES is the main thresholding routine. This function uses preset threshold
%   parameters to determine spike times from a discretized, continuous voltage trace. Waveforms span from t_before to
%   t_after the spike peak; the peaks are returned as timestamps. The channel on which each associated spike is
%   triggered is returned in TRIGGER_CH.
%
%   Note that the actual waveforms themselves are not extracted until CAPTURE_WAVEFORMS() is called.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global data_v Ts palette;

num_channels = size(data_v,2);

if(~exist('threshold','var')), threshold=4; end
if(~exist('t_before','var')), t_before = 300e-6; end
if(~exist('t_after','var')), t_after=500e-6; end

dt=mean(diff(Ts));
freq=round(1/dt);

% Update: We've since removed this and added it in later
% determine deadzone span for a single channel
same_ch_deadzone=str2double(get(handles.same_ch_deadzone_field,'String'))*1e-6;
diff_ch_deadzone=str2double(get(handles.diff_ch_deadzone_field,'String'))*1e-6;
same_ch_deadzone_span=round(same_ch_deadzone/dt);
diff_ch_deadzone_span=round(diff_ch_deadzone/dt);

% run through each channel individually
spikelocs=cell(1,num_channels);
peaks = cell(1,num_channels);
if(nargout>1)
    h=waitbar(0,sprintf('Thresholding %g/%g...',0,num_channels),'color',palette.background);
    set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
end
for i = 1:num_channels
    peaks1=[]; peaks2=[];
    if(nargout>1),waitbar(i/num_channels,h,sprintf('Thresholding %g/%g...',i,num_channels));end
    % top threshold first
    if(threshold(1,i)~=0)
        if(exist('abs_data','var'))
            [peaks1,spikelocs1]=findpeaks2(abs_data(:,i)','minpeakheight',threshold(1,i));
        else
            [peaks1,spikelocs1]=findpeaks2(data_v(:,i)','minpeakheight',threshold(1,i));
        end
    else spikelocs1=[];
    end
    
    % bottom threshold
    if(threshold(2,i)~=0)
        [peaks2,spikelocs2]=findpeaks2(-data_v(:,i)','minpeakheight',-threshold(2,i));
    else spikelocs2=[];
    end
    
    % if we have a string of spikes that are within the deadzone, choose
    % the spike with the greatest amplitude
    [spikelocs{i},ind]=unique([spikelocs1 spikelocs2]);
    if(length(spikelocs{i})>2)
        peaks{i}=[peaks1 peaks2]; peaks{i}=abs(peaks{i}(ind));
        badlocs=diff(spikelocs{i})<same_ch_deadzone_span;
        d=diff([false badlocs]);
        starts=find(d==1); stops=find(d==-1);
        if(length(stops)<length(starts)),stops=[stops starts(end)+1]; end %#ok<AGROW>
        ranges=stops-starts;
        maxpeaks=zeros(1,length(starts));
        for j = 1:length(starts)
            ind=starts(j):starts(j)+ranges(j);
            [~,maxpeaks(j)]=max(peaks{i}(ind));
        end
        maxpeaks=maxpeaks+starts-1;

        b = [false badlocs] | [badlocs false]; b(maxpeaks)=false;
        spikelocs{i}(b)=[];
    end
end
if(nargout>1),waitbar(1,h,'Merging results...'); end

% if we have a string of spikes that are within the cross-channel deadzone,
% choose the spike with the greatest amplitude
spike_nums=cellfun(@length,spikelocs);
tc=cell(1,num_channels);
for i = 1:num_channels
    tc{i}=ones(1,spike_nums(i))*i;
end
s=[spikelocs{:}]; p=[peaks{:}]; trigger_ch=[tc{:}];
[spikelocs,ind]=unique(s); peaks=p(ind); trigger_ch=trigger_ch(ind);
if(length(spikelocs)>2)
    badlocs=diff(spikelocs)<diff_ch_deadzone_span;
    d=diff([false badlocs]);
    starts=find(d==1); stops=find(d==-1);
    if(length(stops)<length(starts)),stops=[stops starts(end)+1]; end
    ranges=stops-starts;
    maxpeaks=zeros(1,length(starts));
    for i = 1:length(starts)
        ind=starts(i):starts(i)+ranges(i);
        [~,maxpeaks(i)]=max(peaks(ind));
    end
    maxpeaks=maxpeaks+starts-1;
    
    b = [false badlocs] | [badlocs false]; b(maxpeaks)=false;
    spikelocs(b)=[];
    trigger_ch(b)=[];
end


% determine number of points per wave
pts_before=round(t_before*freq);
pts_after=round(t_after*freq);

% remove spikes at beginning and end
locs=spikelocs<=pts_before | spikelocs>length(Ts)-pts_after;
spikelocs(locs)=[]; trigger_ch(locs)=[];

% get actual spike times
timestamps = Ts(spikelocs);

delete(h);