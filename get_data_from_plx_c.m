function [data, Ts, gaps] = get_data_from_plx_c(filename, channels, timespan)
%--------------------------------------------------------------------------
% get_spikes_from_plx_c.m - Given a .plx file, grabs the continuous voltage
% trace.  Returns a 1xT vector of timestamps associated with each voltage.
%
% Usage: [waveforms, timestamps] = get_spikes_from_plx_c(filename, channels);
%
% Input:  filename                     * .plx file containing continuous-
%                                        channel recordings
%         channels (optional)          * vector of channels to use.  Uses
%                                        all channels if none specified.
%         timespan (optional)          * time span of data to collect
% Output: data                         * CxT vector of data values,
%                                        specified in amplified volts
%         Ts                           * CxT vector
%         gaps                         * 2xG matrix of gap times. First row
%                                        is starts; second row is stops.
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
num_channels = plx_chan_names(filename);
if(~exist('channels','var'))
    channels=1:num_channels;
end
%channels=unique(channels);
if(num_channels<max(channels))
    error('Invalid Channel: maximum channel specified in file is Channel %g',num_channels);
end

[~,samplecounts]=plx_adchan_samplecounts(filename);
%samplecounts=samplecounts(samplecounts~=0);
%if(~all(~diff(samplecounts))),
%warning('Channels are not aligned or some were disabled partway through the recording.\n');
%    fprintf('Warning: some samples were disabled partway through the recording.\n');
%end
max_samplecounts=max(samplecounts);

[~,adchannels]=plx_ad_chanmap(filename);
channels=adchannels(channels); % convert to zero-based channels
num_channels=length(channels);

% make sure we get the real frequency--sometimes channels are missing
adfreq=0; i=1;
while(adfreq==0 && i<=length(channels))
    [adfreq, ~, t_start] = plx_ad_gap_info(filename, channels(i));
    i=i+1;
end
% we have no data in any channels
if(adfreq==0)
    data=0; Ts=0;
    return;
end
dt=1/adfreq;

% Grab frequency and time count info
if(exist('timespan','var') && ~isempty(timespan))
    % Determine indexing and time vector
    pts_before = max(0,floor((timespan(1)-t_start)*adfreq));
    num_pts = ceil((timespan(2)-timespan(1))*adfreq)+2;
    startCount=pts_before+1;
    endCount=startCount+num_pts-1;
    start_time=(t_start+dt*pts_before);
    Ts=start_time:dt:start_time+dt*(num_pts-1);
    
    % Grab continuous data from file
    % Note: if the data is discontinuous, this might return incorrect results
    % It might be better to get the entire dataspan and chunk it together using
    % t_frag
    data = single(zeros(num_pts,num_channels),'int16');
    for i = 1:num_channels
        [~, ~, data(:,i)] = plx_ad_span(filename, channels(i), startCount,endCount);
        % if it's all -1, continue
        if(all(~(data(:,i)+1)))
            data(:,i)=0;
        end
    end
    data=data';
else
    data = zeros(num_channels,max_samplecounts,'int16');
    for i = 1:num_channels
        current_samplecount=samplecounts(channels(i)+1);
        [~, ~, ts, fn, data(i,1:current_samplecount)] = plx_ad(filename, channels(i));
        
        % we have to return gap information, in case the PLX recording was
        % paused partway through
        gap_start = (ts+dt*fn+dt)';
        gap_stop = [ts(2:end)' ts(end)+fn(end)*dt];
        gaps = [gap_start;gap_stop];
    end
    Ts = (0:max_samplecounts-1)*dt+ts(1);
    
    % note: the output does not include gaps
end