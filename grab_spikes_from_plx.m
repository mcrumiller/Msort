function [waveforms, timestamps, t] = grab_spikes_from_plx(filename, channels, threshold, interp_factor, timespan)
%--------------------------------------------------------------------------
% get_spikes_from_plx.m - Given a .plx file, grabs continuous channel data,
% interpolates it, and filters it.  Spike extraction is then performed, and
% channels are grouped if desired (i.e. voltage traces will be extracted on
% all channels in a group if just one initiates a superthreshold event).
%
% Usage: [waveforms, timestamps, t] = grab_spikes_from_plx(filename, ...
%     channels, threshold, interp_factor, groups);
%
% Input:  filename           * .plx file
%         channels           * list of channels to grab
%         threshold          * threshold for spike extraction in std devs
%         interp_factor      * interpolate data by interp_factor
%         timespan (optional)* timespan of data
% Output: waveforms          * CxWxT matrix of waveforms, zero-padded
%         timestamps         * CxW matrix of spike times
%         t                  * 1xT array of time span of any given spike
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
if(~exist('channels','var')), channels=1; end
if(~exist('threshold','var')), threshold=2; end
if(~exist('interp_factor','var')), interp_factor=4; end
if(~exist('timespan','var')), timespan=[]; end

% extract waveform information
[data,Ts] = get_data_from_plx_c(filename,channels,timespan);

if(isempty(Ts)), waveforms=[]; timestamps=[]; return; end

if(exist('timespan','var') && ~isempty(timespan))
    if(length(timespan)==1), timespan=[0 timespan]; end
    Tsmin=find_index(Ts,timespan(1));
    Tsmax=find_index(Ts,timespan(2));
    Ts=Ts(Tsmin:Tsmax);
    data=data(:,Tsmin:Tsmax);
end
dt=mode(diff(Ts));
freq=1/dt;

% filter data results
data = filter_voltage(data,freq);

% grab spikes
t_before = 300e-6; t_after = 500e-6;
[waveforms, timestamps, t] = threshold_spikes(data, Ts, threshold, t_before-dt, t_after+dt);

% Interpolate spikes
waveforms=interp_waveforms(waveforms,t,interp_factor);
waveforms=permute(waveforms,[3 1 2]);
t2=linspace(t(1),t(end),(length(t)-1)*interp_factor+1);
waveforms=interp1(t,waveforms,t2,'spline');
waveforms=permute(waveforms,[2 3 1]);
t=t2;

% align waveforms
[waveforms, timestamps, t]=align_waveforms(waveforms,timestamps,t,interp_factor);