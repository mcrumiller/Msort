function plot_filtered_voltage_trace(filename, channels, timespan, bandpass, INTERP_FACTOR)
%--------------------------------------------------------------------------
% plot_filtered_voltage_traces.m - given a .plx file and a channel, plots
% the specified span of time.
%
% Usage: plot_voltage_traces(filename, channels, timespan);
%
% Input:  filename                * .plx file
%         channels                * 1xC vector of channels to plot
%         timespan (optional)     * 1x2 vector of start/stop times
%         bandpass (optional)     * bandpass range
% Output: <matlab plot>
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
if(~exist('bandpass','var') || isempty(bandpass))
    bandpass = [150 10000];
end

% Grab data and time points
if(exist('timespan','var'))
    if(flag)
        [data,Ts] = get_data_from_plx_c(filename,channels,[timespan(1) timespan(2)+5]);
    else
        [data,Ts] = get_data_from_plx_c(filename,channels,timespan);
    end
else
    [data,Ts] = get_data_from_plx_c(filename,channels);
end
if(isempty(Ts)), return; end
if(~exist('timespan','var')), timespan=[Ts(1) Ts(end)]; end

sampling_rate=1/(mode(diff(Ts)));

% filter data results
data = filter_voltage(data,sampling_rate,bandpass);
if(exist('INTERP_FACTOR','var') && INTERP_FACTOR>1)
	sampling_rate=sampling_rate*INTERP_FACTOR;
	Ts2=Ts(1):1/sampling_rate:Ts(end);
	data=interp1(Ts,data,Ts2,'pchip');
	Ts=Ts2;
end

if(sum(data==-1)==length(data)), return; end

% Plot results
if(isempty(get(0,'currentfig'))), fig_std; hold on; end
plot(Ts,data','k');
max_v = max(max(abs(data)))*1.01;
axis([timespan(1) timespan(2) -max_v max_v]);
xlabel('Time (s)'); ylabel('Amplified volts');