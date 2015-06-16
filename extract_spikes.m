function [times,waveforms] = extract_spikes(data,Hz,Wn,thresh,plot_on)

%
%   Extracts spikes from continuous data in the form samples x channels
%
%       Inputs: data - in the form samples x channels
%               Hz - sampling frequency in Hz
%               Wn - low and high cuts ([300 6000])
%               thresh - spike detection threshold; can be a scalar, or 
%                         a different threshold can be specified for each
%                         channel
%
%       Outputs: times - spike times relative to start of data
%                waveforms - ms x n x m matrix of waveforms
%                            ms = number of points in 1 ms of data
%                            n = number of channels
%                            m = number of spikes
%
%     code by Josh Siegle (jsiegle@mit.edu)
%

%plot_on = 0;



%N = 3*round((1./min(Wn))./(1./Hz));

%b = fir1(N,Wn./(Hz./2));

%s_filt = filtfilt(b,1,data);
s_filt=data;

if plot_on
    plot(linspace(0,length(data)./Hz,length(data)),s_filt,'k')
end

thresh_cross = [ ];
channel = [ ];

%thresh = 4*median(abs(s_filt)./0.6745); % from Quian Quiroga et al. 2004

%disp(thresh)

for n = 1:size(s_filt,2)
    
    if numel(thresh) > 1
        f = find(s_filt(:,n) > thresh(n));
    else
        f = find(s_filt(:,n) > thresh);
    end
    
    jumps = find(diff(f) > 1);
    jump_start = [1; jumps+1];
    try
    thresh_cross = [thresh_cross; f(jump_start)];
    channel = [channel; ones(size(jump_start)).*n];
    go_on = 1;
    catch
       go_on = 0;
    end
    
end

times = [ ];
waveforms = [ ];


if go_on
% remove overlapping spike inds

[spike_inds,i] = unique(thresh_cross);
channel = channel(i);
[spike_inds,i] = sort(spike_inds);
channel = channel(i);

% remove spike inds separated by less than 1 ms

onems = round(Hz./1000);
f = find(diff(spike_inds) > onems);
spike_inds2 = spike_inds(f);
channel = channel(f);

pt1 = round(onems.*0.1);
pt9 = round(onems.*0.9);
pt3 = round(onems.*0.3);
pt7 = round(onems.*0.7);

% remove first and last 5 spikes!


for m = 5:numel(spike_inds2)-5
    
    candidate_spike = s_filt(spike_inds2(m)-pt1:spike_inds2(m)+pt9,channel(m));
    
    [max_val,peak_ind] = max(candidate_spike);
    
    peak_index(m) = peak_ind - pt1;
    
    times(m-4) = (spike_inds2(m)+peak_index(m)-1)./Hz;
    %disp(peak_index);

    waveforms(:,:,m-4) = s_filt(spike_inds2(m)+peak_index(m)-pt3:spike_inds2(m)+peak_index(m)+pt7,:);

end

    
%times = (spike_inds2(5:end-5)+peak_index(m)-1)./Hz;
if plot_on
    hold on
plot(times,squeeze(max(waveforms)),'.','color',rgb('purple'))
hold off
end

end