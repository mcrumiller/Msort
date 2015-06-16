function [snippets, t] = grab_snippets_from_mat(mat_file, channels, timestamps, t_before, t_after, show_progress)
%--------------------------------------------------------------------------
% grab_snippets_from_plx.m - Given a series of timestamps, grabs the
% specified snippets from plx_file from each channel.  The length of the
% snippet is determined by pts_before and pts_after, which specify the
% number of index points before and after the timestamp to be captured.
%
% Usage: grab_snippets_from_plx(plx_file, timestamps, channels, ...
%    pts_before, pts_after);
%
% Input:  mat_file                     * .mat file
%         num_channels                 * 1xC number of channels
%         pts_before                   * number of index points before
%         pts_after                    * number of index points after
%         timespan (optional)          * timespan of data
%         show_progress (optional)     * display progress bar
% Output: snippets                     * CxWxT matrix of waveform snippets
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
num_wfs = length(timestamps);
num_channels=length(channels);
% load Ts variable
load(mat_file,'Ts_all');
Fs=round(1/(Ts_all(2)-Ts_all(1)));

% to minimize memory usage, process each channel separately
pts_before = round(t_before*Fs); pts_after = round(t_after*Fs);
t=(-pts_before:pts_after)/Fs;
num_pts = pts_before+pts_after+1;
snippets = zeros(num_channels,num_wfs,num_pts,'single');

exp_start = Ts_all(1);
ts_ind = int32(round((timestamps-exp_start)*Fs));
ts_ind(ts_ind<=pts_before)=[];
%ts_ind(ts_ind>=num_points-pts_after)=[];
ts_locs = repmat(ts_ind,num_pts,1);
adder = int32(repmat((-pts_before+1:pts_after+1)',1,num_wfs));
ts_locs=ts_locs+adder; %#ok<NASGU>

% cycle through each channel, grab the snippets
if(show_progress)
    h=waitbar(0,sprintf('Processing Channel 0/%g',num_channels),'name','Grabbing Snippets');
end
for i = 1:num_channels
    if(show_progress)
        waitbar(i/num_channels,h,sprintf('Processing Channel %g/%g',i,num_channels));
        drawnow;
    end
    load(mat_file,sprintf('data_all_%g',channels(i)));
    
    % grab this channel's snippets
    eval(sprintf('snippets(i,:,:) = data_all_%g(ts_locs)'';',channels(i)));
end
if(show_progress), delete(h); end