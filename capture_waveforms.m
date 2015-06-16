function capture_waveforms(handles)
%CAPTURE_WAVEFORMS   Grab waveforms from voltage truce using selected timestamps
%   CAPTURE_WAVEFORMS Uses the current set of timestamps, as well as parameters t_before and t_after, to extract
%   waveforms of duration t_before+t_after centered.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global data_v Ts waveforms;

interp_factor=round(str2double(get(handles.interp_factor_field,'String')));

timestamps=getappdata(handles.output,'timestamps');
if(isempty(timestamps))
    waveforms=[];
    t=[]; %#ok<NASGU>
    return;
end
trigger_ch=getappdata(handles.output,'trigger_ch');
t_before=getappdata(handles.output,'t_before')*1e-6;
t_after=getappdata(handles.output,'t_after')*1e-6;
Fs=getappdata(handles.output,'Fs');
pts_before=round(Fs*t_before)+1;
pts_after=round(Fs*t_after)+1;

% convert timestamps to Ts locations
timestamp_locs = int32(round((timestamps-Ts(1))*Fs)+1);
badlocs=timestamp_locs<pts_before | timestamp_locs>length(Ts)-pts_after;
timestamp_locs(badlocs)=[]; timestamps(badlocs)=[]; trigger_ch(badlocs)=[];
N=length(timestamp_locs);
C = size(data_v,2);
L=pts_before+pts_after+1;
adder=repmat(int32(-pts_before:pts_after)',1,N);

index = repmat(timestamp_locs,L,1)+adder;
index = index(:);
waveforms = permute(reshape(data_v(index,:),[L N C]),[3 2 1]);

% interpolate waveforms
t=(-pts_before:pts_after)/Fs;
[w,t_interp]=interp_waveforms(waveforms,t,interp_factor);
    
% note: doing this may push two waveforms together.
[waveforms_interp,timestamps_interp,t_interp]=align_waveforms(w,timestamps,t_interp,interp_factor-1,trigger_ch);
badlocs=t_interp<-t_before | t_interp>t_after;
t_interp(badlocs)=[]; waveforms_interp(:,:,badlocs)=[];
normalized_waveforms=normalize_waveforms(waveforms);
normalized_waveforms_interp=normalize_waveforms(waveforms_interp);
waveforms=waveforms(:,:,2:end-1); t=t(2:end-1);
normalized_waveforms=normalized_waveforms(:,:,2:end-1);

% update application data
setappdata(handles.output,'original_waveforms',waveforms);
setappdata(handles.output,'original_timestamps',timestamps);
setappdata(handles.output,'original_normalized_waveforms',normalized_waveforms);
setappdata(handles.output,'original_t',t);
setappdata(handles.output,'waveforms_interp',waveforms_interp);
setappdata(handles.output,'timestamps_interp',timestamps_interp);
setappdata(handles.output,'t_interp',t_interp);
setappdata(handles.output,'normalized_waveforms_interp',normalized_waveforms_interp);

% set current data based on check box
if(get(handles.interp_waveforms_checkbox,'Value'))
    waveforms=waveforms_interp;
    normalized_waveforms=normalized_waveforms_interp;
    t=t_interp;
    timestamps=timestamps_interp;
end
setappdata(handles.output,'waveforms',waveforms);
setappdata(handles.output,'normalized_waveforms',normalized_waveforms);
setappdata(handles.output,'timestamps',timestamps);
setappdata(handles.output,'t',t);