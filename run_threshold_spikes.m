function run_threshold_spikes(handles)
%RUN_THRESHOLD_SPIKES   Wrapper function for actual thresholding of spikes
%   RUN_THRESHOLD_SPIKES retrieves thresholding parameters and run sthe thresholding routine, then extracts waveforms.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global waveforms features;
set(handles.status_label,'String','Applying thresholds...'); drawnow;

% rethreshold spikes
t_before=getappdata(handles.output,'t_before')*1e-6;
t_after=getappdata(handles.output,'t_after')*1e-6;

thresholds=getappdata(handles.output,'thresholds');

% threshold spikes
%if(get(handles.square_voltage_checkbox,'Value'))
%    data_squared=getappdata(handles.output,'data_squared');
%    [timestamps,trigger_ch] = M_threshold_spikes(handles, thresholds,t_before,t_after,data_squared);
%else
    [timestamps,trigger_ch] = M_threshold_spikes(handles, thresholds,t_before,t_after);
%end
if(isempty(timestamps))
    set(handles.status_label,'String','No Waveforms found; reduce threshold.');
    set(handles.status_light,'backgroundcolor',[255 178 0]/255);
    waveforms=[];
    setappdata(handles.output,'timestamps',[]);
    return;
end

% aligning may have aligned identical spikes; remove duplicates
setappdata(handles.output,'timestamps',timestamps);
setappdata(handles.output,'trigger_ch',trigger_ch);
capture_waveforms(handles);

num_waveforms=size(waveforms,2);
idx=zeros(1,num_waveforms,'uint8');
features=[];
setappdata(handles.output,'trigger_ch',trigger_ch);

setappdata(handles.output,'idx',idx);
setappdata(handles.output,'idx_old',idx);
fix_clusters(handles);