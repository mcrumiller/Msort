function feature_mask = get_feature_mask(handles)
%GET_FEATURE_MASK   Return mask for selecting out user-specified features
%   GET_FEATURE_MASK() Generates a logical mask which serves to select out the spike features requested by the user.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
% Determine which features to cluster along
PCs=get(handles.num_PCs_slider,'Value');
if(~iscell(PCs)),PCs={PCs}; end
PCs=cell2mat(PCs);

max_PCs=get(handles.num_PCs_slider,'Max');
if(~iscell(max_PCs)),max_PCs={max_PCs}; end
max_PCs=max([max_PCs{:}]);
num_channels=getappdata(handles.output,'num_channels');

% slope valley energy peak amplitude
feat_len=max_PCs+5;

feature_mask=false(feat_len,num_channels);

peaks=get(handles.peak_checkbox,'Value'); if(iscell(peaks)),peaks=cell2mat(peaks); end
valleys=get(handles.valley_checkbox,'Value'); if(iscell(valleys)),valleys=cell2mat(valleys); end
slopes=get(handles.slope_checkbox,'Value'); if(iscell(slopes)),slopes=cell2mat(slopes); end
energies=get(handles.energy_checkbox,'Value'); if(iscell(energies)),energies=cell2mat(energies); end
amplitudes=get(handles.amplitude_checkbox,'Value'); if(iscell(amplitudes)),amplitudes=cell2mat(amplitudes); end

for i=1:num_channels
    ind=1;
    num_PCs=PCs(i);
    feature_mask(1:num_PCs,i)=true; ind=ind+max_PCs;
    if(peaks(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(valleys(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(slopes(i)), feature_mask(ind,i)=true; end
    ind=ind+1; 
    if(energies(i)), feature_mask(ind,i)=true; end
    ind=ind+1;
    if(amplitudes(i)), feature_mask(ind)=true; end
end

selected_channels=getappdata(handles.output,'selected_channels');
if(~isempty(selected_channels) && any(selected_channels))
    feature_mask(:,~selected_channels)=false;
end
feature_mask=reshape(feature_mask,1,[]);