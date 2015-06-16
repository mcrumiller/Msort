function M_extract_features(t, params, handles, PCA_type)
%M_EXTRACT_FEATURES   Extract features from waveforms
%   M_EXTRACT_FEATURES extracts from set of waveform from the features selected by the user in the features panel.
%
%   Input:  waveforms..................CxWxT matrix of waveforms (C channels, W waveforms, T time points)
%           t..........................time series t
%           params (optional)..........parameter file (optional)
%
%   Output: features...................CxWxF matrix of features. Features are concatenated into the last dimension.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global features waveforms;
[num_channels, num_wf, ~] = size(waveforms);
% Deal with parameter file
if(~exist('params','var') || isempty(params))
    params = struct('PCA',3,'slope',0,'pt_amplitude',0,'pt_interval',0);
else
    if(ischar(params))
        params = load_params(params);
    end
end
feature_names=fieldnames(params);
locs = structfun(@(x) x==0,params);
params = rmfield(params,feature_names(locs));
feature_names(locs)=[];
num_features = length(feature_names);
if(any(strcmp('PCA',feature_names)))
    num_PCs=params.PCA;
    feature_length=num_features-1+params.PCA;
else
    num_PCs=0;
    feature_length=num_features;
end

features = zeros(num_channels,num_wf, feature_length);
idx=getappdata(handles.output,'idx');
switch PCA_type
    case 'all'
        selection = true(1,num_wf);
    case 'active'
        selection = idx~=0;
    case 'clusters'
        cluster_number=getappdata(handles.output,'cluster_number');
        u=unique(cluster_number);
        selection=false(1,num_wf);
        for i = 1:length(u)
            selection(idx==u(i))=true;
        end
end

% Extract each type of feature one by one
if(any(strcmp(feature_names,'valleys')) || any(strcmp(feature_names,'peaks')) || any(strcmp(feature_names,'slope')))
    % check 100us away from the valley
    Fs=abs(round(1/(t(2)-t(1))));
    
    x_start=str2double(get(handles.slope_start_field,'String'))*1e-6;
    x_stop=str2double(get(handles.slope_stop_field,'String'))*1e-6;
    dx=x_stop-x_start;
    
    pts_start=round(Fs*x_start); pts_stop=round(Fs*x_stop);
    
    zeroloc=find_index(t,0);    
    valleys=squeeze(waveforms(:,:,zeroloc));
    slope_start=squeeze(waveforms(:,:,zeroloc+pts_start));
    slope_stop=squeeze(waveforms(:,:,zeroloc+pts_stop));
    slopes=(slope_stop-slope_start)./dx;
    
    max_peakspan=300e-6;
    pts_peak=round(Fs*max_peakspan);
    peaks=max(waveforms(:,:,zeroloc:zeroloc+pts_peak),[],3);
    
    energy=sum(waveforms.^2,3);

    nanlocs=isnan(slopes);
    slopes(nanlocs)=0;
end

for c = 1:num_channels
    feature_index=1;
    if(any(strcmp(feature_names,'amplitude')))
        wfs=squeeze(waveforms(c,:,:));
    end
    for f = 1:num_features
        switch(feature_names{f})	% select out the feature
            % run Principal Component Analysis
            case 'PCA'
                % if only active waveforms, we have to run PCA on a subset
                % and then project waveforms onto the PCs
                
                if(num_wf>1)
                    wfs=squeeze(waveforms(c,:,:));
                    % note: this overrides using "active waveforms"
                    if(strcmp(PCA_type,'moving'))
                        timestamps=getappdata(handles.output,'timestamps');
                        window_size=round(str2double(get(handles.PCA_window_field,'String')));
                        [~,scores]=princomp_time(wfs,timestamps,window_size,'pchip');
                    else
                        [coeff,~]=pca(wfs(selection,:));
                        wfs_nomean=wfs-repmat(mean(wfs,1),size(wfs,1),1);
                        scores=wfs_nomean*coeff;
                    end
                else
                    scores=zeros(num_wf,params.PCA);
                end
                
                features(c,:,feature_index+(0:num_PCs-1))=scores(:,1:num_PCs);
                feature_index=feature_index+num_PCs;
            % Maximum point
            case 'peaks'
                features(c,:,feature_index)=peaks(c,:);
                feature_index=feature_index+1;
            % Minimum point
            case 'valleys'
                features(c,:,feature_index)=valleys(c,:);
                feature_index=feature_index+1;
            % Slope from peak to trough
            case 'slope'
                features(c,:,feature_index)=slopes(c,:);
                feature_index=feature_index+1;
            % integral of squared waveform
            case 'energy'
                features(c,:,feature_index)=energy(c,:);
                feature_index=feature_index+1;
            % maximum - minimum span
            case 'amplitude'
                w=squeeze(waveforms(c,:,:));
                features(c,:,feature_index)=max(wfs,[],2)-min(wfs,[],2);
        end
    end
end