function features = extract_features(waveforms, t, params)
%--------------------------------------------------------------------------
% extract_features.m - Given a set of waveforms, extracts the requested
% features and returns them in a matrix.
%
% Usage: features = extract_features(waveforms);
%
% Input:  waveforms                * CxWxT matrix of waveforms (C channels,
%                                    W waveforms, T time points)
%         t                        * time series t
%         params (optional)        * parameter file (optional)
% Output: features                 * CxWxF matrix of feature projections.
%                                    Features are concatenated into the last
%                                    dimension.
%
% List of possible features:
%     - PCA
%     - initial slope
%     - peak-trough amplitude
%     - peak-trough time interval
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
%[num_channels, num_wf, num_pts] = size(waveforms);
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

% Extract each type of feature one by one
if(any(strcmp(feature_names,'valleys')) || any(strcmp(feature_names,'peaks')) || any(strcmp(feature_names,'slope')))
    
    
    
    % check 100us away from the valley
    Fs=abs(round(1/(t(2)-t(1))));
    
    x_start=50e-6; x_stop=100e-6; dx=x_stop-x_start;
    %dx=150e-6; % distance over which we calculate slope
    
    pts_start=round(Fs*x_start); pts_stop=round(Fs*x_stop);
    
    %[valleys,valleylocs]=min(waveforms,[],3);
    zeroloc=find_index(t,0);    
    valleys=squeeze(waveforms(:,:,zeroloc));
    slope_start=squeeze(waveforms(:,:,zeroloc+pts_start));
    slope_stop=squeeze(waveforms(:,:,zeroloc+pts_stop));
    slopes=(slope_stop-slope_start)./dx;
    
    max_peakspan=300e-6;
    pts_peak=round(Fs*max_peakspan);
    peaks=max(waveforms(:,:,zeroloc:zeroloc+pts_peak),[],3);
    
    energy=sum(waveforms.^2,3);
    
    
    
    % convert valleylocs to subind
    %I=repmat((1:size(waveforms,1))',1,size(waveforms,2));
    %J=repmat(1:size(waveforms,2),size(waveforms,1),1);
    %v2_locs=valleylocs+num_pts; v2_locs(v2_locs>length(t))=length(t);
    %ind=sub2ind(size(waveforms),I,J,v2_locs);
    
    %slopes=(peaks-valleys)./(t(peaklocs)-t(valleylocs));
    %slopes=(waveforms(ind)-valleys)./1e-4;
    nanlocs=isnan(slopes);
    slopes(nanlocs)=0;
end

for c = 1:num_channels
    feature_index=1;
    for f = 1:num_features
        switch(feature_names{f})	% select out the feature
            % run Principal Component Analysis
            case 'PCA'
                if(num_wf>1)
                    [~,scores] = princomp(squeeze(waveforms(c,:,:)));
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
        end
    end
end