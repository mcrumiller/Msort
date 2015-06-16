function varargout = get_mean_waveforms(waveforms, idx, include_zero)
%GET_MEAN_WAVEFORMS Return a GxCxWxT matrix of mean waveforms,
% one for each unique group.  Waveforms is a CxWxT matrix of waveforms
% captured from thresholding (output by threshold_spikes.m), and idx is the
% 1xW index matrix that identifies each waveform with a cluster group.
%
% Usage:  waveforms = get_mean_waveforms(waveforms, idx);
%         [waveforms,std_dev] = get_mean_waveforms(waveforms, idx);
%
% Input:  waveforms                   * CxWxT matrix of waveforms
%         idx                         * 1xW vector of waveform group IDs
%                                       Note: IDs = [1, 2, 3, ...]
%         include_zero                * 'idx' count starts at 0, not 1
% Output: mean_wfs                    * GxCxT matrix of mean waveforms
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
unique_groups = unique(idx);
if(exist('include_zero','var') && include_zero)
    if(~any(unique_groups==0))
        unique_groups=[0 unique_groups];
    end
else
    unique_groups(unique_groups==0)=[];
end
G=length(unique_groups);
[C,~,T] = size(waveforms);

mean_wfs = zeros(G,C,T,'single'); % single precision preserves memory

if(exist('include_zero','var') && include_zero)
    wfs=waveforms(:,idx==0,:);
	if(~isempty(wfs)), mean_wfs(1,:,:) = permute(mean(wfs,2),[2 1 3]); end
    offset=1;
else
    offset=0;
end

G=max(unique(unique_groups));

for i = 1:G
	wfs = waveforms(:,idx==i,:);
	mean_wfs(i+offset,:,:) = permute(mean(wfs,2),[2 1 3]);
end

varargout{1}=mean_wfs;

% calculate distance as well
if(nargout==2)
    waveforms=permute(waveforms,[3 1 2]);
    waveforms=reshape(waveforms,T*C,[]);
    mean_wfs=permute(mean_wfs,[3 2 1]);
    mean_wfs=reshape(mean_wfs,T*C,[]);
    if(nargout>1)
        std_dev=zeros(1,G);
        for i = 1:G
            ind=idx==i; num_wfs=sum(ind);
            wfs=waveforms(:,ind);
            mw=repmat(mean_wfs(:,i),1,num_wfs);
            d=sqrt(sum((mw-wfs).^2));
            % note: we want the average deviation from zero, which is equal to the mean
            std_dev(i)=mean(d);
        end
        varargout{2}=std_dev;
    end
end