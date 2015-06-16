function templates = generate_templates(waveforms,idx)
%--------------------------------------------------------------------------
% generate_templates.m - Given a set of waveforms, generates templates by
% taking the mean and normalizing to unity.
%
% Usage: templates = generate_templates(waveforms);
%
% Input:  waveforms                 * CxWxT matrix of waveforms
%         idx                       * 1xW vector of cluster IDs
% Output: templates                 * GxCxT matrix of waveform templates
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[C,~,T]=size(waveforms);
u=unique(idx); u(u==0)=[]; G=length(u);
zerolocs=idx==0;
waveforms(:,zerolocs,:)=[]; idx(zerolocs)=[];
mean_wfs=get_mean_waveforms(waveforms,idx,false); % false means don't include zero
mean_wfs=permute(mean_wfs,[3 2 1]);
mean_wfs=reshape(mean_wfs,C*T,G);
lens=sqrt(sum(mean_wfs.^2,1));
lens=repmat(lens,C*T,1);
mean_wfs=mean_wfs./lens;
mean_wfs=reshape(mean_wfs,[T C G]);
templates=permute(mean_wfs,[3 2 1]);