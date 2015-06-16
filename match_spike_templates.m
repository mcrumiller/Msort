function [idx,dists,confidence] = match_spike_templates(templates, waveforms)
%--------------------------------------------------------------------------
% match_spike_templates.m - Given a set of spike templates and standard
% deviations, calculates the distance of each spike to a particular
% template, assigns it to a cluster, and returns the distance with respect
% to that cluster.
%
% Note: this function does not use adaptive templates.
%
% Usage: [idx,zscores] = match_spike_templates(templates,std_devs,spikes);
%
% Input:  templates           * GxCxT set of spike templates (T templates,
%                               C channels, W time points)
%         std_devs            * 1xT vector of distance standard deviations
%                               for each template
%         waveforms           * CxMxW matrix of M spikes
% Output: idx                 * 1xM vector of M cluster assignments
%         dists               * 1xM vector of distances to each cluster
%         confidence          * 1xM vector of confidence values (0 is high
%                               confidence, 1 is low confidence)
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[G,C,T]=size(templates);
W=size(waveforms,2);
% no templates!
if(G==0)
    idx=zeros(1,W);
    dists=idx; confidence=idx;
    return;
end

% determine amplitudes of each channel
%{
amps=abs(min(templates,[],3));
amps = amps./repmat(sum(amps,2),1,size(amps,2));
amps=repmat(amps,[1 1 T]);
amps=permute(amps,[3 2 1]);
amps=reshape(amps,T*C,G)';
amps=repmat(amps,[1 1 W]);
%}

% determine max 2 channels for each waveform
if(C>1)
    [~,ch]=sort(min(templates,[],3),2,'ascend');
    ch=ch(:,1:2);
end

% set up matrices
waveforms=permute(waveforms,[3 1 2]);
waveforms=reshape(waveforms,T*C,[]);
templates=permute(templates,[3 2 1]);
templates=reshape(templates,T*C,[]);

% normalize waveforms and templates to vectors of unit length
wf_len=sqrt(sum(waveforms.^2));
template_len=sqrt(sum(templates.^2));
waveforms=waveforms./repmat(wf_len,C*T,1);
templates=templates./repmat(template_len,C*T,1);
templates=repmat(templates,[1 1 W]);

% Calculate the distances to each template
distances=zeros(W,G);
%for i = 1:G
%    distances(:,i) = sqrt(sum(  ((waveforms-squeeze(templates(:,i,:))).*squeeze(amps(i,:,:)))  .^2));
%end

for i = 1:G
    % set up templates to only include the top channel
    if(C>1)
        locs=false(T,C); locs(:,ch(i,:))=1;
        distances(:,i) = sqrt(sum((waveforms(locs,:)-squeeze(templates(locs,i,:))).^2));
    else
        distances(:,i) = sqrt(sum((waveforms-squeeze(templates(:,i,:))).^2));
    end
end

min_dist=repmat(min(distances,[],2),1,G);
confidence = min_dist./distances;
confidence = (sum(confidence,2)-1)/(G-1);

% put distances in terms standard deviations
[~,idx]=min(distances,[],2);
I=sub2ind(size(distances),(1:W)',idx);
dists=distances(I);
idx=idx';