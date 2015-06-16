function index = remove_outliers(waveforms, idx)
%--------------------------------------------------------------------------
% remove_outliers.m - Uses PCA for each cluster to determine outliers to be
% removed. A Gaussian PCA projection is assumed, and 6 standard deviations
% are used as the cutoff.
%
% Usage: [waveforms,ind_removed] = remove_outliers(waveforms, idx);
%
% Input:  waveforms                  * CxWxT matrix of waveforms
%         idx                        * 1xW vector of cluster membership IDs
% Output: index                      * 1xW logical vector of waveforms to
%                                      remove.
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
cutoff = 4; % in std devs

[C,W,T]=size(waveforms);
waveforms=reshape(permute(waveforms,[3 1 2]),C*T,[])';
index=false(W,1);

num_groups=length(unique(idx));
for i = 1:num_groups
    locs=find(idx==i);
    wfs=waveforms(locs,:);
    [~,score]=princomp(wfs);
    score=zscore(score(:,[1 2]));
    bad_cells = any(score<-cutoff | score>cutoff,2);
    index(locs(bad_cells))=true;
end