function edges = best_bin_width(X, min_N, max_N, jump)
%BEST_BIN_WIDTH determines bin width for 
% best_bin_width.m - Determines the optimal histogram bin width of a vector
% X.
%
%
% Usage: bin_width = best_bin_width(X, min_binsize, max_binsize);
%
% Input:  X                    * 1xN vector of data
%         min_N (optional)     * minimum number of bins
%         max_N (optional)     * maximum numbe rof bins
%         jump (optional)      * test every jumpth bin number
% Output: num_bins             * optimal number of bins
%
% Reference:
%    Neural Comput. 2007 Jun;19(6):1503-27.
%    A method for selecting the bin size of a time histogram.
%    Shimazaki H1, Shinomoto S.
%    http://www.ncbi.nlm.nih.gov/pubmed/17444758
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
X=X(:)';

if(~exist('min_N','var') || isempty(min_N)), min_N=5; end
if(~exist('max_N','var') || isempty(max_N)), max_N=500; end
if(~exist('jump','var') || isempty(jump)), jump=1; end

test_N=round(min_N):jump:round(max_N);
minval=min(X); maxval=max(X);
D=(maxval-minval)./test_N;
C=zeros(1,length(test_N));

% create histograms for each bin width
for i = 1:length(test_N);
    N=test_N(i);
    edges=linspace(minval,maxval,N);
    ki=histc(X,edges); ki(end)=[];
    k=mean(ki);
    v=sum((ki-k).^2/N);
    C(i)=(2*k-v)/D(i)^2; % cost function
end

% optimal bin width selection
[~,idx]=min(C);
edges=linspace(minval,maxval,test_N(idx));