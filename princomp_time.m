function [coeff, score] = princomp_time(X, t, window_size, method, num_PCs)
%PRINCOMP_TIME   
%   [coeff,score] = princomp_time(X, t, window_size) returns a
%   time-dependent Principal Component Analysis on input data X.
%   
%   The n-by-p data matrix X is split into chunks determined by window_size
%   and the resulting eigenvectors are spline- nterpolated to create smooth
%   functions of time over which the data is projected.
%
%   [coeff,score] = princomp_time(X, t, window_size, method) allows
%   for an additional interpolation argument passed to the INTERP1 function
%   See the interp1 documentation for additional details.
%   
%   Input:  X                     * NxP data matrix. Rows correspond to
%                                   observations, columns to variables.
%           t                     * 1xT vector of associated timestamps for
%                                   each of the N observations.
%           window_size           * length in seconds of each PCA window.
%           method (optional)     * interpolation method used by INTERP1.
%   
%   Output: coeff                 * Eigenvectors return as functions of
%                                   time.
%           score
%
%   Written by Marshall Crumiller
%   email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[M,N]=size(X);
t_start=t(1); t=t-t_start;
if(M~=length(t))
    error('Time vector must equal number of observations');
end

% check input data
if(~exist('method','var') || isempty(method)), method = 'spline'; end
if(~exist('window_size','var') || isempty(window_size))
    % divide into ten even chunks, assuming we have at least 30 data points
    % per chunk.
    if(M<300), T=floor(M/30);
    else T = ceil(M/10);
    end
else
    T=floor(t(end)/window_size);
end

% determine window time midpoints
t_window = (.5:1:T-.5)*window_size;

coeff2=zeros(N,N,T);
score = zeros(size(X));

locs = zeros(2,T);
means2 = zeros(N,T);

for i = 1:T
    locs(1,i)=find(t>=window_size*(i-1),1,'first');
    locs(2,i)=find(t<window_size*i,1,'last');

    if(locs(2,i)<=locs(1,i)),continue; end

    % note: if 'economy' is true, this may return invalid matrix sizes
    coeff2(:,:,i) = pca(X(locs(1,i):locs(2,i),:),'Economy',false);
    means2(:,i)=mean(X(locs(1,i):locs(2,i),:));
end

% first dimension of coeff are the eigenvectors
coeff=zeros(N,N,M);
t_window2 = [t(1) t_window t(end)];
for i = 1:N
    % interpolate points at time
    c = interp1(t_window2, [coeff2(:,i,1)'; squeeze(coeff2(:,i,:))'; coeff2(:,i,end)'], t, method);
    
    % make unit length
    coeff(:,i,:) = (c./repmat(sqrt(sum(c.^2,2)),1,N))';
end

% interpolate mean value as well
means2 = [means2(:,1) means2 means2(:,end)];
means = interp1(t_window2,means2',t,method);
X = X-means;
% project data onto coefficients.
for i = 1:M
    score(i,:) = X(i,:)*squeeze(coeff(:,:,i));
end