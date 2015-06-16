function ind = find_index(X, vals)
%--------------------------------------------------------------------------
% find_index.m - Given a vector X and a set of values v, returns the
% indices of the values in X that most closely matches the values in vals.
%
% Usage: ind = find_index(X, val);
%
% Input:  X                   * vector of values
%         val                 * vector value to match
% Output: ind                 * index location of (closest) val in X
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------

% note: this shit is sexy
%{
lenX=length(X); lenvals = length(vals);
vals=repmat(vals',1,lenX);
X=repmat(X,lenvals,1);
d=abs(X-vals);
m=min(d,[],2);
m=repmat(m,1,lenX);
z=~(d-m);
[I,J] = find(z);
[~,m,~] = unique(I,'first');
ind=J(m);
%}
ind=zeros(1,length(vals));
for i = 1:length(vals)
    d=abs(X-vals(i));
    ind(i)=find(d==min(d),1,'first');
end