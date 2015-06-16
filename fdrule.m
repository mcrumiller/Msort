function [bin_width,edges] = fdrule(x)
%FDRULE   Freedman–Diaconis rule for histogram bin size
%   bin_width = FDRULE(x) returns the optimal bin width for a distribution using the Freedman–Diaconis rule, which is
%   defined as:
%       bin_width = 2*IQR(x) * n^(-1/3)
%   where IQR is the inter-quartile range, x is the vector o fdata, and n is the number of observations.
%   
%   Input:  x........................1xN vector of data
%   Output: bin_width................optimal bin width
%           edges....................1xE edge vector for use with HISTC
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-05-12: Created
%-----------------------------------------------------------------------------------------------------------------------
bin_width=2*iqr(x)*length(x)^(-1/3);
m1=min(x);
N=ceil((max(x)-m1)/bin_width);
edges=(0:N)*bin_width+m1;