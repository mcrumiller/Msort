function [count, u] = count_uniques(x)
%COUNT_UNIQUES number of each unique value in an array
%   [count,u]=count_uniques(x)  returns the number of unique values of each
%   item in the array x, along with their unique values.
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
x=sort(x);
[u,I]=unique(x);
count=diff([I;length(x)+1]);