function [peaks,peaklocs] = find_peaks(v, threshold)
%--------------------------------------------------------------------------
% find_peaks.m - Locates the peaks of a set of functions
%
% Usage: [peaktimes, trace_interp] = find_peaks(v, SCALING_FACTOR, ...
%                                               show_plot, threshold)
% Input:  v(*)                           * MxN matrix of N functions
%         threshold                      * proportion of max peak for each
%                                          histogram; ignore values below
%         show_plot (optional)           * display plot of v with peaks
% Output: peaktimes                      * approximate peak times
%         trace_interp                   * interpolated original vector
%
% (*) If a matrix is given, the peaks are returned as a zero-padded matrix,
% with each row containing the location of the peaks for the corresponding
% input function.
%
% (**) If SCALING_FACTOR > 1, the original vector will be interpolated and
% the peaktimes returned will correspond to those in the interpolated,
% finer vector (returned as trace_interp).
%
% Written by Marshall Crumiller
% Last Updated September 12th, 2010
%    -simplified program, utilized findpeaks.m, vast speedup
%--------------------------------------------------------------------------
if(isvector(v)), v=v(:)'; end
M=size(v,1);

% duplicate threshold if necessary
if(~exist('threshold','var')),threshold=0; end
if(length(threshold)~=M), threshold=repmat(threshold,1,M); end

peaks_cell=cell(M,1);
peaklocs_cell=cell(M,1);
% find peaks for each row
for i = 1:M
    [peaks_cell{i},peaklocs_cell{i}]=findpeaks(double(v(i,:)),'MINPEAKHEIGHT',threshold(i));
end

max_size=max(cellfun(@length,peaks_cell));

peaks=cell2mat(cellfun(@(x) horzcat(x(:)',zeros(1,max_size-length(x))),peaks_cell,'uniformoutput',false));
peaklocs=cell2mat(cellfun(@(x) horzcat(x(:)',zeros(1,max_size-length(x))),peaklocs_cell,'uniformoutput',false));
