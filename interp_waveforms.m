function [waveforms,t] = interp_waveforms(waveforms, t, interp_factor)
%--------------------------------------------------------------------------
% interp_waveforms.m - Interpolates waveforms by a factor of interp_factor.
%
% Usage: [waveforms,t] = interp_waveforms(waveforms, t, interp_factor);
%
% Input:  waveforms          * CxWxT matrix of W waveforms on C channels
%         t                  * 1xT time vector
%         interp_factor      * interpolation factor
% Output: waveforms          * CxWx(T*interp_factor) matrix of waveforms
%         t                  * new time vector
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
t2=linspace(t(1),t(end),(length(t)-1)*interp_factor+1);

siz=size(waveforms);
if(siz(1)==1 && siz(2)==1),
    w=zeros(1,1,length(t2),'single');
    w(1,1,:)=interp1(t,waveforms(:),t2,'spline');
    waveforms=w;
else
    % do some tests for single dimensions
    w=permute(single(waveforms),[3 1 2]);
    waveforms=interp1(t,w,t2,'spline');
    waveforms=permute(waveforms,[2 3 1]);
end
t=t2;