function data = filter_voltage(data, sampling_rate, bandpass, order)
%FILTER_VOLTAGE band- or high-pass filter voltage signal with filter
% filter_voltage.m - Filters a set of data channels using default filter
% settings (600-9000Hz bandpass). Data is converted to single precision.
%
% Usage: data = filter_voltage(data, sampling_rate, bandpass);
%
% Input:  data                    * CxD matrix of raw voltage data
%         sampling_rate           * sampling rate in Hz
%         bandpass (optional)     * bandpass of filter
%         filter_order (optional) * order of Butterworth filter
% Output: data                    * CxD matrix of filtered data
%
% Written by Marshall Crumiller
% email: mcrumiller@gmail.com
%--------------------------------------------------------------------------
if(~exist('bandpass','var') || isempty(bandpass)), bandpass = [600 9000]; end
if(~exist('order','var') || isempty(order)), order=2; end
if(length(bandpass)==1)
    [b,a] = butter(order,bandpass/(sampling_rate/2),'high');
elseif(length(bandpass)==2)
    [b,a] = butter(order,[bandpass(1) bandpass(2)]/(sampling_rate/2),'bandpass');
else
    error('Invalid filter pass.  Should be 1x1 (highpass) or 1x2 (bandpass).');
end

data=single(data);
for i = 1:size(data,1)
    data(i,:) = single(filtfilt(b,a,double(data(i,:))));
end