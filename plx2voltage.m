function conversion_factor = plx2voltage(plx_file, channels)
%--------------------------------------------------------------------------
% plx2voltage.m - Reads .plx header files to determine what transformation
% is required to convert channels from the uint8 format into mV.
%
% % Note: the formula for converting from the uint8 data format to mV is:
%                                 5V
%     V    =   x * ------------------------------------------
%                   (2048)(headstage gain)(amp gain)(NI*gain)
%      where
%      headstage = 20
%      amp = 1000
%      NI = [looked up in plexon file]
% Usage: conversion_factor = plx2voltage(plx_file, channels);
%
% Input:  channels            * 1xC vector of channel numbers
% Output: conversion_vfactor  * 1xC vector of conversion factors
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
num_channels=length(channels);

% amplification factors used
headstage = 20;
%amp = 1000;
% note: amp is included in NI


[~,NI]=plx_adchan_gains(plx_file);
NI=NI(channels)';

conversion_factor=ones(1,num_channels)*5./(2048*headstage*NI);


% The best accurate way is to let Plexon do the conversion internally
%ad_reg=zeros(num_channels,2);
%ad_v=zeros(num_channels,2);
%for i = 1:num_channels
%    [~,~,ad_reg(i,:)]=plx_ad_span(plx_file,channels(i),1,2);
%    [~,~,ad_v(i,:)]=plx_ad_span_v(plx_file,channels(i),1,2);
%end
%conversion_factor=(ad_v(:,1)./ad_reg(:,1))';