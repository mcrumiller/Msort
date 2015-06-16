function fix_feature_menu(handles)
%FIX_FEATURE_MENU   Fix feature menu for plotting
%   FIX_FEATURE_MENU updates the contents of the menu containing features for extracting from waveforms.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
% Determine which features we're going to calculate
num_channels=getappdata(handles.output,'num_channels');
channels=getappdata(handles.output,'channels');
num_PCs = get(handles.num_PCs_slider,'Max');
if(~iscell(num_PCs)), num_PCs={num_PCs}; end
num_PCs=max([num_PCs{:}]);
feature_list = cell(1,num_channels*num_PCs);
feature_names = cell(1,num_channels*num_PCs);

c=round(getappdata(handles.output,'channel_colors')*255);
ind=1;
for i = 1:num_channels
    for j = 1:num_PCs
        feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - PC %g</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i),j);
        feature_names{ind}=sprintf('%g - PC %g',channels(i),j);
        ind=ind+1;
    end
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Peak</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Peak',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Vall</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Valley',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Slope</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Slope',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Ener</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Energy',channels(i)); ind=ind+1;
    feature_list{ind}=sprintf('<html><font color="rgb(%g,%g,%g)"><b>%g - Amp</b></font></html>',c(i,1),c(i,2),c(i,3),channels(i));
    feature_names{ind}=sprintf('%g - Amplitude',channels(i)); ind=ind+1;
end

% note: 'time' means spike time, 'hist' creates a 3D of the other features
% displayed
feature_list = [feature_list 'Time'];
feature_names = [feature_names 'Time'];

set(handles.feature1_menu,'String',feature_list);
set(handles.feature2_menu,'String',[feature_list(:);{'Hist'}]);
set(handles.feature3_menu,'String',[{'<none>'};feature_list(:);{'Hist'}]);
set(handles.feature1_menu,'Value',1);
if(num_channels>1),val=8; else val=2; end
set(handles.feature2_menu,'Value',val);
set(handles.feature3_menu,'Value',1);
setappdata(handles.output,'feature_names',feature_names);