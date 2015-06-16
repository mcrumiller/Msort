function set_theme(handles)
%SET_THEME   Set color theme of Msort
%   SET_THEME generates a color palette based on an available option of themes. Modify this function in order to add
%   additional color themes.
%
%   Note: theme must implement the following values:
%     palette.background
%     palette.gridlines
%     palette.text
%     palette.header
%     palette.button
%     palette.button_text
%     palette.undo
%     palette.green
%     palette.red
%     palette.orange
%     palette.blue
%     palette.ax_background
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global palette;
theme = 'material';

switch(theme)
    case 'material'
        colors=import_colors(theme);
        palette.background=colors.light_grey;
        palette.gridlines=colors.grey;
        palette.text=colors.dark_grey;
        palette.header=colors.blue;
        palette.button=1-(1-colors.blue)/3;
        palette.button_text=colors.dark_grey;
        palette.undo=1-(1-colors.pink)/2;
        palette.view_button=1-(1-colors.indigo)/3;
        palette.green=colors.green;
        palette.red=colors.red;
        palette.orange=colors.orange;
        palette.blue=colors.indigo;
        palette.teal=colors.teal;
        palette.ax_background = 1-(1-colors.light_grey)/2;
        cluster_colors=[colors.red;colors.purple;colors.indigo;colors.teal;...
                        colors.brown;colors.blue_grey;colors.green;colors.orange;...
                        colors.pink;colors.green;colors.indigo;colors.amber;...
                        colors.cyan;colors.deep_purple;colors.light_green;...
                        colors.dark_grey;colors.pink];
    case 'dark'
        palette.background=[0 0 0];
        palette.gridlines=[.5 .5 .5];
        palette.text=[.7 .7 .7];
        palette.header=[.87 .87 .87];
        palette.button=[237 237 237]/255;
        palette.button_text=[.2 .2 .2];
        palette.undo=[236 214 214]/255;
        palette.view_button=[186 212 244]/255;
        palette.green=[0 .7 0];
        palette.red=[.7 0 0];
        palette.orange=[255 127 39]/255;
        palette.blue=[0 0 .7];
        palette.ax_background = [0 0 0];
        cluster_colors=distinguishable_colors(10,'k');
    case 'solarized'
        
end
% set background
set([handles.filename_panel handles.information_panel handles.features_panel handles.feature_control_panel handles.options_panel handles.cluster_panel ...
    handles.group_panel handles.group_cmd_panel handles.params_panel handles.clusters_panel handles.cleanup_panel...
    handles.interface_panel handles.channel_panel handles.status_panel],...
    'backgroundcolor',palette.background,'foregroundcolor',palette.header);
set([get(handles.group_panel,'children')' get(handles.group_cmd_panel,'children')' ...
    get(handles.clusters_panel,'children')' get(handles.cleanup_panel,'children')' ...
    get(handles.interface_panel,'children')' get(handles.information_panel,'children')'...
    get(handles.channel_panel,'children')' get(handles.options_panel,'children')' ...
    get(handles.filename_panel,'children')' get(handles.status_panel,'children')' ...
    get(handles.params_panel,'children')' get(handles.feature_control_panel,'children')'],...
    'backgroundcolor',palette.background,'foregroundcolor',palette.text);

% set buttons
set(findobj(handles.output,'style','pushbutton'),'backgroundcolor',palette.button,'foregroundcolor',palette.button_text);
set([handles.deselect_button handles.undo_button],'backgroundcolor',palette.undo,'foregroundcolor',palette.button_text);

set(handles.output,'color',palette.background);

% determine palette-based cluster colors
palette.cluster_colors=cluster_colors;