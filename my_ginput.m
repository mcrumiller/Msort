function [x, y] = my_ginput(ax, color,return_to_state,display_handle)
%MY_GINPUT Allow user to select threshold with the mouse
%   [x y] = MY_GINPUT allows the user to select a threshold with the mouse.
%   A dotted line tracks the vertical position of the mouse, and a small
%   crosshair indicates the [x,y] position of the mouse.  The results (x,y)
%   are returned as coordinates located within the current axis.
%
%   [x y] = my_ginput(ax, color) allows the user to specify the color of
%   the dotted line.
%
%   Written by Marshall Crumiller
%
%   Updates:
%      2012.04.09 - Improved execution speed, removed unnecessary variables
if(~exist('ax','var') || isempty(ax)), ax=gca; end
if(~exist('color','var') || isempty(color)), color=[.7 0 0]; end
if(~exist('return_to_state','var') || isempty(return_to_state))
    return_to_state=true;
end
if(exist('display_handle','var') && ishandle(display_handle))
    DISPLAY_VAL=true;
else
    DISPLAY_VAL=false;
end

fig=ancestor(ax,'figure');
figure(fig);

ax_lim=axis(ax);
height=ax_lim(4)-ax_lim(3);
Y_height=height*.03;

% get fig position
tmp=get(fig,'units'); set(fig,'units','pixels');
fig_pos = get(fig,'position'); set(fig,'units',tmp);

% get axis position
tmp = get(ax,'units'); set(ax,'units','pixels');
ax_pos = get(ax,'position'); set(ax,'units',tmp);

set(0,'units','pixels');

% determine axis coordinates
g=ax;
x_offset=fig_pos(1); y_offset=fig_pos(2);
while(g~=fig)
    tmp=get(g,'units'); set(g,'units','pixels');
    pos=get(g,'position');
    x_offset=x_offset+pos(1);
    y_offset=y_offset+pos(2);
    set(g,'units',tmp);
    g=get(g,'parent');
end

ax_coord = [x_offset y_offset x_offset y_offset] + [0 0 ax_pos(3) ax_pos(4)];

% determine conversion
x_ratio=(ax_lim(2)-ax_lim(1))/ax_pos(3);
y_ratio=(ax_lim(4)-ax_lim(3))/ax_pos(4);

% set up initial plot
mouse_pt = get(0,'PointerLocation');
x=(mouse_pt(1)-ax_coord(1))*x_ratio+ax_lim(1);
y=(mouse_pt(2)-ax_coord(2))*y_ratio+ax_lim(3);
locsX=[ax_lim(1) ax_lim(2) NaN x x];
locsY=[y y NaN y-Y_height y+Y_height];
hold on; axis manual;
h=plot(locsX,locsY,'--','linewidth',1.5,'Color',color);
h2=plot(locsX,-locsY,'--','linewidth',1,'Color',color,'visible','off');

inside=false;
double=false;

set(fig,'WindowButtonMotionFcn',@c_cursor);
set(fig,'WindowButtonDownFcn',@get_pt);
set(fig,'WindowKeyPressFcn',@doubleswap_down);
set(fig,'WindowKeyReleaseFcn',@doubleswap_up);

% when user moves mouse
    function c_cursor(hObject,eventData)
        
        mouse_pt = get(0,'PointerLocation');
        
        % if we're inside axis
        if(mouse_pt(1)>=ax_coord(1) && mouse_pt(1)<=ax_coord(3) && mouse_pt(2)>=ax_coord(2) && mouse_pt(2)<=ax_coord(4))
            mouse_pt = get(0,'PointerLocation');
            
            % hide mouse cursor if we have to
            if(~inside)
                inside=true;
                set(fig, 'PointerShapeCData', nan(16),'Pointer','custom');
            end
            
            % calculate new point
            x=(mouse_pt(1)-ax_coord(1))*x_ratio+ax_lim(1);
            y=(mouse_pt(2)-ax_coord(2))*y_ratio+ax_lim(3);
            
            locsX=[ax_lim(1) ax_lim(2) NaN x x];
            locsY=[y y NaN y-Y_height y+Y_height];
            set(h,'XData',locsX,'YData',locsY);
            set(h2,'XData',locsX,'YData',-locsY);
            if(DISPLAY_VAL)
                set(display_handle,'String',sprintf('%2.1f std devs',y));
            end
        else
            inside=false;
            set(fig,'pointer','arrow');
        end
    end

% When user clicks mouse, return [x,y] location
    function get_pt(~,~)
        mouse_pt = get(0,'PointerLocation');
        if(mouse_pt(1)>=ax_coord(1) && mouse_pt(1)<=ax_coord(3) && mouse_pt(2)>=ax_coord(2) && mouse_pt(2)<=ax_coord(4))
            x=(mouse_pt(1)-ax_coord(1))*x_ratio+ax_lim(1);
            y=(mouse_pt(2)-ax_coord(2))*y_ratio+ax_lim(3);
            set(fig,'WindowButtonMotionFcn','','WindowButtonDownFcn','','WindowKeyPressFcn','','WindowKeyReleaseFcn','');
            delete([h h2]);
            
            if(double),y=[abs(y) -abs(y)]; end
        end
    end

    % user presses shift key
    function doubleswap_down(hObject,eventData)
        if(strcmp(eventData.Key,'shift'))
            set(h2,'visible','on');
            double=true;
        end
    end

    % user releases shift key
    function doubleswap_up(hObject,eventData)
        if(strcmp(eventData.Key,'shift'))
            set(h2,'visible','off');
            double=false;
        end
    end
waitfor(fig,'WindowButtonDownFcn','');
if(return_to_state), set(fig,'pointer','arrow'); end
end