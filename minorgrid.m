function h_minorgrid = minorgrid(xx, yy, ax, top, lims)
%--------------------------------------------------------------------------
% minorgrid.m - Sets minor grid ticks.  Default is gray [.5 .5 .5].
%
% Usage: h_minorgrid = minorgrid(xx, yy, ax)
%
% Input:  xx                  * X grid locations (leave empty for blank)
%         yy (optional)       * Y grid locations (leave empty for blank)
%         ax (optional)       * axis values
%         top (optional)      * plot grid lines on top of everything else
%                               in plot.  Default is underneath.
%         lims (optional)     * plot with the provided limits
% Output:
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
if(~exist('ax','var'))
    ax=gca;
end

if(exist('lims','var') && ~isempty(lims) && length(lims)==4)
    a=lims;
else
    a=axis(ax);
end
hold(ax,'on');
gray=[.4 .4 .4];

children = get(ax,'children');

% automatically figure out ticks
if((~exist('xx','var') && ~exist('yy','var')) || (isempty(xx) && isempty(yy)))
    num_ticks = 4;
    dx=(a(2)-a(1))/num_ticks; dy=(a(4)-a(3))/num_ticks;
    xx=a(1)+dx*(1:num_ticks-1); yy=a(3) + dy*(1:num_ticks-1);
end

if(~isempty(xx))
    xx=xx(:)'; numX=length(xx);
    x_locs = reshape([xx; xx; NaN(1,numX)],1,[]);
    y_locs = reshape([ones(1,numX)*a(3); ones(1,numX)*a(4); NaN(1,numX)],1,[]);
end
if(exist('yy','var') && ~isempty(yy))
    yy=yy(:)'; numY=length(yy);
    x_locs = [x_locs reshape([ones(1,numY)*a(1); ones(1,numY)*a(2); NaN(1,numY)],1,[])];
    y_locs = [y_locs reshape([yy; yy; NaN(1,numY)],1,[])];
end
%{
if(exist('lims','var') && length(lims)==4)
    x_locs(x_locs==min(x_locs))=lims(1);
    x_locs(x_locs==max(x_locs))=lims(2);
    y_locs(y_locs==min(y_locs))=lims(3);
    y_locs(y_locs==max(y_locs))=lims(4);
end
%}
h_minorgrid = plot(ax,x_locs,y_locs,':','linewidth',1,'color',gray);

% put gridlines into the back
if(~exist('top','var') || top==false)
    set(ax,'Children',[children(:)' h_minorgrid]);
end