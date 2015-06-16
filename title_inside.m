function t=title_inside(title_text,varargin)
%--------------------------------------------------------------------------
% title_inside.m - Places the title inside the plot so as to conserve
% space.
%
% Usage: title_inside(title_text,fontsize);
%
% Input:  title_text            * text to display in title
%         fontsize (optional)   * size of font
% Output: t                     * text handle
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
args=false;
if(exist('varargin','var') && ~isempty(varargin))
    args=true;
    arglist=[];%varargin{1};
    for i = 1:length(varargin)
        c=varargin{i};
        if(isfloat(c))
            if(isvector(c))
                temp_s=sprintf('[%g',c(1));
                for j=2:length(c)
                    temp_s=sprintf('%s %g',temp_s,c(j));
                end
                arglist=sprintf('%s,%s]',arglist,temp_s);
            else
                arglist=sprintf('%s,[%g]',arglist,c);
            end
        else
            arglist=sprintf('%s,''%s''',arglist,c);
        end
    end
end

a = axis(gca);
x=(a(1)+a(2))/2;
height=a(4)-a(3);
margin=height*.018;
y=a(4)-margin;

% Generate text
if(args)
    eval(sprintf('t=text(x,y,title_text,''horizontalalignment'',''center'',''verticalalignment'',''top''%s);',arglist));
else
    t=text(x,y,title_text,'horizontalalignment','center','verticalalignment','top');
end