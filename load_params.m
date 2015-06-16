function params = load_params(filename)
% LOAD_PARAMS read data from configuration file
%
%     params = LOAD_PARAMS(filename) reads the data from a filename and
%     returns the data in a struct. The params file contains information
%     in the following format:
%
%        % comments like this are ignored
%        name1 = value1;
%        name2 = value2;
%
% Input:  filename               * path to configuration file
% Output: params                 * struct containing parameter information.
%                                  using the above example, the data is
%                                  stored as:
%
%        params.name1 = value
%        params.name2 = value
%
% Written by Marshall Crumiller
% email: mcrumiller@gmail.com
%--------------------------------------------------------------------------
if(~exist('filename','var')), filename='params.txt'; end

fid = fopen(filename);
params = [];

while(true)
    C = textscan(fid, '%s = %s\n');
    if(isempty(C{1}))
        break;
    elseif(isempty(C{2}))
        continue;
    else
        key=C{1}{1}; val=C{2}{1};
        if(isnan(str2double(val)))
            if(val(end)==';'),val=val(1:end-1); end
            eval(sprintf('params.%s = ''%s'';',key,val));
        else
            eval(sprintf('params.%s = %g;',key,val));
        end
    end
end

fclose(fid);