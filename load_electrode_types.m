function [names,configurations,geometry]=load_electrode_types(electrode_file)
% LOAD_ELECTRODE_TYPES read data from electrode configuration file
%
%     [names,configuration,geometry = LOAD_ELECTRODE_TYPES(electrode_file)
%     read data from the provided electrode_file, which is a file
%     describing the channel groups and optionally the channel physical
%     locations.
%
% Input:  electrode_file         * path to electrode configuration file
% Output: names                  * 1xE cell array of electrode names
%         configurations         * 1xE cell array of electrode channel
%                                  groupings
%         geometry               * 1xE cell array of electrode geometry
%                                  information
%
% Written by Marshall Crumiller
% email: mcrumiller@gmail.com
%--------------------------------------------------------------------------
f=fopen(electrode_file);

tline=fgetl(f);

names=cell(1);
configurations=cell(1);
geometry = cell(1);
ind=1;
while(ischar(tline))
    % pick up electrode name and configuration
    result=regexp(tline,'(\s)*(?<comment>%)?(\s)*(?<name>.*[^ ])(?= ?= ?)( ?= ?)(?<configuration>{.*});?','names');
    if(~isempty(result) && isempty(result.comment))
        % load electrode geometry
        if(regexp(result.name,'.*geometry.*'))
            eval(sprintf('geometry{%g}=%s;',ind-1,result.configuration));
        % load electrode name and configuration
        else
            names{ind}=result.name;
            eval(sprintf('configurations{%g}=%s;',ind,result.configuration));
            ind=ind+1;
        end
    end
    tline=fgetl(f);
end
fclose(f);

N=length(names); G=length(geometry);
if(G<N)
    geometry(end+1:N)=cell(1,N-G);
end