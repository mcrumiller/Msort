function abf2(abf_file)
%--------------------------------------------------------------------------
% abf2mat2.m - Gets as input a abf file containing continuous channel
% traces. An associated .xml file is generated if it doesn't exist, along
% with a .mat file that contains all of the filtered and unfiltered spike
% data.
%
% If provided, the user can specify which specific shanks should be
% converted to .mat; this saves time if one only wants to convert a subset
% of the channels to .mat.
%
% Usage: abf2mat2(abf_file,selected_shanks);
%
% Input:  abf_file                  * .abf file containing continuous
%                                      recording
%         selected_shanks           * vector of shank numbers
% Output:
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
% get subdirectory and base name of file
expr = regexp(abf_file,'(?<dirname>(.*[\\/])*)(?<basename>.*)(?<ext>\.abf|\.xml)+','names');
if(isempty(expr)), error('Invalid filename.'); end
dirname = expr.dirname;
if(isempty(dirname))
    dirname = pwd;
    if(dirname(end)~='/' || dirname(end)~='\')
        dirname=sprintf('%s\\',dirname);
    end
else dirname=expr.dirname;
end;
basename = expr.basename;

% get file names
abf_file = sprintf('%s%s.abf',dirname,basename);
xml_file = sprintf('%s%s.xml',dirname,basename);
mat_file = sprintf('%s%s.mat',dirname,basename);

if(~exist(mat_file,'file'))
    % create empty mat file
    save(mat_file,'-v7.3','');
end

% Load data from parameter file
if(~exist(abf_file,'file')), error('.abf file does not exist.'); end
if(~exist(xml_file,'file'))
    fprintf('.xml file does not exist.  Generating standard xml...\n');
    generate_standard_ABF_xml(abf_file);
end

% grab shank and channel information from the .abf file
s = xml2struct(xml_file);

num_shanks = length(s.experiment.shank);

% Determine which channels have already been processed in the .mat file
file_vars=whos('-file',mat_file);
vars={file_vars.name}; expr=regexp(vars,'data_(?<channel>\d*$)','names');

[data,~,h]=abfload(abf_file);
Fs = h.lADCResolution; dt=1/Fs;
Ts = (0:length(data)-1)*dt;
data_1 = filter_voltage(data',Fs,[300 1000]);
std_1=std(data);
data_1=data_1./std(data);

save(mat_file,'data_1','Fs','Ts','data_1','std_1');