function generate_standard_ABF_xml(abf_file)
%--------------------------------------------------------------------------
% generate_standard_ABF_xml.m - 
%
% Usage: generate_standard_ABF_xml();
%
% Input:  
% Output: 
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
[data,~,h]=abfload(abf_file);

Fs = h.lADCResolution; dt=1/Fs;

% Process filename
expr = regexp(abf_file,'(?<dirname>(.*[\\/])*)(?<basename>.*)\.abf','names');
if(isempty(expr)), error('Invalid filename.'); end
dirname = expr.dirname;
if(isempty(dirname)), dirname = [pwd '/']; else dirname=expr.dirname; end;
basename = expr.basename;

plx_file=sprintf('%s.abf',basename);
mat_file=sprintf('%s.mat',basename);
xml_file=sprintf('%s%s.xml',dirname,basename);

%s=struct('experiment',struct('group',{}));

s.experiment.plx_file.Text=plx_file;
s.experiment.mat_file.Text=mat_file;
s.experiment.electrode.Text='Single-Channel';
s.experiment.ADRate=sprintf('%g',Fs);
s.experiment.DateTime.Text='0/ 0/0000 00:00:00';
s.experiment.duration.Text=sprintf('%g',length(data)*dt);
s.experiment.shank{1}.channel{1}.threshold.Text='0.05';
s.experiment.shank{1}.Attributes.num='1';
s.experiment.shank{1}.channel{1}.Attributes.num='1';
s.experiment.shank{1}.channel{1}.enabled.Text='true';


% put field names in order
fieldnames = {'plx_file','mat_file','electrode','ADRate','DateTime','duration','shank'};
s.experiment=orderfields(s.experiment,fieldnames);

struct2xml(s,xml_file);