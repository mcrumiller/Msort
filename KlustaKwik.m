function idx = KlustaKwik(features, args)
%--------------------------------------------------------------------------
% KlustaKwik.m - KlustaKwik wrapper for Matlab.  Given a feature matrix,
% runs the KlustaKwik utility to sort, and returns clustering results.
% This wrapper performs the following steps:
%
%   1. Writes features to a temporary file (file.fet.1)
%   2. Runs the KlustaKwik system command.
%   3. Reads the temporarly cluster file (file.clus.1)
%   4. deletes temporary files
%
% Usage: idx = KlustaKwik(features)
%
% Input:  features             * WxD feature matrix; W observations, D
%                                variables (dimensions)
%         args                 * (optional) argument list for KlustaKwik
% Output: idx                  * Wx1 index matrix
%
% Example: idx = KlustaKwik(x,'-MinClusters 2');
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
% Step 1: Make Feature file
% create temp name. Note: must be in current directory
base_name=tempname(pwd);
output_file = sprintf('%s.fet.1', base_name);
f = fopen(output_file,'w');
fprintf(f,'%g\n',size(features,2)); fclose(f);
prec1=min(double(features(:))); prec2=max(double(features(:)));
%precision=ceil(log10(max(ceil([1/prec1 prec2]))));
%precision=max(ceil(log10(max(features(:)))),5);
precision = 4;
dlmwrite(output_file,abs(features),'-append','newline','unix','precision',precision,'delimiter','\t');

% Step 2: run KlustaKwik

% Determine KlustaKwik command
% grab argument list if user provided
%arglist = '-MinClusters 2';
%if(exist('args','var') && ~isempty(args))
%    arglist = sprintf('%s ',args);
%else
%    args=[];
%end
if(~exist('args','var')),args=[]; end

% Execute command
LOGFILE=false; % enable this if you want to generate a log file
if(LOGFILE)
    logfile = sprintf('Klustakwik-%s.txt',strrep(strrep(datestr(now),' ','-'),':','.')); %#ok<UNRCH>
    cmd = sprintf('KlustaKwik "%s" 1 %s -Screen 1 > %s', base_name, args,logfile);
    out=system(cmd);
else
    cmd = sprintf('KlustaKwik "%s" 1 %s -Log 0 -Screen 1 > KlustaKwik.txt', base_name, args);
    out=system(cmd);
end
%cmd = sprintf('KlustaKwik "%s" 1 %s > %s', base_name, args,logfile);
%out=system(cmd);

% Step 3: read Clusters file
cluster_file = sprintf('%s.clu.1',base_name);
f_clus = fopen(cluster_file);
if(f_clus==-1), error('Cannot open cluster file.'); end
idx=fscanf(f_clus,'%g')';
idx=idx(2:end); % ignore first line
idx=idx-min(idx)+1;
fclose(f_clus);

% Step 4: delete cluster and feature file
delete(cluster_file);
delete(output_file);