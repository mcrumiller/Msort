function idx = KlustaKwik_java(features, args, KlustaKwik_jarfile)
%KLUSTAKWIK   KlustaKwik wrapper
%   idx = KLUSTAKWIK(features) runs the KlustaKwik utility to sort a matrix
%   of features and returns clustering results. This wrapper performs the
%   following steps:
%
%      1. Writes features to a temporary file (file.fet.1)
%      2. Runs the java system command.
%      3. Reads the temporarly cluster file (file.clus.1)
%      4. deletes temporary files
%
%   Running this script requires KlustaKwik.jar.  If the jar file is not
%   located in the Matlab Path, you can provide the directory (or full path
%   filename) of the jarfile containing the KlustaKwik code.
%
%   Usage: idx = KlustaKwik(features)
%
%   Input:  features             * WxD feature matrix; W observations, D
%                                  variables (dimensions)
%           args                 * (optional) argument list for KlustaKwik
%           KlustaKwik_jarfile   * (optional) location of KlustaKwik.jar
%   Output: idx                  * Wx1 index matrix
%
%   Example: idx = KlustaKwik(x,'-MinClusters 2');
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2014-05-23: Created
%--------------------------------------------------------------------------
% determine whether java is installed


% determine existence of KlustaKwik.jar
jarfile = 'KlustaKwik.jar';
if(exist('KlustaKwik_jarfile','var')  && ~isempty(KlustaKwik_jarfile))
    expr = regexp(KlustaKwik_jarfile,'(?<basedir>.*[\\/])*(?<jarfile>.*\.jar)*','names');
    if(isempty(expr)), error('Invalid filename.'); end
    basedir = expr.basedir;
    if(~isempty(expr.jarfile)), jarfile = expr.jarfile; end
    jarfile = sprintf('%s%s',basedir,jarfile);
end
if(~exist(jarfile,'file'))
    error('Cannot locate KlustaKwik.jar.  Try providing it as a second argument.');
end
jarfile=which(jarfile);


% Step 1: Make Feature file
% create temp name. Note: must be in current directory
base_name=tempname(pwd);
feature_file = sprintf('%s.fet.1', base_name);
f = fopen(feature_file,'w');
fprintf(f,'%g\n',size(features,2)); fclose(f);
%precision=ceil(log10(max(max(max(double(features))))));
precision = 5;
dlmwrite(feature_file,features,'-append','newline','unix','precision',precision,'delimiter','\t');

% Step 2: run KlustaKwik

% Determine KlustaKwik command
% grab argument list if user provided
arglist = '-MinClusters 2';
if(exist('args','var') && ~isempty(args))
    arglist = sprintf('%s ',args);
end

% Execute command
%LOGFILE=false; % enable this if you want to generate a log file
%if(LOGFILE), logfile = sprintf('Klustakwik-%s.txt',strrep(strrep(datestr(now),' ','-'),':','.'));
%else logfile = 'NUL'; end
logfile='NUL';
cmd = sprintf('java -jar %s -ElecNo 1 -FileBase "%s" -UseFeatures ALL %s > %s',jarfile, base_name, arglist, logfile);
%cmd = sprintf('KlustaKwik "%s" 1 -MinClusters 2 -MaxClusters 6 > %s', base_name, logfile);
system(cmd);

% Step 3: read Clusters file
cluster_file = sprintf('%s.clu.1',base_name);
f_clus = fopen(cluster_file);
if(f_clus==-1), error('Cannot open cluster file.'); end
idx=fscanf(f_clus,'%g')';
idx=idx(2:end); % ignore first line
fclose(f_clus);

% Step 4: delete clusters file and feature file
delete(cluster_file);
delete(feature_file);