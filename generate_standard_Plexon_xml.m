function generate_standard_Plexon_xml(plx_file,electrode)
%--------------------------------------------------------------------------
% generate_standard_Plexon_xml.m - Given a .plx experiment recording file,
% generates a standard template for an associated .xml file.  This template
% knows the Plexon channel mapping per shank (given the 8-shank/64-channel
% probes), and gives each channel a default threshold of 5 standard
% deviations.
%
% Usage: generate_standard_Plexon_xml(plx_file, old);
%
% Input:  plx_file                * .plx file
%         num_shanks              * # shanks used (usually 4 or 8)
%         old (optional)          * flag to swap ch1-32 and ch33-64
% Output: <xml file>
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
global electrode_type;
if(~exist('plx_file','var') || ~exist(plx_file,'file'))
    error('PLX file not found.');
end

% get subdirectory and base filename
expr = regexp(plx_file,'(?<dirname>(.*[\\/])*)(?<basename>.*)\.plx','names');
if(isempty(expr)), error('Invalid filename.'); end
dirname = expr.dirname;
if(isempty(dirname)), dirname = [pwd '/']; else dirname=expr.dirname; end;
basename = expr.basename;

% replace backslashes with forward slashes to make this UNIX-compliant
dirname=strrep(dirname,'\','/');
basename=strrep(basename,'\','/');
plx_file=strrep(plx_file,'\','/');

% Check for associated .xml file
xml_file=sprintf('%s%s.xml',dirname,basename);
if(exist(xml_file,'file'))
    reply = input(sprintf('%s already found.  Overwrite? (y/n) ',xml_file),'s');
    if(strcmpi(reply,'y')~=true), return; end
end

% determine .mat file name
mat_file=sprintf('%s%s.mat',dirname,basename); %#ok<NASGU>

s=struct('experiment',struct('group',{}));
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, duration, DateTime] = plx_information(plx_file);

% Parse date and time
expr = regexp(DateTime,'\D*(?<month>\d*)/\D*(?<day>\d*)/\D*(?<year>\d*)\D*','names');
if(str2double(expr.year)<2012), old=true; end

% determine which shanks are used
fullread=false;
[~, ~, ~, contcounts] = plx_info(plx_file, fullread);
channels=find(contcounts); %#ok<EFIND>
full_channels=1:(ceil(find(contcounts,1,'last')/32)*32);
if(isempty(channels))
    fprintf('Error: channels is empty.  Doing full read on the file. This may take a while...');
    fullread = true;
    [~, ~, ~, contcounts] = plx_info(plx_file, fullread);
    channels=find(contcounts); %#ok<EFIND>
    full_channels=1:(ceil(find(contcounts,1,'last')/32)*32);
    if(isempty(channels))
        error('PLX file has no continuous data.');
    end
    fprintf('done.\n');
end


% generate structure
s.experiment=struct('plx_file',{},'mat_file',{},'ADRate',{},'DateTime',{},'duration',{},'shank',{});
s.experiment(1).plx_file.Text=sprintf('%s.plx',basename);
s.experiment(1).mat_file.Text=sprintf('%s.mat',basename);
[~,freqs] = plx_adchan_freqs(plx_file);
s.experiment(1).ADRate.Text=sprintf('%g',freqs(1));
s.experiment(1).DateTime.Text=DateTime;
s.experiment(1).duration.Text=sprintf('%g',duration);

% Generate default groups
if(~exist('old','var') || isempty(old)), old=false; end

% load electrode types from file
[electrode_names,electrode_configuration]=load_electrode_types('electrodes.cfg');

% Generate uicontrols for electrode selection
if(~exist('electrode','var') || isempty(electrode))
    scrsz = get(0,'ScreenSize'); width = scrsz(3)*.2; height = scrsz(4)*.1;
    left=(scrsz(3)-width)*.5; bot=(scrsz(4)-height)*.5;
    dial=dialog('position',[left bot width height]);
    label_width=.8;
    margin=(1-label_width)/2;
    uicontrol('Style','text','String','Please select electrode type',...
        'fontsize',16,'fontweight','bold','horizontalalignment','left',...
        'Units','normalized','position',[margin .5 label_width .4],...
        'Parent',dial);
    menu_width=.8;
    margin=(1-menu_width)/2;
    menu=uicontrol('Style', 'popup',...
        'String',electrode_names,'fontsize',12,...
        'Units','normalized','Position',[margin 0 menu_width .6],'Parent',dial);
    button_width=.2;
    margin=(1-button_width)/2;
    uicontrol('Style','pushbutton','String','OK','fontsize',16,'Units','normalized',...
        'Position',[margin .05 button_width .25],'Parent',dial,...
        'Callback',{@set_electrode_type,menu});
    waitfor(dial);
else
    electrode_type=electrode;
end

loc = strcmp(electrode_type,electrode_names);
channels=electrode_configuration{loc};

% this is a fix due to a messup in the Kaplan lab involving cables. Nobody
% else has to worry about it
if(strcmp(electrode_type,'Buszaki 8-shank (64ch)') && old)
    channels(1:4)=cellfun(@(x) x-32,channels(1:4),'UniformOutput',false);
    channels(5:8)=cellfun(@(x) x+32,channels(5:8),'UniformOutput',false);
end

[~,gains]=plx_adchan_gains(plx_file);

% create shank+channel structure
num_shanks=length(channels);
for g = 1:num_shanks
    s.experiment.shank{g}.Attributes.num=sprintf('%g',g);
    num_channels=length(channels{g});
    for c = 1:num_channels
        current_channel=channels{g}(c);
        s.experiment.shank{g}.channel{c}.threshold.Text='-5';
        s.experiment.shank{g}.channel{c}.Attributes.num=sprintf('%g',current_channel);
        s.experiment.shank{g}.channel{c}.gain.Text=sprintf('%g',gains(current_channel));
        
        % channel is disabled/does not exist
        if(isempty(find(full_channels==current_channel,1,'first')))
            s.experiment.shank{g}.channel{c}.enabled.Text='false';
        else
            s.experiment.shank{g}.channel{c}.enabled.Text='true';
        end
    end
end

% If we have to swap 1:32 with 33:64
if(old), s.experiment.swap32.Text='1';
else s.experiment.swap32.Text='0';
end

s.experiment.electrode.Text=electrode_type;

% put field names in order
fieldnames = {'plx_file','mat_file','electrode','ADRate','DateTime','duration','swap32','shank'};
s.experiment=orderfields(s.experiment,fieldnames);


% write result to xml file
struct2xml(s,xml_file);

function set_electrode_type(~,~,menu_obj)
global electrode_type;
s=get(menu_obj,'String');
electrode_type=s{get(menu_obj,'Value')};
%setappdata(0,'electrode_type',result);
delete(get(menu_obj,'parent'));