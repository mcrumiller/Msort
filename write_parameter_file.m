function output_file = write_parameter_file(filename, base_name, groups, Fs, num_channels, t, num_features, voltage_range)
%--------------------------------------------------------------------------
% write_parameter_file.m - 
%
% Usage: write_parameter_file();
%
% Input:  
% Output: 
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
output_file = sprintf('%s.xml', base_name);
if(output_file==-1), return; end
f = fopen(output_file,'w');

% get necessary data
%[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, ...
%    SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, ...
%    DateTime] = plx_information(filename); %#ok<NASGU,ASGLU>
gain=1;

fprintf(f,'<parameters>\r\n');
fprintf(f,'  <acquisitionSystem>\r\n');
fprintf(f,'    <nBits>16</nBits>\r\n');
fprintf(f,'    <nChannels>%g</nChannels>\r\n',num_channels);
fprintf(f,'    <samplingRate>%g</samplingRate>\r\n',Fs);
fprintf(f,'    <voltageRange>%g</voltageRange>\r\n', voltage_range); % I think this is the range of data
fprintf(f,'    <amplification>%g</amplification>\r\n',gain); % what if it's different for each channel?
fprintf(f,'    <offset>0</offset>\r\n');
fprintf(f,'  </acquisitionSystem>\r\n');
%fprintf(f,'  <fieldPotentials>\r\n');
%fprintf(f,'    <lfpSamplingRate>1250</lfpSamplingRate>\r\n'); % not sure about this one--do we have LFP?
%fprintf(f,'  </fieldPotentials>\r\n');
fprintf(f,'  <anatomicalDescription>\r\n');
fprintf(f,'    <channelGroups>\r\n');

for i = 1:length(groups)
    group=groups{i};
    fprintf(f,'      <group>\r\n');
    for j = 1:length(group)
        fprintf(f,'        <channel>%g</channel>\r\n',group(j));
    end
    fprintf(f,'      </group>\r\n');
end
fprintf(f,'    </channelGroups>\r\n');

fprintf(f,'  </anatomicalDescription>\r\n');
fprintf(f,'  <spikeDetection>\r\n');
fprintf(f,'    <channelGroups>\r\n');

for i = 1:length(groups)
    group = groups{i};
    fprintf(f,'      <group>\r\n');
    fprintf(f,'        <channels>\r\n');
    for j = 1:length(group)
        fprintf(f,'          <channel>%g</channel>\r\n',group(j));
    end
    fprintf(f,'        </channels>\r\n');
    fprintf(f,'        <nSamples>%g</nSamples>\r\n', length(t)); % points per sample waveform
    fprintf(f,'        <peakSampleIndex>%g</peakSampleIndex>\r\n',find(t==0));
    fprintf(f,'        <nFeatures>%g</nFeatures>\r\n',num_features);    % num features
    fprintf(f,'      </group>\r\n');
end
fprintf(f,'    </channelGroups>\r\n');
fprintf(f,'  </spikeDetection>\r\n');
fprintf(f,'</parameters>');

fclose(f);