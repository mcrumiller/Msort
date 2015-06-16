function output_file = write_spike_file(waveforms, base_name, number)
%--------------------------------------------------------------------------
% write_spike_file.m - Given a set of waveforms extracted from
% get_spikes_from_plx.m, writes the results out to a binary .spk file,
% whose format is described at:
%
%  http://klusters.sourceforge.net/UserManual/data-files.html#data_files
%
% Usage: write_spike_file();
%
% Input:  waveforms             * CxWxT matrix of waveform data
%         base_name             * base filename of output file
% Output: output_file           * name of output file written.  Returns -1
%                                 if failure occurs.
%
% Written by Marshall Crumiller
%--------------------------------------------------------------------------
if(~exist('number','var')), number=1; end
output_file = sprintf('%s.spk.%g', base_name,number);
if(output_file==-1), return; end
f = fopen(output_file,'w');

waveforms = permute(waveforms,[1 3 2]);
fwrite(f,waveforms,'int16');
fclose(f);
