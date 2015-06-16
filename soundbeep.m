function soundbeep(num_beeps, freqs)
%--------------------------------------------------------------------------
% beep.m - Generate beeping sounds.
%
% Usage: beep;                    % default 3 beeps, 440 tone
%        beep(num_beeps,freqs);
%
% Input:  num_beeps               * number of beeps to play
%         freqs                   * frequencies of notes to be played
% Output: <beeping sound>
%
% Example: % produce 3 beeps of an A major7 chord
%          >> soundbeep(3,[440 554.37 659.26 830.61]);
%
% Written by Marshall Crumiller
% email: marshall.crumiller@mssm.edu
%--------------------------------------------------------------------------
if(~exist('num_beeps','var') || isempty(num_beeps)), num_beeps=3; end
if(~exist('freqs','var') || isempty(freqs)), freqs=523.25; end

Fs = 16000;          % sampling rate
beep_duration = .15; % in seconds
pause_duration = .1; % in seconds

% set up timing for beeps and pauses
total_duration=beep_duration*num_beeps + pause_duration*(num_beeps-1);
t=linspace(0,total_duration,Fs*total_duration);
numpoints_beep=beep_duration*Fs;
numpoints_pause=pause_duration*Fs;
silent_locs=repmat(numpoints_beep+1:numpoints_beep+numpoints_pause,num_beeps-1,1);
adder=repmat((0:(num_beeps-2))',1,numpoints_pause)*(numpoints_beep+numpoints_pause);
silent_locs=silent_locs+adder;
silent_locs=reshape(silent_locs',1,[]);
t(silent_locs)=0;

% create sum of frequencies
x=repmat(t,length(freqs),1)*2*pi;
freqs=repmat(freqs',1,length(t));
y=sum(sin(x.*freqs),1); y=y./max(y); % set max = 1

% play the sound
sound(y,Fs);