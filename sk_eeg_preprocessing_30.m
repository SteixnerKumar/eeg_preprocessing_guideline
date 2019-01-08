% This function is the pre-processing protocol for the EEG
% Involves the following parts
% 1.) Reading the file
% 2.) Downsampling
% 3.) cut the data between the start and end triggers
% 4.) Divide the data into two participants
% 5.) Checking for the Bridgeing
% 6.) Filtering
% 7.) Reject bad channels (and interpolate the good ones)
% 8.) Re-referencing
% cut into epoch, baseline correction
% 9.) Artifact correction using ICA (perform rank reduction before becauuse of interpolation)
% 10.)
%


%% basic required parameters
tt.version = 'enemy'; % 'single', 'enemy', 'friend'
tt.session = 1; % session number
tt.sub_A = 30100; %
tt.sub_B = 30300; %
%
plots_wanted = 0; % 1 to check the filter
plots_see_chan = 15; % the channel number to see
fs_new = 512; % new sampling frequency
filter_bp = [0.05 45]; % the band-pass filter cutoff frequencies
n_chan = 128; % for each participants number of channels

% downsampling, band-pass, re-ref, artifact-correction, sectioning, baseline-correction
%% Reading data
fprintf('Loading the data ...\n ');
filename.eegdata = strcat('tt_',tt.version,'_session_',num2str(tt.session),'_vp',num2str(tt.sub_A),'_vp',num2str(tt.sub_B),'.bdf');
data.both = ft_read_data(filename.eegdata);
hdr.both = ft_read_header(filename.eegdata);
events = ft_read_event(filename.eegdata);
% the behavioral data
if ~strcmp(tt.version,'single')
    tt.version = strcat('multi_',tt.version);
end
filename.behavior_A= strcat('tt_',num2str(tt.version),'_sub',num2str(tt.sub_A),'_session',num2str(tt.session),'.mat');
load(filename.behavior_A);behave_data.A.out = out;behave_data.A.ttsk = ttsk;
filename.behavior_B= strcat('tt_',num2str(tt.version),'_sub',num2str(tt.sub_B),'_session',num2str(tt.session),'.mat');
load(filename.behavior_B);behave_data.B.out = out;behave_data.B.ttsk = ttsk;
clear out ttsk
fprintf('Data loaded.\n ');


%% downsampling
% maybe use the anti-aliasing low pass filter before the downsampling process (yet to be imnplemented)
% downsample is needed to save the computation power and time --
% needs to be done on the data, in the header and in the events.

% In the data
fprintf('Downsampling the data by a factor of %d ...\n ',hdr.both.Fs/fs_new);
if ~(floor(hdr.both.Fs/fs_new) == hdr.both.Fs/fs_new)
    error('new sampling rate needs to be a proper divisor of original sampling rate');
end
data_both_temp = NaN(size(data.both,1),ceil(size(data.both,2)/(hdr.both.Fs/fs_new))); % initialization of a temporary variable
fprintf('individual channel-wise out of %d\n',hdr.both.nChans);
for loop_chan = 1:hdr.both.nChans % filter loop as channel-wise downsampling will be done
    data_both_temp(loop_chan,:) = downsample(data.both(loop_chan,:),(hdr.both.Fs/fs_new));
    if (mod(loop_chan,30)==1) && (loop_chan~=1)  % max 30 channels to show in one line
        fprintf('\n ');
    end
    fprintf('%d ',loop_chan);
end
fprintf('\n ');
% fprintf('Downsampling of the data done.\n ');
data.both = data_both_temp;
clear data_both_temp;
% In the events
for loop_markers = 1:length(strcmp('STATUS', {events.type}))
    if strcmp('STATUS', {events(loop_markers).type})
        events(loop_markers).sample = floor(events(loop_markers).sample/(hdr.both.Fs/fs_new));
    end
end
% In the header
hdr.both.nSamples = ceil(hdr.both.nSamples/(hdr.both.Fs/fs_new));
hdr.both.Fs = hdr.both.Fs/(hdr.both.Fs/fs_new);
clear fs_new loop_chan loop_markers;
fprintf('Downsampling of the data done.\n ');


%% Cut between the start and end
fprintf('Cutting the data between the start and end of the experiment on the basis of the defined start and end triggers ...\n');
% to get the trigger values
trig.start = behave_data.A.ttsk.trig.start_experiment; % to cut the data start trigger
trig.end = behave_data.A.ttsk.trig.end_experiment; % to cut the data end trigger
%
events_trigger = events(strcmp({events.type}, 'STATUS'));
byte.second = nan(1,size(events_trigger,2));
for loop_triggers = 1:size(events_trigger,2)
    in_binary = dec2bin(events_trigger(loop_triggers).value,24)-'0';
    binary_usb = in_binary(9:16);
    byte.second(loop_triggers) = bin2dec(num2str(binary_usb));
end
clear in_binary loop_try binary_lsb binary_usb binary_third
% to get the start and the end triggers (tiger_task: start = 100, end  = 101)
temp_start = find(byte.second == trig.start);
temp_end = find(byte.second == trig.end);
if isempty(temp_start)
    temp_start = 1;
    fprintf('no start trigger found hence using the starting of the recording as the start point\n');
end
if isempty(temp_end)
    temp_end = length(events_trigger);
    fprintf('no end trigger found hence using the end of the recording as the end point\n');
end
% to select the triggers only between the start and the end of the experiment
% cut the events
events_trigger = events_trigger(temp_start:temp_end);
% cut the eeg potential data now
data.both = data.both(:,events_trigger(1).sample:events_trigger(end).sample);
% cut the header
hdr.both.nSamples = size(data.both,2);
clear temp_start temp_end loop_triggers byte
%
fprintf('Cutting the data between the start and end of the experiment on the basis of the defined start and end triggers ... Done.\n');


%% Cut into the two participants
fprintf('Cutting the data into the two participants ...\n');
temp_chans_A = nan(hdr.both.nChans,1);
temp_chans_B = nan(hdr.both.nChans,1);
for loop_chans = 1:hdr.both.nChans
    if (str2double(hdr.both.label{loop_chans}(1))==1) % first participant labels start with '1'
        temp_chans_A(loop_chans) = 1;
    elseif (str2double(hdr.both.label{loop_chans}(1))==2) % first participant labels start with '2'
        temp_chans_B(loop_chans) = 1;
    end
end
temp_chans_A = find(temp_chans_A==1);
temp_chans_B = find(temp_chans_B==1);
% setting header file for both participants
hdr.A = hdr.both;
hdr.B = hdr.both;
hdr.A.label = hdr.both.label(temp_chans_A);
hdr.B.label = hdr.both.label(temp_chans_B);
hdr.A.nChans = length(temp_chans_A);
hdr.B.nChans = length(temp_chans_B);
hdr.A.chantype = hdr.both.chantype(temp_chans_A);
hdr.B.chantype = hdr.both.chantype(temp_chans_B);
hdr.A.chanunit = hdr.both.chanunit(temp_chans_A);
hdr.B.chanunit = hdr.both.chanunit(temp_chans_B);
% Cutting the eeg potential data into the participants
data.A = data.both(temp_chans_A,:);
data.B = data.both(temp_chans_B,:);
clear temp_chans_A temp_chans_B loop_chans
fprintf('Cutting the data into the two participants ... Done.\n');

%% Check for Bridging
A.channels_to_remove_bridge = zeros(hdr.A.nChans,1);
B.channels_to_remove_bridge = zeros(hdr.B.nChans,1);
%  [EB_out ,ED_out] = eBridge(ALLEEG(2),{'1-EXG1','1-EXG2','1-EXG3','1-EXG4','1-EXG5','1-EXG6','1-EXG7','1-EXG8','2-EXG1','2-EXG2','2-EXG3','2-EXG4','2-EXG5','2-EXG6','2-EXG7','2-EXG8','Status'});
fprintf('Checking for bridging ...\n');
% arranging data for the use by ebridge algorithm
A.EEG_wrapper.data = data.A;
A.EEG_wrapper.srate = hdr.A.Fs;
A.EEG_wrapper.chanlocs = struct('labels',hdr.A.label);
B.EEG_wrapper.data = data.B;
B.EEG_wrapper.srate = hdr.B.Fs;
B.EEG_wrapper.chanlocs = struct('labels',hdr.B.label);
% sk_ebridge is without plots
fprintf('Checking for bridges in the first participant ...\n');
[A.EB_out ,A.ED_out] = sk_eBridge(A.EEG_wrapper,{'1-EXG1','1-EXG2','1-EXG3','1-EXG4','1-EXG5','1-EXG6','1-EXG7','1-EXG8', ...
    '2-EXG1','2-EXG2','2-EXG3','2-EXG4','2-EXG5','2-EXG6','2-EXG7','2-EXG8','Status'});
fprintf('Checking for bridges in the second participant ...\n');
[B.EB_out ,B.ED_out] = sk_eBridge(B.EEG_wrapper,{'1-EXG1','1-EXG2','1-EXG3','1-EXG4','1-EXG5','1-EXG6','1-EXG7','1-EXG8', ...
    '2-EXG1','2-EXG2','2-EXG3','2-EXG4','2-EXG5','2-EXG6','2-EXG7','2-EXG8','Status'});
% clear A B
fprintf('Checking for bridging ... Done.\n');
% remove the bridged channels
fprintf('Marking the bridged channels for both participants ...\n');
% data.A(A.EB_out.Bridged.Indices,:) = [];
% data.B(B.EB_out.Bridged.Indices,:) = [];
A.channels_to_remove_bridge(A.EB_out.Bridged.Indices,1) = true;
B.channels_to_remove_bridge(B.EB_out.Bridged.Indices,1) = true;
fprintf('Marking the bridged channels for both participants ... Done.\n');

%% Filtering (highpass, then lowpass)
% pre filter plot
if plots_wanted
    figure;
    subplot(2,2,1);plot(data.A(plots_see_chan,:));title('The unfiltered data time');
    [XFreqRange, YAmplitude] = sk_dofft(data.A(plots_see_chan,:), hdr.A.Fs, 4);
    subplot(2,2,2);plot(XFreqRange, YAmplitude);axis([5 100 -0.2 1.5]);title('The unfiltered data freq');
    clear XFreqRange YAmplitude;
end

fprintf('Filtering the EEG potentials for both participants ...\n');
filter.type = 'but'; % the filter type ('but' is for butterworth)
filter.dir = 'twopass'; % the filter direction ('twopass' is default for both direction)
filter.order = 4; % filter order number
filter.lowpass = 45; % in Hz
filter.highpass = 1; % in Hz
fprintf('Using filter: %s, with direction: %s, order: %d, in range: %d - %d Hz.\n',filter.type,filter.dir,filter.order,filter.highpass,filter.lowpass);
[data.Af] = ft_preproc_bandpassfilter(data.A, hdr.A.Fs, [filter.highpass filter.lowpass], filter.order, filter.type, filter.dir);
[data.Bf] = ft_preproc_bandpassfilter(data.B, hdr.B.Fs, [filter.highpass filter.lowpass], filter.order, filter.type, filter.dir);
fprintf('Filtering the EEG potentials for both participants ... Done.\n');

% post filter plot
if plots_wanted
    subplot(2,2,3);plot(data.Af(plots_see_chan,:));title('The filtered data time');
    [XFreqRange, YAmplitude] = sk_dofft(data.Af(plots_see_chan,:), hdr.A.Fs, 4);
    subplot(2,2,4);plot(XFreqRange, YAmplitude);axis([5 100 -0.2 1.5]);title('The filtered data freq');
    clear XFreqRange YAmplitude;
end


%% reject bad channels
% using cleanraw from the eeglab
% a) flatline channels
% b) noisy channels
% c) short-time bursts
% d) incompletely repaird segments from the data
fprintf('Cleaning data ...\n');

% find(A.channels_to_remove)'
% find(B.channels_to_remove)'

% a) flatline channels (based: clean_flatlines)
fprintf('---Checking and marking for flat-line ...\n');
max_flatline_duration = 5; % in seconds
max_allowed_jitter = 20; % as multiples of epsilon
%
A.channels_to_remove_flatline = zeros(hdr.A.nChans,1);
B.channels_to_remove_flatline = zeros(hdr.B.nChans,1);
temp_A.data = data.Af;temp_A.srate = hdr.A.Fs;temp_A.nbchan = hdr.A.nChans;
temp_B.data = data.Bf;temp_B.srate = hdr.B.Fs;temp_B.nbchan = hdr.B.nChans;
%
temp_A = sk_clean_flatlines(temp_A,max_flatline_duration,max_allowed_jitter);
temp_B = sk_clean_flatlines(temp_B,max_flatline_duration,max_allowed_jitter);
%
A.channels_to_remove_flatline = temp_A.removed_channels;
B.channels_to_remove_flatline = temp_B.removed_channels;
clear temp_A temp_B
fprintf('---Checking and marking for flat-line channels ... Done.\n');


% b) noisy channels (remove noisy channels by correlation and line-noise thresholds) (based: clean_channels)
fprintf('---Checking and marking for noisy channels ...\n');
min_corr = 0.45; % minimum correlation
ignored_quantile = 0.1; % check the function description
window_len = 2; % in seconds
max_broken_time = 0.5; % in seconds
linenoise_aware = true; % if there is line noise 
%
A.channels_to_remove_noisy = zeros(hdr.A.nChans,1);
B.channels_to_remove_noisy = zeros(hdr.B.nChans,1);
temp_A.data = data.Af;temp_A.srate = hdr.A.Fs;temp_A.nbchan = hdr.A.nChans;
temp_B.data = data.Bf;temp_B.srate = hdr.B.Fs;temp_B.nbchan = hdr.B.nChans;
%
temp_A = sk_clean_channels_nolocs(temp_A,min_corr,ignored_quantile,window_len,max_broken_time,linenoise_aware);
temp_B = sk_clean_channels_nolocs(temp_B,min_corr,ignored_quantile,window_len,max_broken_time,linenoise_aware);
% flag channels to be implemented
A.channels_to_remove_noisy = temp_A.removed_channels;
B.channels_to_remove_noisy = temp_B.removed_channels;
clear temp_A temp_B
%
fprintf('---Checking and marking for noisy channels ... Done.\n');
%

%%% Here, remove the bad channels before the ASR and the clean windows
fprintf('Removing the bad channels before the ASR ...\n');
A.channels_to_remove =  A.channels_to_remove_bridge | A.channels_to_remove_flatline' | A.channels_to_remove_noisy;
B.channels_to_remove = B.channels_to_remove_bridge | B.channels_to_remove_flatline' | B.channels_to_remove_noisy;
%
data.Af_cut = data.Af(~A.channels_to_remove,:);
data.Bf_cut = data.Bf(~B.channels_to_remove,:);
fprintf('Removing the bad channels before the ASR ... Done.\n');
%%%

% c) short-time bursts (based: clean_asr)
fprintf('---Implementing the artifact sub-space reconstruction (ASR) for short time bursts ...\n');
% use the defaults
temp_A.data = data.Af_cut;temp_A.srate = hdr.A.Fs;temp_A.nbchan = size(data.Af_cut,1);
temp_B.data = data.Bf_cut;temp_B.srate = hdr.B.Fs;temp_B.nbchan = size(data.Bf_cut,1);
%
temp_A = sk_clean_asr(temp_A);
temp_B = sk_clean_asr(temp_B);
%
data.Af_cut_asr =  temp_A.data;
data.Bf_cut_asr =  temp_B.data;
clear temp_A temp_B
fprintf('---Implementing the artifact sub-space reconstruction (ASR) for short time bursts ... Done.\n');

% d) incompletely repaird segments from the data
% use asr data, not filtered data ()
fprintf('---Remove periods with abnormally high-power content from continuous data ...\n');
window_crit = 0.25;
window_crit_tolerances = [-inf 7];
% 
temp_A.data = data.Af_cut_asr;temp_A.srate = hdr.A.Fs;temp_A.nbchan = size(data.Af_cut,1);
temp_B.data = data.Bf_cut_asr;temp_B.srate = hdr.B.Fs;temp_B.nbchan = size(data.Bf_cut,1);
%
temp_A = sk_clean_windows(temp_A,window_crit,window_crit_tolerances); 
temp_B = sk_clean_windows(temp_B,window_crit,window_crit_tolerances); 
%
data.Af_cut_asr_repaired =  temp_A.data;
data.Bf_cut_asr_repaired =  temp_B.data;
clear temp_A temp_B
fprintf('---Remove periods with abnormally high-power content from continuous data ... Done.\n');
%


%% interpolating the missing channels
% Interpolate channels. (not ideal before ICA, but better for re-referencing)
% get EOG channels
A.channels_eog =  zeros(hdr.A.nChans,1);
B.channels_eog =  zeros(hdr.B.nChans,1);
A.channels_eog(contains(hdr.A.label,'EX')) = 1;
B.channels_eog(contains(hdr.B.label,'EX')) = 1;
%
% add fake channels with zeros that need to be interpolated
data.A_temp_interp = zeros(hdr.A.nChans,size(data.Af_cut_asr_repaired,2));
[f_idx,~]=find(A.channels_to_remove);
count_chan = 1;
for loop_fake_chan  =1:hdr.A.nChans
    if ~any(loop_fake_chan == f_idx)
        data.A_temp_interp(loop_fake_chan,:) = data.Af_cut_asr_repaired(count_chan,:);
        count_chan = count_chan+1;
    end
end
clear f_idx count_chan loop_fake_chan
data.B_temp_interp = zeros(hdr.B.nChans,size(data.Bf_cut_asr_repaired,2));
[f_idx,~]=find(B.channels_to_remove);
count_chan = 1;
for loop_fake_chan  =1:hdr.B.nChans
    if ~any(loop_fake_chan == f_idx)
        data.B_temp_interp(loop_fake_chan,:) = data.Bf_cut_asr_repaired(count_chan,:);
        count_chan = count_chan+1;
    end
end
clear f_idx count_chan loop_fake_chan

% remove the channels that have EOG 
data.A_temp_interp = data.A_temp_interp(~A.channels_eog,:);
data.B_temp_interp = data.B_temp_interp(~B.channels_eog,:);
% note the bad channels except the EOG for interpolation
bad_chans_A = A.channels_to_remove(~A.channels_eog)'; [~,bad_chans_A]= find(bad_chans_A);
bad_chans_B = B.channels_to_remove(~B.channels_eog)'; [~,bad_chans_B]= find(bad_chans_B);
% get electrode positions
filename.elec = '128 EQ A on MNI average head with ABCD labels for BESA.sfp';
elec = readlocs(filename.elec);
EEG.chanlocs = elec(4:end);
%
method = 'spherical';
%'
EEG.trials = 1; % no trials
EEG.data = data.A_temp_interp;
EEG.nbchan = size(data.A_temp_interp,1);
EEG.pnts = size(data.Af_cut_asr_repaired,2);
A.EEG_interp = sk_eeg_interp(EEG, bad_chans_A, method); % main function
%
EEG.trials = 1; % no trials
EEG.data = data.B_temp_interp;
EEG.nbchan = size(data.B_temp_interp,1);
EEG.pnts = size(data.Bf_cut_asr_repaired,2);
B.EEG_interp = sk_eeg_interp(EEG, bad_chans_B, method); % main function
%
data.Af_cut_asr_repaired_interp = A.EEG_interp.data;
data.Bf_cut_asr_repaired_interp = B.EEG_interp.data;
clear EEG method

%% Re-referencing
fprintf('re-referenciong the data to the common average reference (CAR)...\n');
mean_data_one = mean(A.EEG_interp.data,1);
mean_data_two = mean(B.EEG_interp.data,1);
A.EEG_interp.data_car = A.EEG_interp.data - repmat(mean_data_one,size(A.EEG_interp.data,1),1);
B.EEG_interp.data_car = B.EEG_interp.data - repmat(mean_data_two,size(B.EEG_interp.data,1),1);
data.Af_cut_asr_repaired_interp_car = A.EEG_interp.data_car;
data.Bf_cut_asr_repaired_interp_car = B.EEG_interp.data_car;
clear mean_data_one mean_data_two
fprintf('re-referenciong the data to the common average reference (CAR)... Done.\n');

%% cut the data here into trials/epochs 
fprintf('Cutting the data into epochs/trials ...\n');
flag_bc_trial = 1; % for baseline correction , 0 not not
if flag_bc_trial
    temp_baseline_correction_time = 200; % time in milliseconds
    temp_baseline_correction_time = ceil((temp_baseline_correction_time/1000)*hdr.both.Fs); % in points
end
% to get the trigger values
events_trigger = events(strcmp({events.type}, 'STATUS'));
byte.second = nan(1,size(events_trigger,2));
for loop_triggers = 1:size(events_trigger,2)
    in_binary = dec2bin(events_trigger(loop_triggers).value,24)-'0';
    binary_usb = in_binary(9:16);
    byte.second(loop_triggers) = bin2dec(num2str(binary_usb));
end
%%%
% trial.type = {'prediction','choice'};
% trial.type{1}%trial.type{2}
% length(trial.type)
% trial.prediction = [behave_data.A.ttsk.trig.start_prediction behave_data.A.ttsk.trig.end_prediction];
%%%
% trial sections
trial.num = length(find(byte.second == behave_data.A.ttsk.trig.start_evidence_own));
trial.prediction = [behave_data.A.ttsk.trig.start_prediction behave_data.A.ttsk.trig.end_prediction];
trial.choice = [behave_data.A.ttsk.trig.start_choice behave_data.A.ttsk.trig.end_choice];
trial.evidence_other = [behave_data.A.ttsk.trig.start_evidence_other_player behave_data.A.ttsk.trig.end_evidence_other_player];
trial.evidence_own = [behave_data.A.ttsk.trig.start_evidence_own behave_data.A.ttsk.trig.end_evidence_own];
trial.response_prediction_A = [behave_data.A.ttsk.trig.start_prediction behave_data.A.ttsk.trig.response_prediction_A];
trial.response_choice_A = [behave_data.A.ttsk.trig.start_choice behave_data.A.ttsk.trig.response_choice_A];
trial.response_prediction_B = [behave_data.A.ttsk.trig.start_prediction behave_data.A.ttsk.trig.response_prediction_B];
trial.response_choice_B = [behave_data.A.ttsk.trig.start_choice behave_data.A.ttsk.trig.response_choice_B];
% get the points in the data
trial.points.prediction = [[events_trigger(byte.second == trial.prediction(1)).sample];[events_trigger(byte.second == trial.prediction(2)).sample]]';
trial.points.choice = [[events_trigger(byte.second == trial.choice(1)).sample];[events_trigger(byte.second == trial.choice(2)).sample]]';
trial.points.evidence_other = [[events_trigger(byte.second == trial.evidence_other(1)).sample];[events_trigger(byte.second == trial.evidence_other(2)).sample]]';
trial.points.evidence_own = [[events_trigger(byte.second == trial.evidence_own(1)).sample];[events_trigger(byte.second == trial.evidence_own(2)).sample]]';
trial.points.response_prediction_A = [[events_trigger(byte.second == trial.response_prediction_A(1)).sample];[events_trigger(byte.second == trial.response_prediction_A(2)).sample]]';
trial.points.response_choice_A = [[events_trigger(byte.second == trial.response_choice_A(1)).sample];[events_trigger(byte.second == trial.response_choice_A(2)).sample]]';
trial.points.response_prediction_B = [[events_trigger(byte.second == trial.response_prediction_B(1)).sample];[events_trigger(byte.second == trial.response_prediction_B(2)).sample]]';
trial.points.response_choice_B = [[events_trigger(byte.second == trial.response_choice_B(1)).sample];[events_trigger(byte.second == trial.response_choice_B(2)).sample]]';
% cutting the data
% rowsXcolumnXtrials
data.trials.prediction.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.prediction(:,2)-trial.points.prediction(:,1))+1,trial.num);
data.trials.choice.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.choice(:,2)-trial.points.choice(:,1))+1,trial.num);
data.trials.evidence_other.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.evidence_other(:,2)-trial.points.evidence_other(:,1))+1,trial.num);
data.trials.evidence_own.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.evidence_own(:,2)-trial.points.evidence_own(:,1))+1,trial.num);
data.trials.response_prediction_A.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.response_prediction_A(:,2)-trial.points.response_prediction_A(:,1))+1,trial.num);
data.trials.response_choice_A.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.response_choice_A(:,2)-trial.points.response_choice_A(:,1))+1,trial.num);
data.trials.response_prediction_B.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.response_prediction_B(:,2)-trial.points.response_prediction_B(:,1))+1,trial.num);
data.trials.response_choice_B.data = nan(size(data.Af_cut_asr_repaired_interp_car,1),max(trial.points.response_choice_B(:,2)-trial.points.response_choice_B(:,1))+1,trial.num);
%
for loop_trials = 1:trial.num
    data.trials.prediction.data(:,1:(trial.points.prediction(loop_trials,2) - trial.points.prediction(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.prediction(loop_trials,1):trial.points.prediction(loop_trials,2));
    data.trials.choice.data(:,1:(trial.points.choice(loop_trials,2) - trial.points.choice(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.choice(loop_trials,1):trial.points.choice(loop_trials,2));
    data.trials.evidence_other.data(:,1:(trial.points.evidence_other(loop_trials,2) - trial.points.evidence_other(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.evidence_other(loop_trials,1):trial.points.evidence_other(loop_trials,2));
    data.trials.evidence_own.data(:,1:(trial.points.evidence_own(loop_trials,2) - trial.points.evidence_own(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.evidence_own(loop_trials,1):trial.points.evidence_own(loop_trials,2));
    data.trials.response_prediction_A.data(:,1:(trial.points.response_prediction_A(loop_trials,2) - trial.points.response_prediction_A(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.response_prediction_A(loop_trials,1):trial.points.response_prediction_A(loop_trials,2));
    data.trials.response_choice_A.data(:,1:(trial.points.response_choice_A(loop_trials,2) - trial.points.response_choice_A(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.response_choice_A(loop_trials,1):trial.points.response_choice_A(loop_trials,2));
    data.trials.response_prediction_B.data(:,1:(trial.points.response_prediction_B(loop_trials,2) - trial.points.response_prediction_B(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.response_prediction_B(loop_trials,1):trial.points.response_prediction_B(loop_trials,2));
    data.trials.response_choice_B.data(:,1:(trial.points.response_choice_B(loop_trials,2) - trial.points.response_choice_B(loop_trials,1)+1),loop_trials) = ...
        data.Af_cut_asr_repaired_interp_car(:,trial.points.response_choice_B(loop_trials,1):trial.points.response_choice_B(loop_trials,2));
end
% cutting the edges for equal trial lengths
data.trials.prediction.data = data.trials.prediction.data(:,1:min(trial.points.prediction(:,2)-trial.points.prediction(:,1))+1,:);
data.trials.choice.data = data.trials.choice.data(:,1:min(trial.points.choice(:,2)-trial.points.choice(:,1))+1,:);
data.trials.evidence_other.data = data.trials.evidence_other.data(:,1:min(trial.points.evidence_other(:,2)-trial.points.evidence_other(:,1))+1,:);
data.trials.evidence_own.data = data.trials.evidence_own.data(:,1:min(trial.points.evidence_own(:,2)-trial.points.evidence_own(:,1))+1,:);
data.trials.response_prediction_A.data = data.trials.response_prediction_A.data(:,1:min(trial.points.response_prediction_A(:,2)-trial.points.response_prediction_A(:,1))+1,:);
data.trials.response_choice_A.data = data.trials.response_choice_A.data(:,1:min(trial.points.response_choice_A(:,2)-trial.points.response_choice_A(:,1))+1,:);
data.trials.response_prediction_B.data = data.trials.response_prediction_B.data(:,1:min(trial.points.response_prediction_B(:,2)-trial.points.response_prediction_B(:,1))+1,:);
data.trials.response_choice_B.data = data.trials.response_choice_B.data(:,1:min(trial.points.response_choice_B(:,2)-trial.points.response_choice_B(:,1))+1,:);
%
fprintf('Cutting the data into epochs/trials ... Done.\n');

%%% Now baseline corect them
if flag_bc_trial
    fprintf('Baseline correction of the epochs/trials ...\n');
    % marking the points
    trial.points.bc.prediction = [trial.points.prediction(:,1)-temp_baseline_correction_time,trial.points.prediction(:,1)];
    trial.points.bc.choice = [trial.points.choice(:,1)-temp_baseline_correction_time,trial.points.choice(:,1)];
    trial.points.bc.evidence_other = [trial.points.evidence_other(:,1)-temp_baseline_correction_time,trial.points.evidence_other(:,1)];
    trial.points.bc.evidence_own = [trial.points.evidence_own(:,1)-temp_baseline_correction_time,trial.points.evidence_own(:,1)];
    trial.points.bc.response_prediction_A = [trial.points.response_prediction_A(:,1)-temp_baseline_correction_time,trial.points.response_prediction_A(:,1)];
    trial.points.bc.response_choice_A = [trial.points.response_choice_A(:,1)-temp_baseline_correction_time,trial.points.response_choice_A(:,1)];
    trial.points.bc.response_prediction_B = [trial.points.response_prediction_B(:,1)-temp_baseline_correction_time,trial.points.response_prediction_B(:,1)];
    trial.points.bc.response_choice_B = [trial.points.response_choice_B(:,1)-temp_baseline_correction_time,trial.points.response_choice_B(:,1)];
    % cutting and getting the mean
    trial.points.bc.mean.prediction = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.choice = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.evidence_other = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.evidence_own = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.response_prediction_A = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.response_choice_A = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.response_prediction_B = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    trial.points.bc.mean.response_choice_B = nan(size(data.Af_cut_asr_repaired_interp_car,1),trial.num);
    for loop_trials = 1:trial.num
        trial.points.bc.mean.prediction(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.prediction(loop_trials,1):trial.points.bc.prediction(loop_trials,2)),2);
        trial.points.bc.mean.choice(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.choice(loop_trials,1):trial.points.bc.choice(loop_trials,2)),2);
        trial.points.bc.mean.evidence_other(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.evidence_other(loop_trials,1):trial.points.bc.evidence_other(loop_trials,2)),2);
        trial.points.bc.mean.evidence_own(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.evidence_own(loop_trials,1):trial.points.bc.evidence_own(loop_trials,2)),2);
        trial.points.bc.mean.response_prediction_A(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.response_prediction_A(loop_trials,1):trial.points.bc.response_prediction_A(loop_trials,2)),2);
        trial.points.bc.mean.response_choice_A(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.response_choice_A(loop_trials,1):trial.points.bc.response_choice_A(loop_trials,2)),2);
        trial.points.bc.mean.response_prediction_B(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.response_prediction_B(loop_trials,1):trial.points.bc.response_prediction_B(loop_trials,2)),2);
        trial.points.bc.mean.response_choice_B(:,loop_trials) = mean(data.Af_cut_asr_repaired_interp_car(:,trial.points.bc.response_choice_B(loop_trials,1):trial.points.bc.response_choice_B(loop_trials,2)),2);
    end
    % substracting the mean from the data
    for loop_trials = 1:trial.num
        data.trials.prediction.data(:,:,loop_trials) = data.trials.prediction.data(:,:,loop_trials) - repmat(trial.points.bc.mean.prediction(:,loop_trials),1,size(data.trials.prediction.data,2));
        data.trials.choice.data(:,:,loop_trials) = data.trials.choice.data(:,:,loop_trials) - repmat(trial.points.bc.mean.choice(:,loop_trials),1,size(data.trials.choice.data,2));
        data.trials.evidence_other.data(:,:,loop_trials) = data.trials.evidence_other.data(:,:,loop_trials) - repmat(trial.points.bc.mean.evidence_other(:,loop_trials),1,size(data.trials.evidence_other.data,2));
        data.trials.evidence_own.data(:,:,loop_trials) = data.trials.evidence_own.data(:,:,loop_trials) - repmat(trial.points.bc.mean.evidence_own(:,loop_trials),1,size(data.trials.evidence_own.data,2));
        data.trials.response_prediction_A.data(:,:,loop_trials) = data.trials.response_prediction_A.data(:,:,loop_trials) - repmat(trial.points.bc.mean.response_prediction_A(:,loop_trials),1,size(data.trials.response_prediction_A.data,2));
        data.trials.response_choice_A.data(:,:,loop_trials) = data.trials.response_choice_A.data(:,:,loop_trials) - repmat(trial.points.bc.mean.response_choice_A(:,loop_trials),1,size(data.trials.response_choice_A.data,2));
        data.trials.response_prediction_B.data(:,:,loop_trials) = data.trials.response_prediction_B.data(:,:,loop_trials) - repmat(trial.points.bc.mean.response_prediction_B(:,loop_trials),1,size(data.trials.response_prediction_B.data,2));
        data.trials.response_choice_B.data(:,:,loop_trials) = data.trials.response_choice_B.data(:,:,loop_trials) - repmat(trial.points.bc.mean.response_choice_B(:,loop_trials),1,size(data.trials.response_choice_B.data,2));
    end
    fprintf('Baseline correction of the epochs/trials ... Done.\n');
end
%%%


% [data.trials.prediction.icaweights,data.trials.prediction.icasphere,data.trials.prediction.icameanvar,data.trials.prediction.icabias,data.trials.prediction.icasigns, ... 
%     data.trials.prediction.icalrates,data.trials.prediction.icadata,data.trials.prediction.icay] = runica(data.trials.prediction.data); % train using defaults 
% figure;erpimage(data.trials.prediction.data(29,:,:),[],1:size(data.trials.prediction.data,2));
% figure;spectopo(data.trials.prediction.data, 0, 512);
% % % EEG.trials = size(data.trials.prediction.data,3); % num trials
% % % EEG.data = data.trials.prediction.data;
% % % EEG.nbchan = size(data.trials.prediction.data,1);
% % % EEG.pnts = size(data.trials.prediction.data,2);EEG.srate = 512;
% % % EEG.xmin = 1;EEG.xmax = EEG.pnts;
% % % pop_erpimage(EEG,1);
% % % %
% % % EEG.trials = size(data.trials.choice.data,3); % no trials
% % % EEG.data = data.trials.choice.data;
% % % EEG.nbchan = size(data.trials.choice.data,1);
% % % EEG.pnts = size(data.trials.choice.data,2);EEG.srate = 512;
% % % EEG.xmin = 1;EEG.xmax = EEG.pnts;
% % % pop_erpimage(EEG,1);


%% remove the interpolated channels before ICA
data.Af_cut_asr_repaired_interp_car_cut = A.EEG_interp.data_car(~A.channels_to_remove(~A.channels_eog),:);
data.Bf_cut_asr_repaired_interp_car_cut = B.EEG_interp.data_car(~B.channels_to_remove(~B.channels_eog),:);

%% now the ICA
fprintf('performing the ICA using the runica algorithm...\n');
% [A.ica.icaweights,A.ica.icasphere] = runica(data.Af_cut_asr_repaired_interp_car_cut); % train using defaults 
[A.ica.icaweights,A.ica.icasphere,A.ica.icameanvar,A.ica.icabias,A.ica.icasigns,A.ica.icalrates,A.ica.icadata,A.ica.icay] = runica(data.Af_cut_asr_repaired_interp_car_cut); % train using defaults 
% [B.ica.icaweights,B.ica.icasphere] = runica(data.Bf_cut_asr_repaired_interp_car_cut); % train using defaults 
[B.ica.icaweights,B.ica.icasphere,B.ica.icameanvar,B.ica.icabias,B.ica.icasigns,B.ica.icalrates,B.ica.icadata,B.ica.icay] = runica(data.Bf_cut_asr_repaired_interp_car_cut); % train using defaults 
fprintf('performing the ICA using the runica algorithm... Done.\n');

%% Auto-detection of bad ICA components
fprintf('auto detecting the bad ICA components...\n');
% calculate the inverse of the ica weights
A.ica.icawinv = pinv( A.ica.icaweights*A.ica.icasphere );
B.ica.icawinv = pinv( B.ica.icaweights*B.ica.icasphere );
A.ica.data = data.Af_cut_asr_repaired_interp_car_cut;
A.ica.nbchan = size(A.ica.icaweights,1);
A.ica.trials = 1;
A.ica.times = 1:size(A.ica.data,2);
A.ica.srate = hdr.A.Fs;
A.ica.icachansind = 1:size(A.ica.icaweights,1);
B.ica.data = data.Bf_cut_asr_repaired_interp_car_cut;
B.ica.nbchan = size(B.ica.icaweights,1);
B.ica.trials = 1;
B.ica.times = 1:size(B.ica.data,2);
B.ica.srate = hdr.B.Fs;
B.ica.icachansind = 1:size(B.ica.icaweights,1);
% 
temp_channels_to_remove = A.channels_to_remove(~A.channels_eog);
A.ica.chanlocs = elec(4:end);A.ica.chanlocs = elec(~temp_channels_to_remove);
A.ica.icaact = A.ica.icaweights*A.ica.icasphere*A.ica.data;
temp_channels_to_remove = B.channels_to_remove(~B.channels_eog);
B.ica.chanlocs = elec(4:end);B.ica.chanlocs = elec(~temp_channels_to_remove);
B.ica.icaact = B.ica.icaweights*B.ica.icasphere*B.ica.data;
%
% automatic detection of ica components
cfg = [];
cfg.autocorr.enable = true;
cfg.focalcomp.enable = true;
% cfg.trialfoc.enable = true;
cfg.SNR.enable = true;
% cfg.resvar.enable = true; % need dipfit to be done before 
% cfg.EOGcorr.enable = true; % need EOG channels to do this
% cfg.chancorr.enable = true;
cfg.ADJUST.enable = true;
% cfg.FASTER.enable = true;
% cfg.MARA.enable = true; % not implememted
cfg.opts.noplot = 1; % '1' for no plots
%
[A.ica_reject, cfg] = sk_eeg_SASICA(A.ica,cfg);
[B.ica_reject, cfg] = sk_eeg_SASICA(B.ica,cfg);
%
% A.ica_reject.reject.gcompreject % to see the components to reject
% B.ica_reject.reject.gcompreject % to see the components to reject
fprintf('auto detecting the bad ICA components... Done.\n');

%% Removing the bad ICA components
fprintf('removing the bad ICA components...\n');
A.ica_clean.setname = 'ica_clean';
A.ica_clean.nbchan = A.ica.nbchan;
A.ica_clean.trials = A.ica.trials;
A.ica_clean.data = A.ica.data;
A.ica_clean.pnts = size(A.ica.times,2);
A.ica_clean.icawinv = A.ica.icawinv;
A.ica_clean.icasphere = A.ica.icasphere;
A.ica_clean.icaweights = A.ica.icaweights;
A.ica_clean.icachansind = A.ica.icachansind;
A.ica_clean.icaact = A.ica.icaweights*A.ica.icasphere*A.ica.icadata;  % Matrix multiplication
A.ica_clean.reject.gcompreject = A.ica_reject.reject.gcompreject;
A.ica_clean.reject.SASICA = A.ica_reject.reject.SASICA;
B.ica_clean.setname = 'ica_clean';
B.ica_clean.nbchan = B.ica.nbchan;
B.ica_clean.trials = B.ica.trials;
B.ica_clean.data = B.ica.data;
B.ica_clean.pnts = size(B.ica.times,2);
B.ica_clean.icawinv = B.ica.icawinv;
B.ica_clean.icasphere = B.ica.icasphere;
B.ica_clean.icaweights = B.ica.icaweights;
B.ica_clean.icachansind = B.ica.icachansind;
B.ica_clean.icaact = B.ica.icaweights*B.ica.icasphere*B.ica.icadata;  % Matrix multiplication
B.ica_clean.reject.gcompreject = B.ica_reject.reject.gcompreject;
B.ica_clean.reject.SASICA = B.ica_reject.reject.SASICA;
% Extract the Component numbers to remove from the output of SASICA in EEG.reject.gcompreject.  Wherever this is 1, that particular ICA component needs to be removed from the data.
A.ica_clean.remove_components = find(A.ica_reject.reject.gcompreject == 1);
B.ica_clean.remove_components = find(B.ica_reject.reject.gcompreject == 1);
% Remove the components from the main data without further ado!
A.ica_clean = sk_pop_subcomp(A.ica_clean, A.ica_clean.remove_components);
B.ica_clean = sk_pop_subcomp(B.ica_clean, B.ica_clean.remove_components);
fprintf('removing the bad ICA components... Done.\n');
% ica_clean.data is the cleaned data




% %%% plotting experiments
% % % % filename = '128 EQ A on MNI average head with ABCD labels for BESA.sfp';
% % % % elec = ft_read_sens(filename);  % for electrodes
% % % % temp_channels_to_remove = A.channels_to_remove(1:128);
% % % % elec.chanpos = elec.chanpos(4:end,:);% elec.chanpos = elec.chanpos(~temp_channels_to_remove,:);
% % % % elec.elecpos = elec.elecpos(4:end,:); %elec.elecpos = elec.elecpos(~temp_channels_to_remove,:);
% % % % elec.label = elec.label(4:end,1); %elec.label = elec.label(~temp_channels_to_remove,:);
% % % % elec.chanunit = elec.chanunit(4:end,1);% elec.chanunit = elec.chanunit(~temp_channels_to_remove,:);
% % % % elec.chantype = elec.chantype(4:end,1); %elec.chantype = elec.chantype(~temp_channels_to_remove,:);
% % % % 
% % % % cfg = [];
% % % % cfg.elec = elec;
% % % % layout = ft_prepare_layout(cfg);
% % % % ft_plot_lay(layout)
% % % % cfg = [];
% % % % cfg.layout = layout;
% % % % 
% % % % % avg_data.avg = data.Af_cut_asr_repaired;
% % % % avg_data.avg = A.EEG_interp.data;
% % % % %  avg_data.avg = A.EEG_interp.data_car
% % % % % avg_data.avg = A.EEG_interp.data_car-A.EEG_interp.data;
% % % % avg_data.fsample = hdr.A.Fs;
% % % % avg_data.label = elec.label;
% % % % avg_data.elec = elec;
% % % % avg_data.time = 0:1/hdr.A.Fs:(size(data.Af_cut_asr_repaired,2)/hdr.A.Fs)-(1/hdr.A.Fs);
% % % % avg_data.dimord = 'chan_time';
% % % % figure; ft_multiplotER(cfg, avg_data);
% % % % axis on % this shows the actual MATLAB axes
% % % % % ft_plot_topo3d(elec.chanpos,data.Af_cut_asr_repaired(1:121,1))
% % % % % cfg.xlim = [5 200];
% % % % %  ft_movieplotER(cfg, avg_data)





% % % % %% Band-pass filtering
% % % % % using an ideal - non realizable (non-causal) filter because the past, present and the
% % % % % future values are available and therefore better filter performance can
% % % % % be achieved (Alternative would be to use a causal realizable butterworth filter)
% % % % data_both_temp = NaN(size(data.both,1),size(data.both,2)); % initialization of a temporary variable
% % % % t = 0:1/hdr.both.Fs:hdr.both.nSamples/hdr.both.Fs-(1/hdr.both.Fs); % time in seconds
% % % % for filter_loop = 1:hdr.both.nChans % filter loop as channel-wise filtering will be done
% % % %     ts_data.both = timeseries(data.both(filter_loop,:),t); % time series to which the ideal fiter is applied
% % % %     mean_data.both = mean(data.both(filter_loop,:)); % mean of a particular channel
% % % %     ts_data.both = idealfilter(ts_data.both,filter_bp,'Pass'); % band pass filter used here.
% % % %     data_both_temp(filter_loop,:) = ts_data.both.data(:) + mean_data.both; % the output of the filter is added with the mean of the particular channel because the ideal filter has substracted the mean as it only works on the zero mean data
% % % %     clear ts_data.both mean_data.both; % clearing the unwanted variables
% % % % end
% % % % data.both= data_both_temp;
% % % % clear data_both_temp filter_loop t;
% % % % 
% % % % if plots_wanted == 1
% % % %     subplot(2,2,3);plot(data.both(plots_see_chan,:));
% % % %     [XFreqRange, YAmplitude] = sk_dofft(data.both(plots_see_chan,:), hdr.both.Fs, 4);
% % % %     subplot(2,2,4);plot(XFreqRange, YAmplitude);axis([5 55 -0.2 1.5]);
% % % %     clear XFreqRange YAmplitude;
% % % % end

% % % % if plots_wanted == 2
% % % %     subplot(4,1,3);plot(data.both(plots_see_chan,:));title('Down-sampled and Band-passed data');
% % % % end

% % % % %% Re-referencing the data
% % % % % Using the common-average-reference (i.e. substracting the mean from the entire data)
% % % % data_one = data.both(1:128,:); % data without the external electrodes
% % % % data_two = data.both(137:264,:); % data without the external electrodes
% % % % mean_data_one = mean(data_one,1);
% % % % mean_data_two = mean(data_two,1);
% % % % 
% % % % % car_data.both = mean(data.both,2);
% % % % car_data.both = (mean_data_one + mean_data_two)/2; % mean of the relevent electrodes
% % % % data.both = data.both - repmat(car_data.both,size(hdr.both.label,1),1);
% % % % clear car_data.both mean_data_one mean_data_two;
% % % % 
% % % % if plots_wanted == 2
% % % %     subplot(4,1,4);plot(data.both(plots_see_chan,:));title('Down-sampled, Band-passed and Re-referenced data');
% % % % end
% % % % 
% % % % %% Divide the data into two participants
% % % % % 1 to 136 - subj_one and 137 to 272 - subj_two. The last one is Status.
% % % % % data_one = data.both(1:128,:); % data without the external electrodes
% % % % % data_two = data.both(137:264,:); % data without the external electrodes
% % % % 
% % % % if plots_wanted == 3
% % % %     figure;subplot(4,1,1);plot(data_one(plots_see_chan,:));title('Subject_1; Down-sampled, Band-passed, Re-referenced'); %axis([256 5024 -200 300]);
% % % %     subplot(4,1,3);plot(data_two(plots_see_chan,:));title('Subject_2; Down-sampled, Band-passed, Re-referenced');%axis([256 5024 -200 300]);
% % % % end
% % % % 
% % % % data_one = data_one';
% % % % data_two = data_two';
% % % % 
% % % % %% Artifact-correction


%% extra functions
% % % % % % % % % % % % % % % % % % % % % function [mu,sig,alpha,beta] = fit_eeg_distribution(X,min_clean_fraction,max_dropout_fraction,quants,step_sizes,beta)
% % % % % % % % % % % % % % % % % % % % % % assign defaults
% % % % % % % % % % % % % % % % % % % % % if ~exist('min_clean_fraction','var') || isempty(min_clean_fraction)
% % % % % % % % % % % % % % % % % % % % %     min_clean_fraction = 0.25; end
% % % % % % % % % % % % % % % % % % % % % if ~exist('max_dropout_fraction','var') || isempty(max_dropout_fraction)
% % % % % % % % % % % % % % % % % % % % %     max_dropout_fraction = 0.1; end
% % % % % % % % % % % % % % % % % % % % % if ~exist('quants','var') || isempty(quants)
% % % % % % % % % % % % % % % % % % % % %     quants = [0.022 0.6]; end
% % % % % % % % % % % % % % % % % % % % % if ~exist('step_sizes','var') || isempty(step_sizes)
% % % % % % % % % % % % % % % % % % % % %     step_sizes = [0.01 0.01]; end
% % % % % % % % % % % % % % % % % % % % % if ~exist('beta','var') || isempty(beta)
% % % % % % % % % % % % % % % % % % % % %     beta = 1.7:0.15:3.5; end
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % sanity checks
% % % % % % % % % % % % % % % % % % % % % if ~isvector(quants) || numel(quants) > 2
% % % % % % % % % % % % % % % % % % % % %     error('Fit quantiles needs to be a 2-element vector (support for matrices deprecated).'); end
% % % % % % % % % % % % % % % % % % % % % if any(quants(:)<0) || any(quants(:)>1)
% % % % % % % % % % % % % % % % % % % % %     error('Unreasonable fit quantiles.'); end
% % % % % % % % % % % % % % % % % % % % % if any(step_sizes<0.0001) || any(step_sizes>0.1)
% % % % % % % % % % % % % % % % % % % % %     error('Unreasonable step sizes.'); end
% % % % % % % % % % % % % % % % % % % % % if any(beta>=7) || any(beta<=1)
% % % % % % % % % % % % % % % % % % % % %     error('Unreasonable shape range.'); end
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % sort data so we can access quantiles directly
% % % % % % % % % % % % % % % % % % % % % X = double(sort(X(:)));
% % % % % % % % % % % % % % % % % % % % % n = length(X);
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % calc z bounds for the truncated standard generalized Gaussian pdf and pdf rescaler
% % % % % % % % % % % % % % % % % % % % % for b=1:length(beta)    
% % % % % % % % % % % % % % % % % % % % %     zbounds{b} = sign(quants-1/2).*gammaincinv(sign(quants-1/2).*(2*quants-1),1/beta(b)).^(1/beta(b)); %#ok<*AGROW>
% % % % % % % % % % % % % % % % % % % % %     rescale(b) = beta(b)/(2*gamma(1/beta(b)));
% % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % determine the quantile-dependent limits for the grid search
% % % % % % % % % % % % % % % % % % % % % lower_min = min(quants);                    % we can generally skip the tail below the lower quantile
% % % % % % % % % % % % % % % % % % % % % max_width = diff(quants);                   % maximum width is the fit interval if all data is clean
% % % % % % % % % % % % % % % % % % % % % min_width = min_clean_fraction*max_width;   % minimum width of the fit interval, as fraction of data
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % get matrix of shifted data ranges
% % % % % % % % % % % % % % % % % % % % % X = X(bsxfun(@plus,(1:round(n*max_width))',round(n*(lower_min:step_sizes(1):lower_min+max_dropout_fraction))));
% % % % % % % % % % % % % % % % % % % % % X1 = X(1,:); X = bsxfun(@minus,X,X1);
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % opt_val = Inf;
% % % % % % % % % % % % % % % % % % % % % % for each interval width...
% % % % % % % % % % % % % % % % % % % % % for m = round(n*(max_width:-step_sizes(2):min_width))
% % % % % % % % % % % % % % % % % % % % %     % scale and bin the data in the intervals
% % % % % % % % % % % % % % % % % % % % %     nbins = round(3*log2(1+m/2));
% % % % % % % % % % % % % % % % % % % % %     H = bsxfun(@times,X(1:m,:),nbins./X(m,:));
% % % % % % % % % % % % % % % % % % % % %     logq = log(histc(H,[0:nbins-1,Inf]) + 0.01);
% % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % %     % for each shape value...
% % % % % % % % % % % % % % % % % % % % %     for b=1:length(beta)
% % % % % % % % % % % % % % % % % % % % %         bounds = zbounds{b};
% % % % % % % % % % % % % % % % % % % % %         % evaluate truncated generalized Gaussian pdf at bin centers
% % % % % % % % % % % % % % % % % % % % %         x = bounds(1)+(0.5:(nbins-0.5))/nbins*diff(bounds);
% % % % % % % % % % % % % % % % % % % % %         p = exp(-abs(x).^beta(b))*rescale(b); p=p'/sum(p);
% % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % %         % calc KL divergences
% % % % % % % % % % % % % % % % % % % % %         kl = sum(bsxfun(@times,p,bsxfun(@minus,log(p),logq(1:end-1,:)))) + log(m);
% % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % %         % update optimal parameters
% % % % % % % % % % % % % % % % % % % % %         [min_val,idx] = min(kl);
% % % % % % % % % % % % % % % % % % % % %         if min_val < opt_val
% % % % % % % % % % % % % % % % % % % % %             opt_val = min_val;
% % % % % % % % % % % % % % % % % % % % %             opt_beta = beta(b);
% % % % % % % % % % % % % % % % % % % % %             opt_bounds = bounds;
% % % % % % % % % % % % % % % % % % % % %             opt_lu = [X1(idx) X1(idx)+X(m,idx)];
% % % % % % % % % % % % % % % % % % % % %         end
% % % % % % % % % % % % % % % % % % % % %     end
% % % % % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % recover distribution parameters at optimum
% % % % % % % % % % % % % % % % % % % % % alpha = (opt_lu(2)-opt_lu(1))/diff(opt_bounds);
% % % % % % % % % % % % % % % % % % % % % mu = opt_lu(1)-opt_bounds(1)*alpha;
% % % % % % % % % % % % % % % % % % % % % beta = opt_beta;
% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % calculate the distribution's standard deviation from alpha and beta
% % % % % % % % % % % % % % % % % % % % % sig = sqrt((alpha^2)*gamma(3/beta)/gamma(1/beta));
% % % % % % % % % % % % % % % % % % % % % end


% % % function [signal,sample_mask] = clean_windows(signal,max_bad_channels,zthresholds,window_len,window_overlap,max_dropout_fraction,min_clean_fraction,truncate_quant,step_sizes,shape_range)
% % % if ~exist('max_bad_channels','var') || isempty(max_bad_channels) max_bad_channels = 0.2; end
% % % if ~exist('zthresholds','var') || isempty(zthresholds) zthresholds = [-3.5 5]; end
% % % if ~exist('window_len','var') || isempty(window_len) window_len = 1; end
% % % if ~exist('window_overlap','var') || isempty(window_overlap) window_overlap = 0.66; end
% % % if ~exist('max_dropout_fraction','var') || isempty(max_dropout_fraction) max_dropout_fraction = 0.1; end
% % % if ~exist('min_clean_fraction','var') || isempty(min_clean_fraction) min_clean_fraction = 0.25; end
% % % if ~exist('truncate_quant','var') || isempty(truncate_quant) truncate_quant = [0.022 0.6]; end
% % % if ~exist('step_sizes','var') || isempty(step_sizes) step_sizes = [0.01 0.01]; end
% % % if ~exist('shape_range','var') || isempty(shape_range) shape_range = 1.7:0.15:3.5; end
% % % if ~isempty(max_bad_channels) && max_bad_channels > 0 && max_bad_channels < 1 %#ok<*NODEF>
% % %     max_bad_channels = round(size(signal.data,1)*max_bad_channels); end
% % % 
% % % signal.data = double(signal.data);
% % % [C,S] = size(signal.data);
% % % N = window_len*signal.srate;
% % % wnd = 0:N-1;
% % % offsets = round(1:N*(1-window_overlap):S-N);
% % % 
% % % fprintf('Determining time window rejection thresholds...');
% % % % for each channel...
% % % for c = C:-1:1
% % %     % compute RMS amplitude for each window...
% % %     X = signal.data(c,:).^2;
% % %     X = sqrt(sum(X(bsxfun(@plus,offsets,wnd')))/N);
% % %     % robustly fit a distribution to the clean EEG part
% % %     [mu,sig] = fit_eeg_distribution(X, ...
% % %         min_clean_fraction, max_dropout_fraction, ...
% % %         truncate_quant, step_sizes,shape_range);
% % %     % calculate z scores relative to that
% % %     wz(c,:) = (X - mu)/sig;
% % % end
% % % disp('done.');
% % % 
% % % % sort z scores into quantiles
% % % swz = sort(wz);
% % % % determine which windows to remove
% % % remove_mask = false(1,size(swz,2));
% % % if max(zthresholds)>0
% % %     remove_mask(swz(end-max_bad_channels,:) > max(zthresholds)) = true; end
% % % if min(zthresholds)<0
% % %     remove_mask(swz(1+max_bad_channels,:) < min(zthresholds)) = true; end
% % % removed_windows = find(remove_mask);
% % % 
% % % % find indices of samples to remove
% % % removed_samples = repmat(offsets(removed_windows)',1,length(wnd))+repmat(wnd,length(removed_windows),1);
% % % % mask them out
% % % sample_mask = true(1,S); 
% % % sample_mask(removed_samples(:)) = false;
% % % fprintf('Keeping %.1f%% (%.0f seconds) of the data.\n',100*(mean(sample_mask)),nnz(sample_mask)/signal.srate);
% % % % determine intervals to retain
% % % % retain_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
% % % % retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
% % % signal.data = signal.data(:,sample_mask);
% % % end
% % % 
% % % 
% % % function signal = clean_asr(signal,cutoff,windowlen,stepsize,maxdims,ref_maxbadchannels,ref_tolerances,ref_wndlen,usegpu)
% % % if ~exist('cutoff','var') || isempty(cutoff) cutoff = 5; end
% % % if ~exist('windowlen','var') || isempty(windowlen) windowlen = max(0.5,1.5*signal.nbchan/signal.srate); end
% % % if ~exist('stepsize','var') || isempty(stepsize) stepsize = []; end
% % % if ~exist('maxdims','var') || isempty(maxdims) maxdims = 0.66; end
% % % if ~exist('ref_maxbadchannels','var') || isempty(ref_maxbadchannels) ref_maxbadchannels = 0.075; end
% % % if ~exist('ref_tolerances','var') || isempty(ref_tolerances) ref_tolerances = [-3.5 5.5]; end
% % % if ~exist('ref_wndlen','var') || isempty(ref_wndlen) ref_wndlen = 1; end
% % % if ~exist('usegpu','var') || isempty(usegpu) usegpu = false; end
% % % signal.data = double(signal.data);
% % % if isnumeric(ref_maxbadchannels) && isnumeric(ref_tolerances) && isnumeric(ref_wndlen)
% % %     disp('Finding a clean section of the data...');
% % %     try
% % %         ref_section = clean_windows(signal,ref_maxbadchannels,ref_tolerances,ref_wndlen); 
% % %     catch e
% % %         disp('An error occurred while trying to identify a subset of clean calibration data from the recording.');
% % %         disp('If this is because do not have EEGLAB loaded or no Statistics toolbox, you can generally');
% % %         disp('skip this step by passing in ''off'' as the ReferenceMaxBadChannels parameter.');
% % %         disp('Error details: ');
% % %         hlp_handleerror(e,1);
% % %         disp('Falling back to using the entire data for calibration.')
% % %         ref_section = signal;
% % %     end
% % % elseif strcmp(ref_maxbadchannels,'off') || strcmp(ref_tolerances,'off') || strcmp(ref_wndlen,'off')
% % %     disp('Using the entire data for calibration (reference parameters set to ''off'').')
% % %     ref_section = signal;
% % % elseif ischar(ref_maxbadchannels) && isvarname(ref_maxbadchannels)
% % %     disp('Using a user-supplied data set in the workspace.');
% % %     ref_section = evalin('base',ref_maxbadchannels);
% % % elseif all(isfield(ref_maxbadchannels,{'data','srate','chanlocs'}))
% % %     disp('Using a user-supplied clean section of data.');
% % %     ref_section = ref_maxbadchannels; 
% % % else
% % %     error('Unsupported value for argument ref_maxbadchannels.');
% % % end
% % % disp('Estimating calibration statistics; this may take a while...');
% % % if exist('hlp_diskcache','file')
% % %     state = hlp_diskcache('filterdesign',@asr_calibrate,ref_section.data,ref_section.srate,cutoff);
% % % else
% % %     state = asr_calibrate(ref_section.data,ref_section.srate,cutoff);
% % % end
% % % clear ref_section;
% % % if isempty(stepsize)
% % %     stepsize = floor(signal.srate*windowlen/2); end
% % % sig = [signal.data bsxfun(@minus,2*signal.data(:,end),signal.data(:,(end-1):-1:end-round(windowlen/2*signal.srate)))];
% % % [signal.data,state] = asr_process(sig,signal.srate,state,windowlen,windowlen/2,stepsize,maxdims,[],usegpu);
% % % signal.data(:,1:size(state.carry,2)) = [];
% % % end
% % % 
% % % 
% % % function state = asr_calibrate(X,srate,cutoff,blocksize,B,A,window_len,window_overlap,max_dropout_fraction,min_clean_fraction)
% % % [C,S] = size(X);
% % % if nargin < 3 || isempty(cutoff)
% % %     cutoff = 5; end
% % % if nargin < 4 || isempty(blocksize)
% % %     blocksize = 10; end
% % % blocksize = max(blocksize,ceil((C*C*S*8*3*2)/hlp_memfree));
% % % if nargin < 6 || isempty(A) || isempty(B)
% % %     try
% % %         [B,A] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
% % %     catch e %#ok<NASGU>
% % %         switch srate
% % %             case 100
% % %                 [B,A] = deal([0.9314233528641650 -1.0023683814963549 -0.4125359862018213  0.7631567476327510  0.4160430392910331 -0.6549131038692215 -0.0372583518046807  0.1916268458752655  0.0462411971592346],[1.0000000000000000 -0.4544220180303844 -1.0007038682936749  0.5374925521337940  0.4905013360991340 -0.4861062879351137 -0.1995986490699414  0.1830048420730026  0.0457678549234644]);
% % %             case 128
% % %                 [B,A] = deal([1.1027301639165037 -2.0025621813611867  0.8942119516481342  0.1549979524226999  0.0192366904488084  0.1782897770278735 -0.5280306696498717  0.2913540603407520 -0.0262209802526358],[1.0000000000000000 -1.1042042046423233 -0.3319558528606542  0.5802946221107337 -0.0010360013915635  0.0382167091925086 -0.2609928034425362  0.0298719057761086  0.0935044692959187]);
% % %             case 200
% % %                 [B,A] = deal([1.4489483325802353 -2.6692514764802775  2.0813970620731115 -0.9736678877049534  0.1054605060352928 -0.1889101692314626  0.6111331636592364 -0.3616483013075088  0.1834313060776763],[1.0000000000000000 -0.9913236099393967  0.3159563145469344 -0.0708347481677557 -0.0558793822071149 -0.2539619026478943  0.2473056615251193 -0.0420478437473110  0.0077455718334464]);
% % %             case 256
% % %                 [B,A] = deal([1.7587013141770287 -4.3267624394458641  5.7999880031015953 -6.2396625463547508  5.3768079046882207 -3.7938218893374835  2.1649108095226470 -0.8591392569863763  0.2569361125627988],[1.0000000000000000 -1.7008039639301735  1.9232830391058724 -2.0826929726929797  1.5982638742557307 -1.0735854183930011  0.5679719225652651 -0.1886181499768189  0.0572954115997261]);
% % %             case 300
% % %                 [B,A] = deal([1.9153920676433143  -5.7748421104926795   9.1864764859103936 -10.7350356619363630   9.6423672437729007  -6.6181939699544277   3.4219421494177711  -1.2622976569994351   0.2968423019363821],[1.0000000000000000 -2.3143703322055491  3.2222567327379434 -3.6030527704320621  2.9645154844073698 -1.8842615840684735  0.9222455868758080 -0.3103251703648485  0.0634586449896364]);
% % %             case 500
% % %                 [B,A] = deal([2.3133520086975823 -11.9471223009159130  29.1067166493384340 -43.7550171007238190  44.3385767452216370 -30.9965523846388000  14.6209883020737190  -4.2743412400311449   0.5982553583777899],[1.0000000000000000  -4.6893329084452580  10.5989986701080210 -14.9691518101365230  14.3320358399731820  -9.4924317069169977   4.2425899618982656  -1.1715600975178280   0.1538048427717476]);
% % %             case 512
% % %                 [B,A] = deal([2.3275475636130865 -12.2166478485960430  30.1632789058248850 -45.8009842020820410  46.7261263011068880 -32.7796858196767220  15.4623349612560630  -4.5019779685307473   0.6242733481676324],[1.0000000000000000  -4.7827378944258703  10.9780696236622980 -15.6795187888195360  15.1281978667576310 -10.0632079834518220   4.5014690636505614  -1.2394100873286753   0.1614727510688058]);
% % %             otherwise
% % %                 error('repair_bursts:NoYulewalk','The yulewalk() function was not found and there is no pre-computed spectral filter for your sampling rate. If you would like to use the default spectral filter please try to resample to one of the supported rates (100,128,200,256,300,500,512) or get the appropriate toobox license (you can also disable the spectral weighting feature or supply your own precalculated IIR filter coefficients).');
% % %         end
% % %     end
% % % end
% % % if nargin < 8 || isempty(window_len)
% % %     window_len = 0.5; end
% % % if nargin < 9 || isempty(window_overlap)
% % %     window_overlap = 0.66; end
% % % if nargin < 10 || isempty(max_dropout_fraction)
% % %     max_dropout_fraction = 0.1; end
% % % if nargin < 11 || isempty(min_clean_fraction)
% % %     min_clean_fraction = 0.25; end
% % % X(~isfinite(X(:))) = 0;
% % % [X,iirstate] = filter(B,A,double(X),[],2); X = X';
% % % if any(~isfinite(X(:)))
% % %     error('The IIR filter diverged on your data. Please try using either a more conservative filter or removing some bad sections/channels from the calibration data.'); end
% % % U = zeros(length(1:blocksize:S),C*C);
% % % for k=1:blocksize
% % %     range = min(S,k:blocksize:(S+k-1));
% % %     U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1)),size(U));
% % % end
% % % M = sqrtm(real(reshape(block_geometric_median(U/blocksize),C,C)));
% % % N = round(window_len*srate);
% % % fprintf('Determining per-component thresholds...');
% % % [V,D] = eig(M); %#ok<NASGU>
% % % X = abs(X*V);
% % % for c = C:-1:1
% % %     rms = X(:,c).^2;
% % %     rms = sqrt(sum(rms(bsxfun(@plus,round(1:N*(1-window_overlap):S-N),(0:N-1)')))/N);
% % %     [mu(c),sig(c)] = fit_eeg_distribution(rms,min_clean_fraction,max_dropout_fraction);
% % % end
% % % T = diag(mu + cutoff*sig)*V';
% % % disp('done.');
% % % state = struct('M',M,'T',T,'B',B,'A',A,'cov',[],'carry',[],'iir',iirstate,'last_R',[],'last_trivial',true);
% % % end
% % % 
% % % 
% % % function y = block_geometric_median(X,blocksize,varargin)
% % % if nargin < 2 || isempty(blocksize)
% % %     blocksize = 1; end
% % % 
% % % if blocksize > 1
% % %     [o,v] = size(X);       % #observations & #variables
% % %     r = mod(o,blocksize);  % #rest in last block
% % %     b = (o-r)/blocksize;   % #blocks
% % %     if r > 0
% % %         X = [reshape(sum(reshape(X(1:(o-r),:),blocksize,b*v)),b,v); sum(X((o-r+1):end,:))*(blocksize/r)];
% % %     else
% % %         X = reshape(sum(reshape(X,blocksize,b*v)),b,v);
% % %     end
% % % end
% % % 
% % % try
% % %     y = gather(geometric_median(gpuArray(X),varargin{:}))/blocksize;
% % % catch
% % %     y = geometric_median(X,varargin{:})/blocksize;
% % % end
% % % end
% % % 
% % % 
% % % function y = geometric_median(X,tol,y,max_iter)
% % % if ~exist('tol','var') || isempty(tol)
% % %     tol = 1.e-5; end
% % % if ~exist('y','var') || isempty(y)
% % %     y = median(X); end
% % % if ~exist('max_iter','var') || isempty(max_iter)
% % %     max_iter = 500; end
% % % 
% % % for i=1:max_iter
% % %     invnorms = 1./sqrt(sum(bsxfun(@minus,X,y).^2,2));
% % %     [y,oldy] = deal(sum(bsxfun(@times,X,invnorms)) / sum(invnorms),y);
% % %     if norm(y-oldy)/norm(y) < tol
% % %         break; end
% % % end
% % % end
% % % 
% % % 
% % % function result = hlp_memfree
% % % result = java.lang.management.ManagementFactory.getOperatingSystemMXBean().getFreePhysicalMemorySize();
% % % end
% % % 
% % % function [X,Zf] = moving_average(N,X,Zi)
% % % if nargin <= 2 || isempty(Zi)
% % %     Zi = zeros(size(X,1),N); end
% % % Y = [Zi X]; M = size(Y,2);
% % % I = [1:M-N; 1+N:M];
% % % S = [-ones(1,M-N); ones(1,M-N)]/N;
% % % X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% % % X = X(:,2:2:end);
% % % if nargout > 1
% % %     Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end
% % % end
% % % 
% % % 
% % % function [outdata,outstate] = asr_process(data,srate,state,windowlen,lookahead,stepsize,maxdims,maxmem,usegpu)
% % % if nargin < 4 || isempty(windowlen) 
% % %     windowlen = 0.5; end
% % % windowlen = max(windowlen,1.5*size(data,1)/srate);
% % % if nargin < 5 || isempty(lookahead)
% % %     lookahead = windowlen/2; end
% % % if nargin < 6 || isempty(stepsize)
% % %     stepsize = 32; end
% % % if nargin < 7 || isempty(maxdims)
% % %     maxdims = 0.66; end
% % % if nargin < 9 || isempty(usegpu)
% % %     usegpu = false; end
% % % if nargin < 8 || isempty(maxmem)
% % %     if usegpu
% % %         dev = gpuDevice(); maxmem = dev.FreeMemory/2^20; 
% % %     else
% % %         maxmem = hlp_memfree/(2^21);
% % %     end
% % % end
% % % if maxdims < 1
% % %     maxdims = round(size(data,1)*maxdims); end
% % % if isempty(data)
% % %     outdata = data; outstate = state; return; end
% % % [C,S] = size(data);
% % % N = round(windowlen*srate);
% % % P = round(lookahead*srate);
% % % [T,M,A,B] = deal(state.T,state.M,state.A,state.B);
% % % if isempty(state.carry)
% % %     state.carry = repmat(2*data(:,1),1,P) - data(:,1+mod(((P+1):-1:2)-1,S)); end
% % % data = [state.carry data];
% % % data(~isfinite(data(:))) = 0;
% % % splits = ceil((C*C*S*8*8 + C*C*8*S/stepsize + C*S*8*2 + S*8*5) / (maxmem*1024*1024 - C*C*P*8*3));
% % % if splits > 1
% % %     fprintf('Now cleaning data in %i blocks',splits); end
% % % for i=1:splits
% % %     range = 1+floor((i-1)*S/splits) : min(S,floor(i*S/splits));
% % %     if ~isempty(range)
% % %         [X,state.iir] = filter(B,A,double(data(:,range+P)),state.iir,2);
% % %         if usegpu && length(range) > 1000
% % %             try X = gpuArray(X); catch,end; end
% % %         [Xcov,state.cov] = moving_average(N,reshape(bsxfun(@times,reshape(X,1,C,[]),reshape(X,C,1,[])),C*C,[]),state.cov);
% % %         update_at = min(stepsize:stepsize:(size(Xcov,2)+stepsize-1),size(Xcov,2));
% % %         if isempty(state.last_R)
% % %             update_at = [1 update_at]; 
% % %             state.last_R = eye(C);
% % %         end
% % %         Xcov = reshape(Xcov(:,update_at),C,C,[]);
% % %         if usegpu
% % %             Xcov = gather(Xcov); end
% % %         last_n = 0;
% % %         for j=1:length(update_at)
% % %             [V,D] = eig(Xcov(:,:,j));
% % %             [D,order] = sort(reshape(diag(D),1,C)); V = V(:,order);
% % %             keep = D<sum((T*V).^2) | (1:C)<(C-maxdims);
% % %             trivial = all(keep);
% % %             if ~trivial
% % %                 R = real(M*pinv(bsxfun(@times,keep',V'*M))*V');
% % %             else
% % %                 R = eye(C);
% % %             end
% % %             n = update_at(j);
% % %             if ~trivial || ~state.last_trivial
% % %                 subrange = range((last_n+1):n);
% % %                 blend = (1-cos(pi*(1:(n-last_n))/(n-last_n)))/2;
% % %                 data(:,subrange) = bsxfun(@times,blend,R*data(:,subrange)) + bsxfun(@times,1-blend,state.last_R*data(:,subrange));
% % %             end
% % %             [last_n,state.last_R,state.last_trivial] = deal(n,R,trivial);
% % %         end
% % %     end
% % %     if splits > 1
% % %         fprintf('.'); end
% % % end
% % % if splits > 1
% % %     fprintf('\n'); end
% % % state.carry = [state.carry data(:,(end-P+1):end)];
% % % state.carry = state.carry(:,(end-P+1):end);
% % % outdata = data(:,1:(end-P));
% % % if usegpu
% % %     state.iir = gather(state.iir);
% % %     state.cov = gather(state.cov);
% % % end
% % % outstate = state;
% % % end
% % % 
% % % 
