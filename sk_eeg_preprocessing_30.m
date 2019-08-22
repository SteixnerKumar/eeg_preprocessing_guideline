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

%% add paths
% warning!!!!!!!!
% please change these accourding to your system setup
% eeglab
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b');
addpath(genpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/functions'));
removepath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/functions/octavefunc');
% eeglab plugins
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/plugins/');
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/plugins/amica1.5/');
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/plugins/clean_rawdata0.34/');
addpath(genpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/eeglab14_1_2b/plugins/tmullen-cleanline-696a7181b7d0/'));
% fieldtrip
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/fieldtrip-lite-20170706/fieldtrip-20170706/');
ft_defaults;
% SASICA toolbox
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/SASICA-master');
% sk (personal) utilities
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/sk_utilities');
% behavioral data folder (if needed)
addpath('/projects/crunchie/kumar/kumar_working_directory/matlab/work/sk_codes/EEG_data_analysis_01/behave_data');
% eeg data folder (if needed)
clear;
clc;

%% basic required parameters
tt.version = 'enemy'; % 'single', 'enemy', 'friend'
tt.session = 3;%2; % session number
tt.sub_A = 30121;%30124; %
tt.sub_B = 30321;%30324; %
%
wanted.save = 0;
wanted.photo_diode_trigger_adjustment = 1; % '1' for yes, '0' for no
wanted.downsampling = 1; % '1' for yes, '0' for no
wanted.fs_new = 256;%512; % new sampling frequency
wanted.bridging_check = 1; % '1' for yes, '0' for no
wanted.filtering= 1; % '1' for yes, '0' for no (Please filter, else in ASR there is no clean section of data to choose from, hence error)
wanted.rereferencing = 1; % '1' for yes, '0' for no
wanted.clean_data = 1; % '1' for yes, '0' for no
wanted.ica = 1; % '1' for yes, '0' for no
wanted.ica_algorithm = 'amica'; % wanted.ica_algorithm = 'runica' or 'amica'
wanted.trials = 1;
wanted.plots = 0; % 1 to check the filter
plots_see_chan = 64;%15; % the channel number to see


%% Logging
diary(strcat('log_tt',strcat('_',datestr(now,'YYYYmmDD'),'_',datestr(now,'HHMMSS'),'_'),...
    tt.version,'_session_',num2str(tt.session),'_vp',num2str(tt.sub_A),'_vp',num2str(tt.sub_B),'.txt'));
% print the default parameters for the log
clc;
disp('The file parameters : ');disp(tt);fprintf('\n');
disp('The code parameters -  ''1'' for yes, ''0'' for no : ');disp(wanted);fprintf('\n');
try
    tic;
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
    
    %% tiggers timing adjustment according to the photo diodes
    if wanted.photo_diode_trigger_adjustment
        fprintf('Adjusting the data according to the delays of the screens ...\n');
        % the diode_avg.default structure is important as it decides the
        % triggers that are adjusted, please make sure the names are exactly as in the behaviour files in ttsk.trig structure
        diode_avg.default.start_prediction = 28;
        diode_avg.default.start_choice = 27;
        diode_avg.default.start_evidence_other_player = 143;
        diode_avg.default.start_evidence_own = 185;
        %         diode_avg.default.end_choice = 10;
        % to get the trigger values
        events_trigger = events(strcmp({events.type}, 'STATUS'));
        byte.second = nan(1,size(events_trigger,2));
        for loop_triggers = 1:size(events_trigger,2)
            in_binary = dec2bin(events_trigger(loop_triggers).value,24)-'0';
            binary_usb = in_binary(9:16);
            byte.second(loop_triggers) = bin2dec(num2str(binary_usb));
        end
        % now here the avg of all the triggers
        % find all the positions of 210,220,230,240 in byte.second
        trig_cat = fieldnames(diode_avg.default);
        for loop_trig_cat = 1:numel(fieldnames(diode_avg.default))
            eval(string(strcat('trig_diode_location_def.',trig_cat(loop_trig_cat),' = find(byte.second==behave_data.A.ttsk.trig.',trig_cat(loop_trig_cat),') +1;')));
            eval(string(strcat('trig_diode_location.',trig_cat(loop_trig_cat),' = trig_diode_location_def.',trig_cat(loop_trig_cat),'(byte.second(trig_diode_location_def.',trig_cat(loop_trig_cat),')== 0);')));
            eval(string(strcat('trig_diode_location_def.',trig_cat(loop_trig_cat),' = trig_diode_location_def.',trig_cat(loop_trig_cat),'-1;')));
            eval(string(strcat('diode_sum.',trig_cat(loop_trig_cat),' = 0;')));
            eval(string(strcat('for loop_trig_diode_location = 1:length(trig_diode_location.',trig_cat(loop_trig_cat),');diode_sum.',trig_cat(loop_trig_cat),' = diode_sum.',trig_cat(loop_trig_cat)...
                ,' + (events_trigger(trig_diode_location.',trig_cat(loop_trig_cat),'(loop_trig_diode_location)).sample - events_trigger(trig_diode_location.',trig_cat(loop_trig_cat)...
                ,'(loop_trig_diode_location)-1).sample);end')));
            eval(string(strcat('diode_avg.',trig_cat(loop_trig_cat),' = floor(diode_sum.',trig_cat(loop_trig_cat),'/length(trig_diode_location.',trig_cat(loop_trig_cat),'));')));
        end
        % adjust the delay if trigger found
        for loop_trig_cat = 1:numel(fieldnames(diode_avg.default))
            eval(string(strcat('if ~isnan(diode_avg.',trig_cat(loop_trig_cat),');for loop_trig_diode_location = 1: length(trig_diode_location_def.',trig_cat(loop_trig_cat),');events_trigger(trig_diode_location_def.',...
                trig_cat(loop_trig_cat),'(loop_trig_diode_location)).sample = events_trigger(trig_diode_location_def.',trig_cat(loop_trig_cat),'(loop_trig_diode_location)).sample',...
                ' + diode_avg.',trig_cat(loop_trig_cat),';end;','fprintf(''adjusted the delay of : %s screen by %d points\n'',string(trig_cat(loop_trig_cat)),diode_avg.',trig_cat(loop_trig_cat),');',...
                'else;fprintf(''did not find the delay of : %s screen as no trigger found, therefore adjusting by average default value of: %d points\n'',string(trig_cat(loop_trig_cat)),diode_avg.default.',trig_cat(loop_trig_cat),');',...
                'for loop_trig_diode_location = 1: length(trig_diode_location_def.',trig_cat(loop_trig_cat),');events_trigger(trig_diode_location_def.',trig_cat(loop_trig_cat)...
                ,'(loop_trig_diode_location)).sample = events_trigger(trig_diode_location_def.',trig_cat(loop_trig_cat),'(loop_trig_diode_location)).sample'...
                ,' + diode_avg.default.',trig_cat(loop_trig_cat),';end;end;')));
        end
        %
        % Now important to put the event triggers back into the events for use
        % by other sections of this code.
        events(strcmp({events.type}, 'STATUS')) = events_trigger;
        %
        clear loop_trig_diode_location loop_triggers trig_diode_location_def trig_diode_location diode_sum diode_avg
        clear diode_sum trig_diode_location
        clear in_binary loop_try binary_lsb binary_usb binary_third
        fprintf('Adjusting the data according to the delays of the screens ... Done.\n');
    else
        fprintf('Adjusting the data according to the delays of the screens not Done.\n ');
    end
    
    %% downsampling
    % maybe use the anti-aliasing low pass filter before the downsampling process (yet to be imnplemented)
    % downsample is needed to save the computation power and time --
    % needs to be done on the data, in the header and in the events.
    
    % In the data
    if wanted.downsampling
        fprintf('Downsampling the data by a factor of %d ...\n ',hdr.both.Fs/wanted.fs_new);
        if ~(floor(hdr.both.Fs/wanted.fs_new) == hdr.both.Fs/wanted.fs_new)
            error('new sampling rate needs to be a proper divisor of original sampling rate');
        end
        data_both_temp = NaN(size(data.both,1),ceil(size(data.both,2)/(hdr.both.Fs/wanted.fs_new))); % initialization of a temporary variable
        fprintf('individual channel-wise out of %d\n',hdr.both.nChans);
        for loop_chan = 1:hdr.both.nChans % filter loop as channel-wise downsampling will be done
            data_both_temp(loop_chan,:) = downsample(data.both(loop_chan,:),(hdr.both.Fs/wanted.fs_new));
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
                events(loop_markers).sample = floor(events(loop_markers).sample/(hdr.both.Fs/wanted.fs_new));
            end
        end
        % In the header
        hdr.both.nSamples = ceil(hdr.both.nSamples/(hdr.both.Fs/wanted.fs_new));
        hdr.both.Fs = hdr.both.Fs/(hdr.both.Fs/wanted.fs_new);
        clear wanted.fs_new loop_chan loop_markers;
        fprintf('Downsampling of the data done.\n ');
    elseif wanted.downsampling==0
        fprintf('Downsampling not Done.\n ');
    end
    
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
    if wanted.bridging_check
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
    else
        fprintf('Bridging check not done.\n');
    end
    
    
    %% Filtering (highpass, then lowpass)
    % implementation of the line noise filter
    % cleanline(EEG, 'LineFrequencies',[50 100]);
    % [EEG.A, Sorig, Sclean, f, amps, freqs, g] = cleanline('EEG',EEG.A, 'LineFrequencies',[50 100]);
    % pre filter plot
    if wanted.filtering
        fprintf('Filtering the line noise out for both participants ...\n');
        if wanted.plots
            figure;
            subplot(3,2,1);plot(data.A(plots_see_chan,:));title('The unfiltered data time');
            [XFreqRange, YAmplitude] = sk_dofft(data.A(plots_see_chan,:), hdr.A.Fs, 4);
            subplot(3,2,2);plot(XFreqRange, YAmplitude);axis([5 150 -0.2 1.5]);title('The unfiltered data freq');
            clear XFreqRange YAmplitude;
        end
        EEG.A.data = data.A;
        EEG.A.icawinv = [];
        EEG.A.srate = hdr.both.Fs;
        EEG.A.nbchan = size(data.A,1);
        EEG.A.trials = 1; % no trials
        EEG.A.pnts = size(data.A,2);
        EEG.B.data = data.B;
        EEG.B.icawinv = [];
        EEG.B.srate = hdr.both.Fs;
        EEG.B.nbchan = size(data.B,1);
        EEG.B.trials = 1; % no trials
        EEG.B.pnts = size(data.B,2);
        % [EEG.A, ~, ~, ~, ~, , ~] = cleanline('EEG',EEG.A, 'LineFrequencies',[50 100]);
        [EEG.A, ~, ~, ~, ~, ~, ~] = cleanline('EEG',EEG.A, 'LineFrequencies',[50 100],'ScanForLines',1,'LineAlpha',1,'Bandwidth',1);
        [EEG.B, ~, ~, ~, ~, ~, ~] = cleanline('EEG',EEG.B, 'LineFrequencies',[50 100],'ScanForLines',1,'LineAlpha',1,'Bandwidth',1);
        data.Af = EEG.A.data;
        data.Bf = EEG.B.data;
        clear EEG
        % post line filter plot
        if wanted.plots
            subplot(3,2,3);plot(data.Af(plots_see_chan,:));title('The filtered data time');
            [XFreqRange, YAmplitude] = sk_dofft(data.Af(plots_see_chan,:), hdr.A.Fs, 4);
            subplot(3,2,4);plot(XFreqRange, YAmplitude);axis([5 150 -0.2 1.5]);title('The filtered data freq');
            clear XFreqRange YAmplitude;
        end
        fprintf('Filtering the line noise out for both participants ... Done.\n');
        %
        %
        fprintf('Using low-pass filtering for both participants ...\n');
        filter.type = 'but'; % the filter type ('but' is for butterworth)
        filter.dir = 'twopass'; % the filter direction ('twopass' is default for both direction)
        filter.order = 4; % filter order number
        filter.lowpass = 90; % in Hz
        filter.highpass = 1; % in Hz
        fprintf('Using filter: %s, with direction: %s, order: %d, in range: %d - %d Hz.\n',filter.type,filter.dir,filter.order,filter.highpass,filter.lowpass);
        [data.Af] = ft_preproc_bandpassfilter(data.Af, hdr.A.Fs, [filter.highpass filter.lowpass], filter.order, filter.type, filter.dir);
        [data.Bf] = ft_preproc_bandpassfilter(data.Bf, hdr.B.Fs, [filter.highpass filter.lowpass], filter.order, filter.type, filter.dir);
        fprintf('Using low-pass filtering for both participants ... Done.\n');
        % post all filter plot
        if wanted.plots
            subplot(3,2,5);plot(data.Af(plots_see_chan,:));title('The filtered data time');
            [XFreqRange, YAmplitude] = sk_dofft(data.Af(plots_see_chan,:), hdr.A.Fs, 4);
            subplot(3,2,6);plot(XFreqRange, YAmplitude);axis([5 150 -0.2 1.5]);title('The filtered data freq');
            clear XFreqRange YAmplitude;
        end
    else
        fprintf('Data filtering not done.\n');
        data.Af = data.A;
        data.Bf = data.B;
    end
    
    %% reject bad channels
    % using cleanraw from the eeglab
    % a) flatline channels
    % b) noisy channels
    % c) short-time bursts
    % d) incompletely repaird segments from the data
    
    A.channels_to_remove = zeros(hdr.A.nChans,1);
    B.channels_to_remove = zeros(hdr.B.nChans,1);
    if wanted.clean_data
        
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
        
    else
        fprintf('Data cleaning not done.\n');
        data.Af_cut_asr_repaired = data.Af;
        data.Bf_cut_asr_repaired = data.Bf;
    end
    
    
    %% interpolating the missing channels
    % Interpolate channels. (not ideal before ICA, but better for re-referencing)
    % get EOG channels
    A.channels_eog =  zeros(hdr.A.nChans,1);
    B.channels_eog =  zeros(hdr.B.nChans,1);
    A.channels_eog(contains(hdr.A.label,'EX')) = 1;
    B.channels_eog(contains(hdr.B.label,'EX')) = 1;
    %
    fprintf('---Interpolating missing channels for the re-referencing part ....\n');
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
    fprintf('---Interpolating missing channels for the re-referencing part ... Done.\n');
    
    %% Re-referencing
    if wanted.rereferencing
        fprintf('re-referenciong the data to the common average reference (CAR)...\n');
        mean_data_one = mean(A.EEG_interp.data,1);
        mean_data_two = mean(B.EEG_interp.data,1);
        A.EEG_interp.data_car = A.EEG_interp.data - repmat(mean_data_one,size(A.EEG_interp.data,1),1);
        B.EEG_interp.data_car = B.EEG_interp.data - repmat(mean_data_two,size(B.EEG_interp.data,1),1);
        data.Af_cut_asr_repaired_interp_car = A.EEG_interp.data_car;
        data.Bf_cut_asr_repaired_interp_car = B.EEG_interp.data_car;
        clear mean_data_one mean_data_two
        fprintf('re-referenciong the data to the common average reference (CAR)... Done.\n');
    else
        fprintf('Re-referenciong not done.\n');
        data.Af_cut_asr_repaired_interp_car = data.Af_cut_asr_repaired_interp;
        data.Bf_cut_asr_repaired_interp_car = data.Bf_cut_asr_repaired_interp;
        A.EEG_interp.data_car = A.EEG_interp.data;
        B.EEG_interp.data_car = B.EEG_interp.data;
    end
    
    
    %% remove the interpolated channels before ICA
    data.Af_cut_asr_repaired_interp_car_cut = A.EEG_interp.data_car(~A.channels_to_remove(~A.channels_eog),:);
    data.Bf_cut_asr_repaired_interp_car_cut = B.EEG_interp.data_car(~B.channels_to_remove(~B.channels_eog),:);
    
    %% now the ICA
    if wanted.ica
        %
        if strcmp(wanted.ica_algorithm,'amica')
            fprintf('performing the ICA using the AMICA algorithm...\n');
            [A.ica.icaweights,A.ica.icasphere,A.ica.mods] = runamica15(data.Af_cut_asr_repaired_interp_car_cut);
            [B.ica.icaweights,B.ica.icasphere,B.ica.mods] = runamica15(data.Bf_cut_asr_repaired_interp_car_cut);
            fprintf('performing the ICA using the AMICA algorithm... Done.\n');
        elseif strcmp(wanted.ica_algorithm,'runica')
            fprintf('performing the ICA using the runica algorithm...\n');
            [A.ica.icaweights,A.ica.icasphere,A.ica.icameanvar,A.ica.icabias,A.ica.icasigns,A.ica.icalrates,A.ica.icadata,A.ica.icay] = runica(data.Af_cut_asr_repaired_interp_car_cut); % train using defaults
            [B.ica.icaweights,B.ica.icasphere,B.ica.icameanvar,B.ica.icabias,B.ica.icasigns,B.ica.icalrates,B.ica.icadata,B.ica.icay] = runica(data.Bf_cut_asr_repaired_interp_car_cut); % train using defaults
            fprintf('performing the ICA using the runica algorithm... Done.\n');
        else
            error('Please choose either amica or runica algorithm in the ICA section');
        end
        
        
        %% Auto-detection of bad ICA components
        fprintf('auto detecting the bad ICA components...\n');
        % calculate the inverse of the ica weights
        A.ica.icawinv = pinv(A.ica.icaweights*A.ica.icasphere);
        B.ica.icawinv = pinv(B.ica.icaweights*B.ica.icasphere);
        % other parameters
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
        if strcmp(wanted.ica_algorithm,'amica')
            A.ica_clean.icaact = A.ica.icaweights*A.ica.icasphere*A.ica.data;  % Matrix multiplication
            B.ica_clean.icaact = B.ica.icaweights*B.ica.icasphere*B.ica.data;  % Matrix multiplication
        elseif strcmp(wanted.ica_algorithm,'runica')
            A.ica_clean.icaact = A.ica.icaweights*A.ica.icasphere*A.ica.icadata;  % Matrix multiplication
            B.ica_clean.icaact = B.ica.icaweights*B.ica.icasphere*B.ica.icadata;  % Matrix multiplication
        else
            error('Please choose either amica or runica algorithm in the ICA section');
        end
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
        
    else
        fprintf('ICA not done.\n');
    end
    
    %% cut the data here into trials/epochs
    if wanted.trials && wanted.ica
        % for participant A
        sktt = [];
        sktt.Fs = hdr.both.Fs;
        sktt.events = events;
        sktt.data = A.ica_clean.data;
        sktt.flag_bc_trial = 1;
        sktt.trig = behave_data.A.ttsk.trig;
        [temp_data, ~] = sk_trials_creation(sktt);
        A.ica_clean_trials = temp_data.trials;
        clear ttsk temp_data
        %
        % for participant B
        sktt = [];
        sktt.Fs = hdr.both.Fs;
        sktt.events = events;
        sktt.data = B.ica_clean.data;
        sktt.flag_bc_trial = 1;
        % temporarily interchange the A and B trigger numbers for cutting the
        % data according to the participant B
        temp_trig_B = behave_data.B.ttsk.trig;
        temp_trig_B.response_prediction_A = behave_data.B.ttsk.trig.response_prediction_B;
        temp_trig_B.response_prediction_B = behave_data.B.ttsk.trig.response_prediction_A;
        temp_trig_B.response_choice_A = behave_data.B.ttsk.trig.response_prediction_B;
        temp_trig_B.response_choice_B = behave_data.B.ttsk.trig.response_prediction_A;
        sktt.trig = temp_trig_B;
        [temp_data, ~] = sk_trials_creation(sktt);
        B.ica_clean_trials = temp_data.trials;
        clear ttsk temp_data temp_trig_B
    else
        fprintf('data not cut into trials.\n');
    end
    %
    %% Save the relevent data at this point
    %
    if wanted.save
        fprintf('------ saving the ICA cleaned session ------ and ------ ICA cleaned data in trials. ------ \n');
        if wanted.ica
            temp_save_A.data = A.ica_clean;
            temp_save_B.data = B.ica_clean;
        end
        if wanted.trials
            temp_save_A.trials = A.ica_clean_trials;
            temp_save_B.trials = B.ica_clean_trials;
        end
        temp_save_A.original = data.A;
        temp_save_B.original = data.B;
        temp_save_A.hdr = hdr.A;
        temp_save_B.hdr = hdr.B;
        temp_save_A.events = events;
        temp_save_B.events = events;
        temp_save_A.settings = wanted;
        temp_save_B.settings = wanted;
        temp_save_A.behaviour = behave_data.A;
        temp_save_B.behaviour = behave_data.B;
        temp_save_A.elec = elec;
        temp_save_B.elec = elec;
        temp_save_A.bad_chans = bad_chans_A;
        temp_save_B.bad_chans = bad_chans_B;
        temp_save_A.name = strcat('tt_preprocessed_',strcat('_',datestr(now,'YYYYmmDD'),'_',datestr(now,'HHMMSS'),'_'),...
            tt.version,'_session_',num2str(tt.session),'_vp',num2str(tt.sub_A),'.mat');
        temp_save_B.name = strcat('tt_preprocessed_',strcat('_',datestr(now,'YYYYmmDD'),'_',datestr(now,'HHMMSS'),'_'),...
            tt.version,'_session_',num2str(tt.session),'_vp',num2str(tt.sub_B),'.mat');
        save(temp_save_A.name,'temp_save_A', '-v7.3');
        save(temp_save_B.name,'temp_save_B', '-v7.3');
        fprintf('------ saving the ICA cleaned session ------ and ------ ICA cleaned data in trials. ------ Done.\n');
    else
        fprintf('data not saved.\n');
    end
    %
    %%
    fprintf('\n');
    toc;
    diary('off');
catch exception
    fprintf('\nlogging stopped because of error.\n');
    toc;
    diary('off');
    rethrow(exception);
end

%%
