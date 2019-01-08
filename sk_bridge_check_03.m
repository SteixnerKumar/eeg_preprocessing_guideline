% This function is created from the paper
% Identifying electrode bridging from electrical distance distributions:
% A survey of publicly-available EEG data using a new method
%
% ED (electrical difference)
% Formula:
% in the paper (temporal variance; Tenke and Kayser, 2001; Neuroscan Inc., 1993; 1995)
%
% author: Saurabh Kumar (s.kumar@uke.de)
%

% compare with ebridge
%  [EB_out ,ED_out] = eBridge(ALLEEG(2),{'1-EXG1','1-EXG2','1-EXG3','1-EXG4','1-EXG5','1-EXG6','1-EXG7','1-EXG8','2-EXG1','2-EXG2','2-EXG3','2-EXG4','2-EXG5','2-EXG6','2-EXG7','2-EXG8','Status'});

%% The parameters needed
epoch_length = 1; % (in seconds)
filename = 'tt_enemy_session_2_vp30101_vp30301.bdf'; % the file to be loaded
fs_new = 512; % new sampling frequency

%% reading data
data_both = ft_read_data(filename);
hdr_both = ft_read_header(filename);
events = ft_read_event(filename);

%% downsampling
% In the data
if ~(floor(hdr_both.Fs/fs_new) == hdr_both.Fs/fs_new)
    error('new sampling rate needs to be a proper divisor of original sampling rate');
end
data_both_temp = NaN(size(data_both,1),size(data_both,2)/(hdr_both.Fs/fs_new)); % initialization of a temporary variable
for loop_chan = 1:hdr_both.nChans % filter loop as channel-wise downsampling will be done
    data_both_temp(loop_chan,:) = downsample(data_both(loop_chan,:),(hdr_both.Fs/fs_new));
end
data_both = data_both_temp;
clear data_both_temp;
% In the events
for loop_markers = 1:length(strcmp('STATUS', {events.type}))
    if strcmp('STATUS', {events(loop_markers).type})
        events(loop_markers).sample = events(loop_markers).sample/(hdr_both.Fs/fs_new);
    end
end
% In the header
hdr_both.nSamples = hdr_both.nSamples/(hdr_both.Fs/fs_new);
hdr_both.Fs = hdr_both.Fs/(hdr_both.Fs/fs_new);
clear fs_new loop_chan loop_markers;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temporarily only taking the first 100000 samples (for proof of concept)
% data_both = data_both(:,1:102400); % 512*200
% data_both = data_both(:,1:51200); % 512*100
% data_both = data_both(:,1:25600); % 512*50
% data_both = data_both(:,1:5120); % 512*10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ignore the external channels and the trigger channel
data_both_new=nan(size(data_both,1),size(data_both,2));
for loop_chan = 1:size(data_both,1)
    if isempty(strfind(hdr_both.label{loop_chan,1},'EX')) && isempty(strfind(hdr_both.label{loop_chan,1},'Status'))
        data_both_new(loop_chan,:) = data_both(loop_chan,:);
    end
end
data_both_new(isnan(data_both_new(:,1)),:) = []; % deleting the nan (external electrodes)
data_both = data_both_new;
clear data_both_new



%% To get data_format : channels X epochs X sample-points -------:channels X channels X epochs
num_epochs = floor(size(data_both,2)/(epoch_length*hdr_both.Fs));
% channels X epochs X sample-points
data_both_new = nan(size(data_both,1),num_epochs,(epoch_length*hdr_both.Fs));
for loop_chan = 1:size(data_both,1)
    for loop_epoch = 1:num_epochs
        % epochs
        data_both_new(loop_chan,loop_epoch,1:(epoch_length*hdr_both.Fs)) = data_both(loop_chan,((epoch_length*hdr_both.Fs)*(loop_epoch-1)+1):((epoch_length*hdr_both.Fs)*loop_epoch));
    end
end
%

% calculate the ED
[jBarHandle,pb_fig] = sk_progressbar(1,((size(data_both_new,1))-1));
pause(0.5);
ED = nan(size(data_both_new,1),size(data_both_new,1),num_epochs);
for loop_chan_a = 1:((size(data_both_new,1))-1)
    for loop_chan_b = (loop_chan_a+1):(size(data_both_new,1))
        % potential difference between the channel
        diffE = squeeze(data_both_new((loop_chan_a),:,:) - data_both_new((loop_chan_b),:,:))';
        % electrical difference according to the formula (temporal variance; Tenke and Kayser, 2001; Neuroscan Inc., 1993; 1995)
        ED((loop_chan_a),(loop_chan_b),:) = var(diffE,1,1);
    end
    javaMethodEDT('setValue', jBarHandle, loop_chan_a); % update the progress bar
end
pause(0.5);
close(pb_fig);
clear diffE loop_chan_a loop_chan_b jBarHandle pb_fig loop_chan loop_epoch
%
%%
% load('ed')

%% Multiply by scale factor 100/median
ED = transpose(ED(~isnan(ED)));
ED = ED*(100/median(ED(:)));

%% plot vc
% min(ED(:)):max(ED(:))
nbins = (max(ED(:))-min(ED(:)))/0.25;
%
[counts, bins] = hist(ED(:),round(nbins));
plot(bins, counts); % get a line plot of the histogram
xlim([0 500]); ylim([0 100000])

%%
% save('ed')

% [counts, bins] = hist(ED(:),round(nbins));
% size(ED)

