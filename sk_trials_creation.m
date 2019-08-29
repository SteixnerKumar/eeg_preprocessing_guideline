function [data, trial] = sk_trials_creation(sktt)
%
% ***************************************************************
% To cut the data into trials/epochs
% (the structure will be channals X trial-length X trial-number)
% ***************************************************************
% YET to IMPLEMENT THE FORCED LISTEN prediction and choice trial correction
% Usage:
% function [data, trial] = sk_trials_creation(sktt)
%
% The data input to this function needs to be arranged in the perspective
% of the main participant being A, the participant who's data is cut. The
% partner participant is called as B. This is an important segregation to
% be noted else the data may be falsly cut.
% So we mainly need to look at the input trigger values because it
% specifically has the A and the B part. Based on this information the
% entire algorithm is based. (the main decision will be to either
% temporarily flip the triggers for A and B for this function or not !)
%
% The trials are taken as :
%       1.) between the start of the prediction phase to the end
%       2.) between the start of the choice phase to the end
%       3.) between the start of the own_evidence phase to the end
%       4.) between the start of the partner_evidence phase to the end
%       5.) between the start of the prediction phase to the end of own response
%       6.) between the start of the choice phase to the end of own response
%       7.) between the start of the prediction phase to the end of partner response
%       8.) between the start of the choice phase to the end of partner response
%
% _______________________________________________________________
% OUTPUT
% data : the cut data that was put in in the defined structure (the
%           structure will be channals X trial-length X trial-number)
% trial : the exact cut points for each trials
% _______________________________________________________________
%
% INPUT
% sktt. --
% Fs : The sampling frequency.
% data : The main data that needs to be cut into trials.
% events : The trigger data points instances, as to know at which points
%                       the data need to be cut. In the function its mainly used as
%                       'events(strcmp({sktt.events.type}, 'STATUS'))', meaning the type STATUS
%                       is used to identify the trigger rows. In short the
%                       events is a struct array with atleast the 'type'
%                       and 'value' fields having the type 'STATUS'.
% trig : the trigger values that were used. This has all the trigger
%                       values that were used. In this function the trigger values used are:
%                               01.) start_experiment
%                               02.) end_experiment
%                               03.) start_prediction
%                               04.) end_prediction
%                               05.) start_choice
%                               06.) end_choice
%                               07.) response_prediction_A
%                               08.) response_choice_A
%                               09.) response_prediction_B
%                               10.) response_choice_B
%                                       Where, A is the main participant
%                                       who's data is cut, B is the partner
%                                       participant.
% flag_bc_trial : '1', '0', as to whether the baseline correction needs to
%                       be performed (only substracting the mean of the while trial is
%                       implemented in case if '1' is selected)
%
%
% creater: Saurabh kumar
%

% for one participant only (based on the partner ofcourse)
fprintf('Cutting the data into epochs/trials ...\n');
if sktt.flag_bc_trial
    temp_baseline_correction_time = 200; % time in milliseconds
    temp_baseline_correction_time = ceil((temp_baseline_correction_time/1000)*sktt.Fs); % in points
end
% to get the trigger values
sktt.events_trigger = sktt.events(strcmp({sktt.events.type}, 'STATUS'));
byte.second = nan(1,size(sktt.events_trigger,2));
for loop_triggers = 1:size(sktt.events_trigger,2)
    in_binary = dec2bin(sktt.events_trigger(loop_triggers).value,24)-'0';
    binary_usb = in_binary(9:16);
    byte.second(loop_triggers) = bin2dec(num2str(binary_usb));
end

% trial sections
trial.num = length(find(byte.second == sktt.trig.start_evidence_own));
trial.prediction = [sktt.trig.start_prediction sktt.trig.end_prediction];
trial.choice = [sktt.trig.start_choice sktt.trig.end_choice];
trial.evidence_other = [sktt.trig.start_evidence_other_player sktt.trig.end_evidence_other_player];
trial.evidence_own = [sktt.trig.start_evidence_own sktt.trig.end_evidence_own];
trial.response_prediction_A = [sktt.trig.start_prediction sktt.trig.response_prediction_A];
trial.response_choice_A = [sktt.trig.start_choice sktt.trig.response_choice_A];
trial.response_prediction_B = [sktt.trig.start_prediction sktt.trig.response_prediction_B];
trial.response_choice_B = [sktt.trig.start_choice sktt.trig.response_choice_B];
trial.num_response_prediction_A = length(find(byte.second == sktt.trig.response_prediction_A));
trial.num_response_prediction_B = length(find(byte.second == sktt.trig.response_prediction_B));
trial.num_response_choice_A = length(find(byte.second == sktt.trig.response_choice_A));
trial.num_response_choice_B = length(find(byte.second == sktt.trig.response_choice_B));
% get the points in the data
trial.points.prediction = [[sktt.events_trigger(byte.second == trial.prediction(1)).sample];[sktt.events_trigger(byte.second == trial.prediction(2)).sample]]';
trial.points.choice = [[sktt.events_trigger(byte.second == trial.choice(1)).sample];[sktt.events_trigger(byte.second == trial.choice(2)).sample]]';
trial.points.evidence_other = [[sktt.events_trigger(byte.second == trial.evidence_other(1)).sample];[sktt.events_trigger(byte.second == trial.evidence_other(2)).sample]]';
trial.points.evidence_own = [[sktt.events_trigger(byte.second == trial.evidence_own(1)).sample];[sktt.events_trigger(byte.second == trial.evidence_own(2)).sample]]';
if length(find(byte.second == trial.response_prediction_A(2))) ~= length(find(byte.second == trial.response_prediction_A(1))) % forced predictions
    temp_one = find(byte.second == trial.response_prediction_A(1));
    temp_two = find(byte.second == trial.response_prediction_A(2));
    count_three = 1;
    temp_three =nan(1,length(temp_two));
    for loop_one_two_outer = 1:length(temp_two)
        temp_temp = nan(1,length(temp_one));
        for loop_one_two = 1:length(temp_one)
            temp_temp(loop_one_two) = abs(temp_one(loop_one_two) -  temp_two(loop_one_two_outer));
        end
        temp_temp_loc = temp_temp == min(temp_temp);
        temp_three(count_three) = temp_one(temp_temp_loc);temp_two(loop_one_two_outer);
        count_three = count_three +1 ;
    end
    clear count_three loop_one_two_outer loop_one_two temp_temp temp_temp_loc temp_one
    temp_one = zeros(1,length(byte.second));
    temp_one(temp_three) = 1;
    trial.points.response_prediction_A = [[sktt.events_trigger(logical(temp_one)).sample];[sktt.events_trigger(byte.second == trial.response_prediction_A(2)).sample]]';
    clear temp_one temp_two temp_three
else
    trial.points.response_prediction_A = [[sktt.events_trigger(byte.second == trial.response_prediction_A(1)).sample];[sktt.events_trigger(byte.second == trial.response_prediction_A(2)).sample]]';
end
if length(find(byte.second == trial.response_choice_A(2))) ~= length(find(byte.second == trial.response_choice_A(1))) % forced listens
    temp_one = find(byte.second == trial.response_choice_A(1));
    temp_two = find(byte.second == trial.response_choice_A(2));
    count_three = 1;
    temp_three =nan(1,length(temp_two));
    for loop_one_two_outer = 1:length(temp_two)
        temp_temp = nan(1,length(temp_one));
        for loop_one_two = 1:length(temp_one)
            temp_temp(loop_one_two) = abs(temp_one(loop_one_two) -  temp_two(loop_one_two_outer));
        end
        temp_temp_loc = temp_temp == min(temp_temp);
        temp_three(count_three) = temp_one(temp_temp_loc);temp_two(loop_one_two_outer);
        count_three = count_three +1 ;
    end
    clear count_three loop_one_two_outer loop_one_two temp_temp temp_temp_loc temp_one
    temp_one = zeros(1,length(byte.second));
    temp_one(temp_three) = 1;
    trial.points.response_choice_A = [[sktt.events_trigger(logical(temp_one)).sample];[sktt.events_trigger(byte.second == trial.response_choice_A(2)).sample]]';
    clear temp_one temp_two temp_three
else
    trial.points.response_choice_A = [[sktt.events_trigger(byte.second == trial.response_choice_A(1)).sample];[sktt.events_trigger(byte.second == trial.response_choice_A(2)).sample]]';
end
if length(find(byte.second == trial.response_prediction_B(2))) ~= length(find(byte.second == trial.response_prediction_B(1))) % forced predictions
    temp_one = find(byte.second == trial.response_prediction_B(1));
    temp_two = find(byte.second == trial.response_prediction_B(2));
    count_three = 1;
    temp_three =nan(1,length(temp_two));
    for loop_one_two_outer = 1:length(temp_two)
        temp_temp = nan(1,length(temp_one));
        for loop_one_two = 1:length(temp_one)
            temp_temp(loop_one_two) = abs(temp_one(loop_one_two) -  temp_two(loop_one_two_outer));
        end
        temp_temp_loc = temp_temp == min(temp_temp);
        temp_three(count_three) = temp_one(temp_temp_loc);temp_two(loop_one_two_outer);
        count_three = count_three +1 ;
    end
    clear count_three loop_one_two_outer loop_one_two temp_temp temp_temp_loc temp_one
    temp_one = zeros(1,length(byte.second));
    temp_one(temp_three) = 1;
    trial.points.response_prediction_B = [[sktt.events_trigger(logical(temp_one)).sample];[sktt.events_trigger(byte.second == trial.response_prediction_B(2)).sample]]';
    clear temp_one temp_two temp_three
else
    trial.points.response_prediction_B = [[sktt.events_trigger(byte.second == trial.response_prediction_B(1)).sample];[sktt.events_trigger(byte.second == trial.response_prediction_B(2)).sample]]';
end
if length(find(byte.second == trial.response_choice_B(2))) ~= length(find(byte.second == trial.response_choice_B(1))) % forced listens
    temp_one = find(byte.second == trial.response_choice_B(1));
    temp_two = find(byte.second == trial.response_choice_B(2));
    count_three = 1;
    temp_three =nan(1,length(temp_two));
    for loop_one_two_outer = 1:length(temp_two)
        temp_temp = nan(1,length(temp_one));
        for loop_one_two = 1:length(temp_one)
            temp_temp(loop_one_two) = abs(temp_one(loop_one_two) -  temp_two(loop_one_two_outer));
        end
        temp_temp_loc = temp_temp == min(temp_temp);
        temp_three(count_three) = temp_one(temp_temp_loc);temp_two(loop_one_two_outer);
        count_three = count_three +1 ;
    end
    clear count_three loop_one_two_outer loop_one_two temp_temp temp_temp_loc temp_one
    temp_one = zeros(1,length(byte.second));
    temp_one(temp_three) = 1;
    trial.points.response_choice_B = [[sktt.events_trigger(logical(temp_one)).sample];[sktt.events_trigger(byte.second == trial.response_choice_B(2)).sample]]';
    clear temp_one temp_two temp_three
else
    trial.points.response_choice_B = [[sktt.events_trigger(byte.second == trial.response_choice_B(1)).sample];[sktt.events_trigger(byte.second == trial.response_choice_B(2)).sample]]';
end
% cutting the data
% rowsXcolumnXtrials
data.trials.prediction.data = nan(size(sktt.data,1),max(trial.points.prediction(:,2)-trial.points.prediction(:,1))+1,trial.num);
data.trials.choice.data = nan(size(sktt.data,1),max(trial.points.choice(:,2)-trial.points.choice(:,1))+1,trial.num);
data.trials.evidence_other.data = nan(size(sktt.data,1),max(trial.points.evidence_other(:,2)-trial.points.evidence_other(:,1))+1,trial.num);
data.trials.evidence_own.data = nan(size(sktt.data,1),max(trial.points.evidence_own(:,2)-trial.points.evidence_own(:,1))+1,trial.num);
data.trials.response_prediction_A.data = nan(size(sktt.data,1),max(trial.points.response_prediction_A(:,2)-trial.points.response_prediction_A(:,1))+1,trial.num);
data.trials.response_choice_A.data = nan(size(sktt.data,1),max(trial.points.response_choice_A(:,2)-trial.points.response_choice_A(:,1))+1,trial.num);
data.trials.response_prediction_B.data = nan(size(sktt.data,1),max(trial.points.response_prediction_B(:,2)-trial.points.response_prediction_B(:,1))+1,trial.num);
data.trials.response_choice_B.data = nan(size(sktt.data,1),max(trial.points.response_choice_B(:,2)-trial.points.response_choice_B(:,1))+1,trial.num);
%
for loop_trials = 1:trial.num
    data.trials.prediction.data(:,1:(trial.points.prediction(loop_trials,2) - trial.points.prediction(loop_trials,1)+1),loop_trials) = ...
        sktt.data(:,trial.points.prediction(loop_trials,1):trial.points.prediction(loop_trials,2));
    data.trials.choice.data(:,1:(trial.points.choice(loop_trials,2) - trial.points.choice(loop_trials,1)+1),loop_trials) = ...
        sktt.data(:,trial.points.choice(loop_trials,1):trial.points.choice(loop_trials,2));
    data.trials.evidence_other.data(:,1:(trial.points.evidence_other(loop_trials,2) - trial.points.evidence_other(loop_trials,1)+1),loop_trials) = ...
        sktt.data(:,trial.points.evidence_other(loop_trials,1):trial.points.evidence_other(loop_trials,2));
    data.trials.evidence_own.data(:,1:(trial.points.evidence_own(loop_trials,2) - trial.points.evidence_own(loop_trials,1)+1),loop_trials) = ...
        sktt.data(:,trial.points.evidence_own(loop_trials,1):trial.points.evidence_own(loop_trials,2));
end
for loop_trials_other = 1:trial.num_response_prediction_A
    data.trials.response_prediction_A.data(:,1:(trial.points.response_prediction_A(loop_trials_other,2) - trial.points.response_prediction_A(loop_trials_other,1)+1),loop_trials_other) = ...
        sktt.data(:,trial.points.response_prediction_A(loop_trials_other,1):trial.points.response_prediction_A(loop_trials_other,2));
end
for loop_trials_other = 1:trial.num_response_prediction_B
    data.trials.response_prediction_B.data(:,1:(trial.points.response_prediction_B(loop_trials_other,2) - trial.points.response_prediction_B(loop_trials_other,1)+1),loop_trials_other) = ...
        sktt.data(:,trial.points.response_prediction_B(loop_trials_other,1):trial.points.response_prediction_B(loop_trials_other,2));
end
for loop_trials_other = 1:trial.num_response_choice_A
    data.trials.response_choice_A.data(:,1:(trial.points.response_choice_A(loop_trials_other,2) - trial.points.response_choice_A(loop_trials_other,1)+1),loop_trials_other) = ...
        sktt.data(:,trial.points.response_choice_A(loop_trials_other,1):trial.points.response_choice_A(loop_trials_other,2));
end
for loop_trials_other = 1:trial.num_response_choice_B
    data.trials.response_choice_B.data(:,1:(trial.points.response_choice_B(loop_trials_other,2) - trial.points.response_choice_B(loop_trials_other,1)+1),loop_trials_other) = ...
        sktt.data(:,trial.points.response_choice_B(loop_trials_other,1):trial.points.response_choice_B(loop_trials_other,2));
end
%
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
if sktt.flag_bc_trial
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
    trial.points.bc.mean.prediction = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.choice = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.evidence_other = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.evidence_own = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.response_prediction_A = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.response_choice_A = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.response_prediction_B = nan(size(sktt.data,1),trial.num);
    trial.points.bc.mean.response_choice_B = nan(size(sktt.data,1),trial.num);
    for loop_trials = 1:trial.num
        trial.points.bc.mean.prediction(:,loop_trials) = mean(sktt.data(:,trial.points.bc.prediction(loop_trials,1):trial.points.bc.prediction(loop_trials,2)),2);
        trial.points.bc.mean.choice(:,loop_trials) = mean(sktt.data(:,trial.points.bc.choice(loop_trials,1):trial.points.bc.choice(loop_trials,2)),2);
        trial.points.bc.mean.evidence_other(:,loop_trials) = mean(sktt.data(:,trial.points.bc.evidence_other(loop_trials,1):trial.points.bc.evidence_other(loop_trials,2)),2);
        trial.points.bc.mean.evidence_own(:,loop_trials) = mean(sktt.data(:,trial.points.bc.evidence_own(loop_trials,1):trial.points.bc.evidence_own(loop_trials,2)),2);
    end
    for loop_trials_other = 1:trial.num_response_prediction_A
        trial.points.bc.mean.response_prediction_A(:,loop_trials_other) = mean(sktt.data(:,trial.points.bc.response_prediction_A(loop_trials_other,1):trial.points.bc.response_prediction_A(loop_trials_other,2)),2);
    end
    for loop_trials_other = 1:trial.num_response_prediction_B
        trial.points.bc.mean.response_prediction_B(:,loop_trials_other) = mean(sktt.data(:,trial.points.bc.response_prediction_B(loop_trials_other,1):trial.points.bc.response_prediction_B(loop_trials_other,2)),2);
    end
    for loop_trials_other = 1:trial.num_response_choice_A
        trial.points.bc.mean.response_choice_A(:,loop_trials_other) = mean(sktt.data(:,trial.points.bc.response_choice_A(loop_trials_other,1):trial.points.bc.response_choice_A(loop_trials_other,2)),2);
    end
    for loop_trials_other = 1:trial.num_response_choice_B
        trial.points.bc.mean.response_choice_B(:,loop_trials_other) = mean(sktt.data(:,trial.points.bc.response_choice_B(loop_trials_other,1):trial.points.bc.response_choice_B(loop_trials_other,2)),2);
    end
    % substracting the mean from the data
    for loop_trials = 1:trial.num
        data.trials.prediction.data(:,:,loop_trials) = data.trials.prediction.data(:,:,loop_trials) - repmat(trial.points.bc.mean.prediction(:,loop_trials),1,size(data.trials.prediction.data,2));
        data.trials.choice.data(:,:,loop_trials) = data.trials.choice.data(:,:,loop_trials) - repmat(trial.points.bc.mean.choice(:,loop_trials),1,size(data.trials.choice.data,2));
        data.trials.evidence_other.data(:,:,loop_trials) = data.trials.evidence_other.data(:,:,loop_trials) - repmat(trial.points.bc.mean.evidence_other(:,loop_trials),1,size(data.trials.evidence_other.data,2));
        data.trials.evidence_own.data(:,:,loop_trials) = data.trials.evidence_own.data(:,:,loop_trials) - repmat(trial.points.bc.mean.evidence_own(:,loop_trials),1,size(data.trials.evidence_own.data,2));
    end
    for loop_trials_other = 1:trial.num_response_prediction_A
        data.trials.response_prediction_A.data(:,:,loop_trials_other) = data.trials.response_prediction_A.data(:,:,loop_trials_other) - repmat(trial.points.bc.mean.response_prediction_A(:,loop_trials_other),1,size(data.trials.response_prediction_A.data,2));
    end
    for loop_trials_other = 1:trial.num_response_prediction_B
        data.trials.response_prediction_B.data(:,:,loop_trials_other) = data.trials.response_prediction_B.data(:,:,loop_trials_other) - repmat(trial.points.bc.mean.response_prediction_B(:,loop_trials_other),1,size(data.trials.response_prediction_B.data,2));
    end
    for loop_trials_other = 1:trial.num_response_choice_A
        data.trials.response_choice_A.data(:,:,loop_trials_other) = data.trials.response_choice_A.data(:,:,loop_trials_other) - repmat(trial.points.bc.mean.response_choice_A(:,loop_trials_other),1,size(data.trials.response_choice_A.data,2));
    end
    for loop_trials_other = 1:trial.num_response_choice_B
        data.trials.response_choice_B.data(:,:,loop_trials_other) = data.trials.response_choice_B.data(:,:,loop_trials_other) - repmat(trial.points.bc.mean.response_choice_B(:,loop_trials_other),1,size(data.trials.response_choice_B.data,2));
    end
    fprintf('Baseline correction of the epochs/trials ... Done.\n');
end
%%%


