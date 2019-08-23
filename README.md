# Contents
## [General idea](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#tiger_eeg_preprocessing)
## [Versions and toolboxes](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#the-default-versions-used-for-this-codebase)
## [EEG pre-processing protocol](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#eeg-pre-processing-protocol-1)

# tiger_eeg_preprocessing
The general idea of why this code was created.
The EEG pre-processing code for the tiger task. It consists of the Hyper-scanning EEG data that needs to be preprocessed.
###### [Back to Contents](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#Contents)

# The default versions used for this codebase
--> MATLAB version and toolboxes versions

+ MATLAB Version: 9.3.0.713579 (R2017b)

+ Operating System: Linux 4.9.0-8-amd64 #1 SMP Debian 4.9.110-3+deb9u5 (2018-09-30) x86_64

+ Java Version: Java 1.8.0_121-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

+ Signal Processing Toolbox                             Version 7.5         (R2017b)

+ fieldtrip-lite-20170706

+ eeglab14_1_2b

+ SASICA - Copyright (C) 2014  Maximilien Chaumon

+ AMICA 1.5 (c) Jason Palmer, University of California San Diego, 2015. (info: 'chmod 777 ...', admin permissions to the relevent file required to run)

+ eBridge - EB.Info.Version = '1.0.01' - Copyright © 2013-2014 by Daniel Alschuler

+ CleanLine (tmullen-cleanline-696a7181b7d0) Feb 8, 2012, Tim Mullen, SCCN/INC/UCSD Copyright (C) 2011
###### [Back to Contents](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#Contents)

# EEG pre-processing protocol
+ Reading the data (EEG data and the Behavioral data).

+ Adjust the specified triggers automatically according to the photo diode trigger values.

+ Downsampling the data from 1024 to 512 Hz or any preference but has to be powers of 2 exponent (0<power<10). [optional step, can be turned on ('1') or off ('0')]

+ Cut the data between the most extreme points (recognized by the beginning and end trigger).

+ Cut the data into the 2 participants, so we work on these two different datasets seperately.

+ Check for the bridging of the data for the individual participants and mark the channels. [optional step, can be turned on ('1') or off ('0')]

+ Filter the data. Using a line filter to get rid of the line noise at the EU standard 50Hz (can be changed to 60Hz if required). Then a the Butterworth bandpass FIR filter 4th order [1-90 Hz] and baseline correct it using the average of the session. [optional step, can be turned on ('1') or off ('0')]

+ Clean the data at this point. [optional step, can be turned on ('1') or off ('0')]

  + Check for flatlines and mark them up. (This is an automated artifact rejection function which ensures that the data contains no flat-lined channels. If it is found then the channel is marked.)

  + Check for noisy channels and mark them up. (This is an automated artifact rejection function which ensures that the data contains no channels that record only noise for extended periods of time. If channels with control signals are contained in the data these are usually also removed. The criterion is based on correlation: if a channel is decorrelated from all others (pairwise correlation < a given threshold), excluding a given fraction of most correlated channels -- and if this holds on for a sufficiently long fraction of the data set -- then the channel is marked.)

  + Reject the bad channels marked so far.

  + Implementation of the artifact subspace reconstruction for short-time burst algorithm (ASR). (This is an automated artifact rejection function that ensures that the data contains no events that have abnormally strong power; the subspaces on which those events occur are reconstructed (interpolated) based on the rest of the EEG signal during these time periods. The basic principle is to first find a section of data that represents clean "reference" EEG and to compute statistics on there. Then, the function goes over the whole data in a sliding window and finds the subspaces in which there is activity that is more than a few standard deviations away from the reference EEG (this threshold is a tunable parameter). Once the function has found the bad subspaces it will treat them as missing data and reconstruct their content using a mixing matrix that was calculated on the clean data.)

  + Remove incompletely repaired segments from the data. (This function cuts segments from the data which contain high-power artifacts. Specifically, only windows are retained which have less than a certain fraction of "bad" channels, where a channel is bad in a window if its power is above or below a given upper/lower threshold (in standard deviations from a robust estimate of the EEG power distribution in the channel).)

+ Interpolate the missing channels that were rejected so far.  This is to minimize a potential bias towards a particular brain region or hemisphere. Also here remove the extra EOG channels (before the re-referencing) .

+ Re-referencing the data to common average reference (CAR). [optional step, can be turned on ('1') or off ('0')]

+ perform baseline correction.

+ remove the interpolated channels before the ICA is performed. This is done to avoid any rank deficiency that came with the inperpolation of the missing/rejected channels.

+ Independent component analysis (ICA) taking the entire session. The trials individually do not have enough data points to perform an ICA independently. There is an option to perform the ICA using the 'RUNICA' or the 'AMICA' algorithm; By default the 'AMICA' option is selected. [optional step, can be turned on ('1') or off ('0')]

+ Automatic detection of the artifacts in the component space using the Guided Selection of ICA components for Artifact rejection (SASICA) algorithm.

+ Cut the data into the trials  (prediction, choice, etc...). (optional mean baseline correction possibility) [optional step, can be turned on ('1') or off ('0')]

+ save the data for each participant until this step. [optional step, can be turned on ('1') or off ('0')]

+ Option here is to work in the component space here / or / convert back using the inverse ICAweights to the sensor space.

+ Further statistics and analysis to be implemented from here on.
###### [Back to Contents](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#Contents)

***author: Saurabh Kumar*** 
