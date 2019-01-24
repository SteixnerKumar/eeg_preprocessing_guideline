# Contents
## [General idea](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#tiger_eeg_preprocessing)
## [Versions and toolboxes](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#the-default-versions-used-for-this-codebase)
## [EEG pre-processing protocol](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#EEG-pre-processing-protocol)

# tiger_eeg_preprocessing
The general idea of why this code was created.
The EEG pre-processing code for the tiger task. It consists of the Hyper-scanning EEG data that needs to be preprocessed.
[TOP](https://github.com/saurabhsay/tiger_eeg_preprocessing/blob/master/README.md#Contents)

# The default versions used for this codebase
--> MATLAB version and toolboxes versions
MATLAB Version: 9.3.0.713579 (R2017b)

Operating System: Linux 4.9.0-8-amd64 #1 SMP Debian 4.9.110-3+deb9u5 (2018-09-30) x86_64

Java Version: Java 1.8.0_121-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

Signal Processing Toolbox                             Version 7.5         (R2017b)

fieldtrip-lite-20170706

eeglab14_1_2b

SASICA - Copyright (C) 2014  Maximilien Chaumon

AMICA 1.5 (c) Jason Palmer, University of California San Diego, 2015.

eBridge - EB.Info.Version = '1.0.01' - Copyright © 2013-2014 by Daniel Alschuler

CleanLine (tmullen-cleanline-696a7181b7d0) Feb 8, 2012, Tim Mullen, SCCN/INC/UCSD Copyright (C) 2011

# EEG pre-processing protocol
1.) Reading the data (EEG data and the Behavioral data).

2.) Downsampling the data from 1024 to 512 Hz.

3.) Cut the data between the most extreme points (recognized by the beginning and end trigger).

4.) Cut the data into the 2 participants, so we work on these two different datasets seperately.

5.) Check for the bridging of the data for the individual participants and mark the channels.

5.5) Yet to include the automatic adjustment of the stimulus onset timings here based on the photo-diode triggers.

6.) Filter the data, Currently using the Butterworth bandpass FIR filter 4th order [1-45 Hz] and baseline correct it using the average of the session.

7.) Check for flatlines and mark them up.

8.) Check for noisy channels and mark them up.

9.) Reject the bad channels marked so far.

10.) Perform the ASR (artifact subspace reconstruction for short-time burst).

11.) Remove incompletely repaired segments from the data.

12.) Interpolate the missing channels that were rejected so far.

13.) Re-referencing the data to CAR (common average reference).

14.) Cut the data into the trials  (prediction, choice, etc...).

15.) perform baseline correction.

16.) remove the interpolated channels before the ICA is performed.

17.) ICA either trial-wise or the session. (AMICA algorithm)

18.) ICA automatic detection of artifacts implemented using the SASICA algorithm.

19.) Can work in the component space here / or / convert back the inverse ICAweights to come back into the sensor space.

20.) save the data for each participant until this step.

21.) Further statistics and analysis to be implemented from here on.


