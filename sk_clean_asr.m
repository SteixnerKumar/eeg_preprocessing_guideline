function signal = sk_clean_asr(signal,cutoff,windowlen,stepsize,maxdims,ref_maxbadchannels,ref_tolerances,ref_wndlen,usegpu)
% Run the ASR method on some high-pass filtered recording.
% Signal = clean_asr(Signal,StandardDevCutoff,WindowLength,BlockSize,MaxDimensions,ReferenceMaxBadChannels,RefTolerances,ReferenceWindowLength)
%
% This is an automated artifact rejection function that ensures that the data contains no events
% that have abnormally strong power; the subspaces on which those events occur are reconstructed 
% (interpolated) based on the rest of the EEG signal during these time periods.
%
% The basic principle is to first find a section of data that represents clean "reference" EEG and
% to compute statistics on there. Then, the function goes over the whole data in a sliding window
% and finds the subspaces in which there is activity that is more than a few standard deviations
% away from the reference EEG (this threshold is a tunable parameter). Once the function has found
% the bad subspaces it will treat them as missing data and reconstruct their content using a mixing
% matrix that was calculated on the clean data.
%
% Notes: 
%   This function by default attempts to use the Statistics toolbox in order to automatically
%   extract calibration data for use by ASR from the given recording. This step is automatically
%   skipped if no Statistics toolbox is present (then the entire recording will be used for
%   calibration, which is fine for mildly contaminated data -- see ReferenceMaxBadChannels below).
%
% In:
%   Signal : continuous data set, assumed to be *zero mean*, e.g., appropriately high-passed (e.g.
%            >0.5Hz or with a 0.5Hz - 1.0Hz transition band)
%
%   Cutoff : Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance
%            is larger than this threshold relative to the calibration data are considered missing
%            data and will be removed. The most aggressive value that can be used without losing
%            much EEG is 3. For new users it is recommended to at first visually inspect the difference 
%            between the original and cleaned data to get a sense of the removed content at various 
%            levels. A quite conservative value is 5. Default: 5.
%
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these better to your data.
%
%   WindowLength : Length of the statistcs window, in seconds. This should not be much longer 
%                  than the time scale over which artifacts persist, but the number of samples in
%                  the window should not be smaller than 1.5x the number of channels. Default:
%                  max(0.5,1.5*Signal.nbchan/Signal.srate);
%
%   StepSize : Step size for processing. The reprojection matrix will be updated every this many
%              samples and a blended matrix is used for the in-between samples. If empty this will
%              be set the WindowLength/2 in samples. Default: []
%
%   MaxDimensions : Maximum dimensionality to reconstruct. Up to this many dimensions (or up to this 
%                   fraction of dimensions) can be reconstructed for a given data segment. This is
%                   since the lower eigenvalues are usually not estimated very well. Default: 2/3.
%
%   ReferenceMaxBadChannels : If a number is passed in here, the ASR method will be calibrated based
%                             on sufficiently clean data that is extracted first from the recording
%                             that is then processed with ASR. This number is the maximum tolerated
%                             fraction of "bad" channels within a given time window of the recording
%                             that is considered acceptable for use as calibration data. Any data
%                             windows within the tolerance range are then used for calibrating the
%                             threshold statistics. Instead of a number one may also directly pass
%                             in a data set that contains calibration data (for example a minute of
%                             resting EEG) or the name of a data set in the workspace.
%
%                             If this is set to 'off', all data is used for calibration. This will
%                             work as long as the fraction of contaminated data is lower than the
%                             the breakdown point of the robust statistics in the ASR calibration
%                             (50%, where 30% of clearly recognizable artifacts is a better estimate
%                             of the practical breakdown point).
%
%                             A lower value makes this criterion more aggressive. Reasonable range:
%                             0.05 (very aggressive) to 0.3 (quite lax). If you have lots of little
%                             glitches in a few channels that don't get entirely cleaned you might
%                             want to reduce this number so that they don't go into the calibration
%                             data. Default: 0.075.
%                             
%
%   ReferenceTolerances : These are the power tolerances outside of which a channel in a
%                         given time window is considered "bad", in standard deviations relative to
%                         a robust EEG power distribution (lower and upper bound). Together with the
%                         previous parameter this determines how ASR calibration data is be
%                         extracted from a recording. Can also be specified as 'off' to achieve the
%                         same effect as in the previous parameter. Default: [-3.5 5.5].
%
%   ReferenceWindowLength : Granularity at which EEG time windows are extracted
%                           for calibration purposes, in seconds. Default: 1.
%
%   UseGPU : Whether to run on the GPU. This makes sense for offline processing if you have a a card with
%            enough memory and good double-precision performance (e.g., NVIDIA GTX Titan or K20). 
%            Note that for this to work you need to a) have the Parallel Computing toolbox and b) remove 
%            the dummy gather.m file from the path. Default: false
%
% Out:
%   Signal : data set with local peaks removed
%
% Examples:
%   % use the defaults
%   eeg = clean_asr(eeg);
%
%   % use a more aggressive threshold
%   eeg = clean_asr(eeg,2.5);
%
%   % disable subset selection of calibration data (use all data instead)
%   eeg = clean_asr(eeg,[],[],[],[],'off');
%
%   % use a custom calibration measurement (e.g., EEGLAB dataset containing a baseline recording)
%   eeg = clean_asr(eeg,[],[],[],[],mybaseline);
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-10-15

% Copyright (C) Christian Kothe, SCCN, 2012, ckothe@ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if ~exist('cutoff','var') || isempty(cutoff) cutoff = 5; end
if ~exist('windowlen','var') || isempty(windowlen) windowlen = max(0.5,1.5*signal.nbchan/signal.srate); end
if ~exist('stepsize','var') || isempty(stepsize) stepsize = []; end
if ~exist('maxdims','var') || isempty(maxdims) maxdims = 0.66; end
if ~exist('ref_maxbadchannels','var') || isempty(ref_maxbadchannels) ref_maxbadchannels = 0.075; end
if ~exist('ref_tolerances','var') || isempty(ref_tolerances) ref_tolerances = [-3.5 5.5]; end
if ~exist('ref_wndlen','var') || isempty(ref_wndlen) ref_wndlen = 1; end
if ~exist('usegpu','var') || isempty(usegpu) usegpu = false; end

signal.data = double(signal.data);

% first determine the reference (calibration) data
if isnumeric(ref_maxbadchannels) && isnumeric(ref_tolerances) && isnumeric(ref_wndlen)
    disp('Finding a clean section of the data...');
    try
        ref_section = clean_windows(signal,ref_maxbadchannels,ref_tolerances,ref_wndlen); 
    catch e
        disp('An error occurred while trying to identify a subset of clean calibration data from the recording.');
        disp('If this is because do not have EEGLAB loaded or no Statistics toolbox, you can generally');
        disp('skip this step by passing in ''off'' as the ReferenceMaxBadChannels parameter.');
        disp('Error details: ');
        hlp_handleerror(e,1);
        disp('Falling back to using the entire data for calibration.')
        ref_section = signal;
    end
elseif strcmp(ref_maxbadchannels,'off') || strcmp(ref_tolerances,'off') || strcmp(ref_wndlen,'off')
    disp('Using the entire data for calibration (reference parameters set to ''off'').')
    ref_section = signal;
elseif ischar(ref_maxbadchannels) && isvarname(ref_maxbadchannels)
    disp('Using a user-supplied data set in the workspace.');
    ref_section = evalin('base',ref_maxbadchannels);
elseif all(isfield(ref_maxbadchannels,{'data','srate','chanlocs'}))
    disp('Using a user-supplied clean section of data.');
    ref_section = ref_maxbadchannels; 
else
    error('Unsupported value for argument ref_maxbadchannels.');
end

% calibrate on the reference data
disp('Estimating calibration statistics; this may take a while...');
if exist('hlp_diskcache','file')
    state = hlp_diskcache('filterdesign',@asr_calibrate,ref_section.data,ref_section.srate,cutoff);
else
    state = asr_calibrate(ref_section.data,ref_section.srate,cutoff);
end
clear ref_section;

if isempty(stepsize)
    stepsize = floor(signal.srate*windowlen/2); end

% extrapolate last few samples of the signal
sig = [signal.data bsxfun(@minus,2*signal.data(:,end),signal.data(:,(end-1):-1:end-round(windowlen/2*signal.srate)))];
% process signal using ASR
[signal.data,state] = asr_process(sig,signal.srate,state,windowlen,windowlen/2,stepsize,maxdims,[],usegpu);
% shift signal content back (to compensate for processing delay)
signal.data(:,1:size(state.carry,2)) = [];

end




function result = hlp_memfree
% Get the amount of free physical memory, in bytes
result = java.lang.management.ManagementFactory.getOperatingSystemMXBean().getFreePhysicalMemorySize();
end

function y = geometric_median(X,tol,y,max_iter)
% Calculate the geometric median for a set of observations (mean under a Laplacian noise distribution)
% This is using Weiszfeld's algorithm.
%
% In:
%   X : the data, as in mean
%   tol : tolerance (default: 1.e-5)
%   y : initial value (default: median(X))
%   max_iter : max number of iterations (default: 500)
%
% Out:
%   g : geometric median over X

if ~exist('tol','var') || isempty(tol)
    tol = 1.e-5; end
if ~exist('y','var') || isempty(y)
    y = median(X); end
if ~exist('max_iter','var') || isempty(max_iter)
    max_iter = 500; end

for i=1:max_iter
    invnorms = 1./sqrt(sum(bsxfun(@minus,X,y).^2,2));
    [y,oldy] = deal(sum(bsxfun(@times,X,invnorms)) / sum(invnorms),y);
    if norm(y-oldy)/norm(y) < tol
        break; end
end
end

function y = block_geometric_median(X,blocksize,varargin)
% Calculate a blockwise geometric median (faster and less memory-intensive 
% than the regular geom_median function).
%
% This statistic is not robust to artifacts that persist over a duration that
% is significantly shorter than the blocksize.
%
% In:
%   X : the data (#observations x #variables)
%   blocksize : the number of successive samples over which a regular mean 
%               should be taken
%   tol : tolerance (default: 1.e-5)
%   y : initial value (default: median(X))
%   max_iter : max number of iterations (default: 500)
%
% Out:
%   g : geometric median over X
%
% Notes:
%   This function is noticably faster if the length of the data is divisible by the block size.
%   Uses the GPU if available.
% 

if nargin < 2 || isempty(blocksize)
    blocksize = 1; end

if blocksize > 1
    [o,v] = size(X);       % #observations & #variables
    r = mod(o,blocksize);  % #rest in last block
    b = (o-r)/blocksize;   % #blocks
    if r > 0
        X = [reshape(sum(reshape(X(1:(o-r),:),blocksize,b*v)),b,v); sum(X((o-r+1):end,:))*(blocksize/r)];
    else
        X = reshape(sum(reshape(X,blocksize,b*v)),b,v);
    end
end

try
    y = gather(geometric_median(gpuArray(X),varargin{:}))/blocksize;
catch
    y = geometric_median(X,varargin{:})/blocksize;
end
end


function state = asr_calibrate(X,srate,cutoff,blocksize,B,A,window_len,window_overlap,max_dropout_fraction,min_clean_fraction)
% Calibration function for the Artifact Subspace Reconstruction (ASR) method.
% State = asr_calibrate(Data,SamplingRate,Cutoff,BlockSize,FilterB,FilterA,WindowLength,WindowOverlap,MaxDropoutFraction,MinCleanFraction)
%
% The input to this data is a multi-channel time series of calibration data. In typical uses the
% calibration data is clean resting EEG data of ca. 1 minute duration (can also be longer). One can
% also use on-task data if the fraction of artifact content is below the breakdown point of the
% robust statistics used for estimation (50% theoretical, ~30% practical). If the data has a
% proportion of more than 30-50% artifacts then bad time windows should be removed beforehand. This
% data is used to estimate the thresholds that are used by the ASR processing function to identify
% and remove artifact components.
%
% The calibration data must have been recorded for the same cap design from which data for cleanup
% will be recorded, and ideally should be from the same session and same subject, but it is possible
% to reuse the calibration data from a previous session and montage to the extent that the cap is
% placed in the same location (where loss in accuracy is more or less proportional to the mismatch
% in cap placement).
%
% The calibration data should have been high-pass filtered (for example at 0.5Hz or 1Hz using a
% Butterworth IIR filter).
%
% In:
%   Data : Calibration data [#channels x #samples]; *zero-mean* (e.g., high-pass filtered) and
%          reasonably clean EEG of not much less than 30 seconds length (this method is typically
%          used with 1 minute or more).
%
%   SamplingRate : Sampling rate of the data, in Hz.
%
%
%   The following are optional parameters (the key parameter of the method is the RejectionCutoff):
%
%   RejectionCutoff: Standard deviation cutoff for rejection. Data portions whose variance is larger
%                    than this threshold relative to the calibration data are considered missing
%                    data and will be removed. The most aggressive value that can be used without
%                    losing too much EEG is 2.5. A quite conservative value would be 5. Default: 5.
%
%   Blocksize : Block size for calculating the robust data covariance and thresholds, in samples;
%               allows to reduce the memory and time requirements of the robust estimators by this 
%               factor (down to Channels x Channels x Samples x 16 / Blocksize bytes). Default: 10
%
%   FilterB, FilterA : Coefficients of an IIR filter that is used to shape the spectrum of the signal
%                      when calculating artifact statistics. The output signal does not go through
%                      this filter. This is an optional way to tune the sensitivity of the algorithm
%                      to each frequency component of the signal. The default filter is less
%                      sensitive at alpha and beta frequencies and more sensitive at delta (blinks)
%                      and gamma (muscle) frequencies. Default: 
%                      [b,a] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
%
%   WindowLength : Window length that is used to check the data for artifact content. This is 
%                  ideally as long as the expected time scale of the artifacts but short enough to 
%				   allow for several 1000 windows to compute statistics over. Default: 0.5.
%
%   WindowOverlap : Window overlap fraction. The fraction of two successive windows that overlaps.
%                   Higher overlap ensures that fewer artifact portions are going to be missed (but
%                   is slower). Default: 0.66
%
%   MaxDropoutFraction : Maximum fraction of windows that can be subject to signal dropouts 
%                        (e.g., sensor unplugged), used for threshold estimation. Default: 0.1
%
%   MinCleanFraction : Minimum fraction of windows that need to be clean, used for threshold
%                      estimation. Default: 0.25
%
%
% Out:
%   State : initial state struct for asr_process
%
% Notes:
%   This can run on a GPU with large memory and good double-precision performance for faster processing 
%   (e.g., on an NVIDIA GTX Titan or K20), but requires that the Parallel Computing toolbox is
%   installed.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-31

% asr_calibrate_version<1.03> -- for the cache

% UC Copyright Notice
% This software is Copyright (C) 2013 The Regents of the University of California. All Rights Reserved.
% 
% Permission to copy, modify, and distribute this software and its documentation for educational,
% research and non-profit purposes, without fee, and without a written agreement is hereby granted,
% provided that the above copyright notice, this paragraph and the following three paragraphs appear
% in all copies.
% 
% Permission to make commercial use of this software may be obtained by contacting:
% Technology Transfer Office
% 9500 Gilman Drive, Mail Code 0910
% University of California
% La Jolla, CA 92093-0910
% (858) 534-5815
% invent@ucsd.edu 
% 
% This software program and documentation are copyrighted by The Regents of the University of
% California. The software program and documentation are supplied "as is", without any accompanying
% services from The Regents. The Regents does not warrant that the operation of the program will be
% uninterrupted or error-free. The end-user understands that the program was developed for research
% purposes and is advised not to rely exclusively on the program for any reason.
% 
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
% CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS.

[C,S] = size(X);

if nargin < 3 || isempty(cutoff)
    cutoff = 5; end
if nargin < 4 || isempty(blocksize)
    blocksize = 10; end
blocksize = max(blocksize,ceil((C*C*S*8*3*2)/hlp_memfree));
if nargin < 6 || isempty(A) || isempty(B)
    try
        % try to use yulewalk to design the filter (Signal Processing toolbox required)
        [B,A] = yulewalk(8,[[0 2 3 13 16 40 min(80,srate/2-1)]*2/srate 1],[3 0.75 0.33 0.33 1 1 3 3]);
    catch e %#ok<NASGU>
        % yulewalk not available (maybe no toolbox installed) -- use precomputed filter
        % coefficients depending on sampling rate
        switch srate
            case 100
                [B,A] = deal([0.9314233528641650 -1.0023683814963549 -0.4125359862018213  0.7631567476327510  0.4160430392910331 -0.6549131038692215 -0.0372583518046807  0.1916268458752655  0.0462411971592346],[1.0000000000000000 -0.4544220180303844 -1.0007038682936749  0.5374925521337940  0.4905013360991340 -0.4861062879351137 -0.1995986490699414  0.1830048420730026  0.0457678549234644]);
            case 128
                [B,A] = deal([1.1027301639165037 -2.0025621813611867  0.8942119516481342  0.1549979524226999  0.0192366904488084  0.1782897770278735 -0.5280306696498717  0.2913540603407520 -0.0262209802526358],[1.0000000000000000 -1.1042042046423233 -0.3319558528606542  0.5802946221107337 -0.0010360013915635  0.0382167091925086 -0.2609928034425362  0.0298719057761086  0.0935044692959187]);
            case 200
                [B,A] = deal([1.4489483325802353 -2.6692514764802775  2.0813970620731115 -0.9736678877049534  0.1054605060352928 -0.1889101692314626  0.6111331636592364 -0.3616483013075088  0.1834313060776763],[1.0000000000000000 -0.9913236099393967  0.3159563145469344 -0.0708347481677557 -0.0558793822071149 -0.2539619026478943  0.2473056615251193 -0.0420478437473110  0.0077455718334464]);
            case 256
                [B,A] = deal([1.7587013141770287 -4.3267624394458641  5.7999880031015953 -6.2396625463547508  5.3768079046882207 -3.7938218893374835  2.1649108095226470 -0.8591392569863763  0.2569361125627988],[1.0000000000000000 -1.7008039639301735  1.9232830391058724 -2.0826929726929797  1.5982638742557307 -1.0735854183930011  0.5679719225652651 -0.1886181499768189  0.0572954115997261]);
            case 300
                [B,A] = deal([1.9153920676433143  -5.7748421104926795   9.1864764859103936 -10.7350356619363630   9.6423672437729007  -6.6181939699544277   3.4219421494177711  -1.2622976569994351   0.2968423019363821],[1.0000000000000000 -2.3143703322055491  3.2222567327379434 -3.6030527704320621  2.9645154844073698 -1.8842615840684735  0.9222455868758080 -0.3103251703648485  0.0634586449896364]);
            case 500
                [B,A] = deal([2.3133520086975823 -11.9471223009159130  29.1067166493384340 -43.7550171007238190  44.3385767452216370 -30.9965523846388000  14.6209883020737190  -4.2743412400311449   0.5982553583777899],[1.0000000000000000  -4.6893329084452580  10.5989986701080210 -14.9691518101365230  14.3320358399731820  -9.4924317069169977   4.2425899618982656  -1.1715600975178280   0.1538048427717476]);
            case 512
                [B,A] = deal([2.3275475636130865 -12.2166478485960430  30.1632789058248850 -45.8009842020820410  46.7261263011068880 -32.7796858196767220  15.4623349612560630  -4.5019779685307473   0.6242733481676324],[1.0000000000000000  -4.7827378944258703  10.9780696236622980 -15.6795187888195360  15.1281978667576310 -10.0632079834518220   4.5014690636505614  -1.2394100873286753   0.1614727510688058]);
            otherwise
                error('repair_bursts:NoYulewalk','The yulewalk() function was not found and there is no pre-computed spectral filter for your sampling rate. If you would like to use the default spectral filter please try to resample to one of the supported rates (100,128,200,256,300,500,512) or get the appropriate toobox license (you can also disable the spectral weighting feature or supply your own precalculated IIR filter coefficients).');
        end
    end
end
if nargin < 8 || isempty(window_len)
    window_len = 0.5; end
if nargin < 9 || isempty(window_overlap)
    window_overlap = 0.66; end
if nargin < 10 || isempty(max_dropout_fraction)
    max_dropout_fraction = 0.1; end
if nargin < 11 || isempty(min_clean_fraction)
    min_clean_fraction = 0.25; end

X(~isfinite(X(:))) = 0;

% apply the signal shaping filter and initialize the IIR filter state
[X,iirstate] = filter(B,A,double(X),[],2); X = X';
if any(~isfinite(X(:)))
    error('The IIR filter diverged on your data. Please try using either a more conservative filter or removing some bad sections/channels from the calibration data.'); end

% calculate the sample covariance matrices U (averaged in blocks of blocksize successive samples)
U = zeros(length(1:blocksize:S),C*C);
for k=1:blocksize
    range = min(S,k:blocksize:(S+k-1));
    U = U + reshape(bsxfun(@times,reshape(X(range,:),[],1,C),reshape(X(range,:),[],C,1)),size(U));
end

% get the mixing matrix M
M = sqrtm(real(reshape(block_geometric_median(U/blocksize),C,C)));

% window length for calculating thresholds
N = round(window_len*srate);

% get the threshold matrix T
fprintf('Determining per-component thresholds...');
[V,D] = eig(M); %#ok<NASGU>
X = abs(X*V);
for c = C:-1:1
    % compute RMS amplitude for each window...
    rms = X(:,c).^2;
    rms = sqrt(sum(rms(bsxfun(@plus,round(1:N*(1-window_overlap):S-N),(0:N-1)')))/N);
    % fit a distribution to the clean part
    [mu(c),sig(c)] = fit_eeg_distribution(rms,min_clean_fraction,max_dropout_fraction);
end
T = diag(mu + cutoff*sig)*V';
disp('done.');

% initialize the remaining filter state
state = struct('M',M,'T',T,'B',B,'A',A,'cov',[],'carry',[],'iir',iirstate,'last_R',[],'last_trivial',true);
end

function [X,Zf] = moving_average(N,X,Zi)
% Run a moving-average filter along the second dimension of the data.
% [X,Zf] = moving_average(N,X,Zi)
%
% In:
%   N : filter length in samples
%   X : data matrix [#Channels x #Samples]
%   Zi : initial filter conditions (default: [])
%
% Out:
%   X : the filtered data
%   Zf : final filter conditions
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2012-01-10

if nargin <= 2 || isempty(Zi)
    Zi = zeros(size(X,1),N); end

% pre-pend initial state & get dimensions
Y = [Zi X]; M = size(Y,2);
% get alternating index vector (for additions & subtractions)
I = [1:M-N; 1+N:M];
% get sign vector (also alternating, and includes the scaling)
S = [-ones(1,M-N); ones(1,M-N)]/N;
% run moving average
X = cumsum(bsxfun(@times,Y(:,I(:)),S(:)'),2);
% read out result
X = X(:,2:2:end);

if nargout > 1
    Zf = [-(X(:,end)*N-Y(:,end-N+1)) Y(:,end-N+2:end)]; end
end


function [outdata,outstate] = asr_process(data,srate,state,windowlen,lookahead,stepsize,maxdims,maxmem,usegpu)
% Processing function for the Artifact Subspace Reconstruction (ASR) method.
% [Data,State] = asr_process(Data,SamplingRate,State,WindowLength,LookAhead,StepSize,MaxDimensions,MaxMemory,UseGPU)
%
% This function is used to clean multi-channel signal using the ASR method. The required inputs are 
% the data matrix, the sampling rate of the data, and the filter state (as initialized by
% asr_calibrate). If the data is used on successive chunks of data, the output state of the previous 
% call to asr_process should be passed in.
%
% In:
%   Data : Chunk of data to process [#channels x #samples]. This is a chunk of data, assumed to be
%          a continuation of the data that was passed in during the last call to asr_process (if
%          any). The data should be *zero-mean* (e.g., high-pass filtered the same way as for
%          asr_calibrate).
%   
%   SamplingRate : sampling rate of the data in Hz (e.g., 250.0)
%
%   State : initial filter state (determined by asr_calibrate or from previous call to asr_process)
%
%   WindowLength : Length of the statistcs window, in seconds (e.g., 0.5). This should not be much
%                  longer than the time scale over which artifacts persist, but the number of samples 
%                  in the window should not be smaller than 1.5x the number of channels. Default: 0.5
%
%   LookAhead : Amount of look-ahead that the algorithm should use. Since the processing is causal,
%               the output signal will be delayed by this amount. This value is in seconds and should
%               be between 0 (no lookahead) and WindowLength/2 (optimal lookahead). The recommended
%               value is WindowLength/2. Default: WindowLength/2
%
%   StepSize : The statistics will be updated every this many samples. The larger this is, the faster 
%              the algorithm will be. The value must not be larger than WindowLength*SamplingRate.
%              The minimum value is 1 (update for every sample) while a good value is 1/3 of a second.
%              Note that an update is always performed also on the first and last sample of the data
%              chunk. Default: 32
%
%   MaxDimensions : Maximum dimensionality of artifacts to remove. Up to this many dimensions (or up 
%                   to this fraction of dimensions) can be removed for a given data segment. If the
%                   algorithm needs to tolerate extreme artifacts a higher value than the default
%                   may be used (the maximum fraction is 1.0). Default 0.66
%
%   MaxMemory : The maximum amount of memory used by the algorithm when processing a long chunk with
%               many channels, in MB. The recommended value is at least 256. To run on the GPU, use
%               the amount of memory available to your GPU here (needs the parallel computing toolbox).
%               default: min(5000,1/2 * free memory in MB). Using smaller amounts of memory leads to
%               longer running times.
%
%   UseGPU : Whether to run on the GPU. This makes sense for offline processing if you have a a card
%            with enough memory and good double-precision performance (e.g., NVIDIA GTX Titan or
%            K20). Note that for this to work you need to have the Parallel Computing toolbox.
%            Default: false
%
% Out:
%   Data : cleaned data chunk (same length as input but delayed by LookAhead samples)
%
%   State : final filter state (can be passed in for subsequent calls)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-08-31

% UC Copyright Notice
% This software is Copyright (C) 2013 The Regents of the University of California. All Rights Reserved.
% 
% Permission to copy, modify, and distribute this software and its documentation for educational,
% research and non-profit purposes, without fee, and without a written agreement is hereby granted,
% provided that the above copyright notice, this paragraph and the following three paragraphs appear
% in all copies.
% 
% Permission to make commercial use of this software may be obtained by contacting:
% Technology Transfer Office
% 9500 Gilman Drive, Mail Code 0910
% University of California
% La Jolla, CA 92093-0910
% (858) 534-5815
% invent@ucsd.edu 
% 
% This software program and documentation are copyrighted by The Regents of the University of
% California. The software program and documentation are supplied "as is", without any accompanying
% services from The Regents. The Regents does not warrant that the operation of the program will be
% uninterrupted or error-free. The end-user understands that the program was developed for research
% purposes and is advised not to rely exclusively on the program for any reason.
% 
% IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
% CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS.

if nargin < 4 || isempty(windowlen) 
    windowlen = 0.5; end
windowlen = max(windowlen,1.5*size(data,1)/srate);
if nargin < 5 || isempty(lookahead)
    lookahead = windowlen/2; end
if nargin < 6 || isempty(stepsize)
    stepsize = 32; end
if nargin < 7 || isempty(maxdims)
    maxdims = 0.66; end
if nargin < 9 || isempty(usegpu)
    usegpu = false; end
if nargin < 8 || isempty(maxmem)
    if usegpu
        dev = gpuDevice(); maxmem = dev.FreeMemory/2^20; 
    else
        maxmem = hlp_memfree/(2^21);
    end
end
if maxdims < 1
    maxdims = round(size(data,1)*maxdims); end
if isempty(data)
    outdata = data; outstate = state; return; end

[C,S] = size(data);
N = round(windowlen*srate);
P = round(lookahead*srate);
[T,M,A,B] = deal(state.T,state.M,state.A,state.B);

% initialize prior filter state by extrapolating available data into the past (if necessary)
if isempty(state.carry)
    state.carry = repmat(2*data(:,1),1,P) - data(:,1+mod(((P+1):-1:2)-1,S)); end

data = [state.carry data];
data(~isfinite(data(:))) = 0;

% split up the total sample range into k chunks that will fit in memory
splits = ceil((C*C*S*8*8 + C*C*8*S/stepsize + C*S*8*2 + S*8*5) / (maxmem*1024*1024 - C*C*P*8*3));
if splits > 1
    fprintf('Now cleaning data in %i blocks',splits); end

for i=1:splits
    range = 1+floor((i-1)*S/splits) : min(S,floor(i*S/splits));
    if ~isempty(range)
        % get spectrally shaped data X for statistics computation (range shifted by lookahead)
        [X,state.iir] = filter(B,A,double(data(:,range+P)),state.iir,2);
        % move it to the GPU if applicable
        if usegpu && length(range) > 1000
            try X = gpuArray(X); catch,end; end
        % compute running mean covariance (assuming a zero-mean signal)
        [Xcov,state.cov] = moving_average(N,reshape(bsxfun(@times,reshape(X,1,C,[]),reshape(X,C,1,[])),C*C,[]),state.cov);
        % extract the subset of time points at which we intend to update
        update_at = min(stepsize:stepsize:(size(Xcov,2)+stepsize-1),size(Xcov,2));
        % if there is no previous R (from the end of the last chunk), we estimate it right at the first sample
        if isempty(state.last_R)
            update_at = [1 update_at]; 
            state.last_R = eye(C);
        end
        Xcov = reshape(Xcov(:,update_at),C,C,[]);
        if usegpu
            Xcov = gather(Xcov); end
        % do the reconstruction in intervals of length stepsize (or shorter if at the end of a chunk)
        last_n = 0;
        for j=1:length(update_at)
            % do a PCA to find potential artifact components
            [V,D] = eig(Xcov(:,:,j));
            [D,order] = sort(reshape(diag(D),1,C)); V = V(:,order);
            % determine which components to keep (variance below directional threshold or not admissible for rejection)
            keep = D<sum((T*V).^2) | (1:C)<(C-maxdims);
            trivial = all(keep);
            % update the reconstruction matrix R (reconstruct artifact components using the mixing matrix)
            if ~trivial
                R = real(M*pinv(bsxfun(@times,keep',V'*M))*V');
            else
                R = eye(C);
            end
            % apply the reconstruction to intermediate samples (using raised-cosine blending)
            n = update_at(j);
            if ~trivial || ~state.last_trivial
                subrange = range((last_n+1):n);
                blend = (1-cos(pi*(1:(n-last_n))/(n-last_n)))/2;
                data(:,subrange) = bsxfun(@times,blend,R*data(:,subrange)) + bsxfun(@times,1-blend,state.last_R*data(:,subrange));
            end
            [last_n,state.last_R,state.last_trivial] = deal(n,R,trivial);
        end
    end
    if splits > 1
        fprintf('.'); end
end
if splits > 1
    fprintf('\n'); end

% carry the look-ahead portion of the data over to the state (for successive calls)
state.carry = [state.carry data(:,(end-P+1):end)];
state.carry = state.carry(:,(end-P+1):end);

% finalize outputs
outdata = data(:,1:(end-P));
if usegpu
    state.iir = gather(state.iir);
    state.cov = gather(state.cov);
end
outstate = state;
end


function [signal,sample_mask] = clean_windows(signal,max_bad_channels,zthresholds,window_len,window_overlap,max_dropout_fraction,min_clean_fraction,truncate_quant,step_sizes,shape_range)
if ~exist('max_bad_channels','var') || isempty(max_bad_channels) max_bad_channels = 0.2; end
if ~exist('zthresholds','var') || isempty(zthresholds) zthresholds = [-3.5 5]; end
if ~exist('window_len','var') || isempty(window_len) window_len = 1; end
if ~exist('window_overlap','var') || isempty(window_overlap) window_overlap = 0.66; end
if ~exist('max_dropout_fraction','var') || isempty(max_dropout_fraction) max_dropout_fraction = 0.1; end
if ~exist('min_clean_fraction','var') || isempty(min_clean_fraction) min_clean_fraction = 0.25; end
if ~exist('truncate_quant','var') || isempty(truncate_quant) truncate_quant = [0.022 0.6]; end
if ~exist('step_sizes','var') || isempty(step_sizes) step_sizes = [0.01 0.01]; end
if ~exist('shape_range','var') || isempty(shape_range) shape_range = 1.7:0.15:3.5; end
if ~isempty(max_bad_channels) && max_bad_channels > 0 && max_bad_channels < 1 %#ok<*NODEF>
    max_bad_channels = round(size(signal.data,1)*max_bad_channels); end

signal.data = double(signal.data);
[C,S] = size(signal.data);
N = window_len*signal.srate;
wnd = 0:N-1;
offsets = round(1:N*(1-window_overlap):S-N);

fprintf('Determining time window rejection thresholds...');
% for each channel...
for c = C:-1:1
    % compute RMS amplitude for each window...
    X = signal.data(c,:).^2;
    X = sqrt(sum(X(bsxfun(@plus,offsets,wnd')))/N);
    % robustly fit a distribution to the clean EEG part
    [mu,sig] = fit_eeg_distribution(X, ...
        min_clean_fraction, max_dropout_fraction, ...
        truncate_quant, step_sizes,shape_range);
    % calculate z scores relative to that
    wz(c,:) = (X - mu)/sig;
end
disp('done.');

% sort z scores into quantiles
swz = sort(wz);
% determine which windows to remove
remove_mask = false(1,size(swz,2));
if max(zthresholds)>0
    remove_mask(swz(end-max_bad_channels,:) > max(zthresholds)) = true; end
if min(zthresholds)<0
    remove_mask(swz(1+max_bad_channels,:) < min(zthresholds)) = true; end
removed_windows = find(remove_mask);

% find indices of samples to remove
removed_samples = repmat(offsets(removed_windows)',1,length(wnd))+repmat(wnd,length(removed_windows),1);
% mask them out
sample_mask = true(1,S); 
sample_mask(removed_samples(:)) = false;
fprintf('Keeping %.1f%% (%.0f seconds) of the data.\n',100*(mean(sample_mask)),nnz(sample_mask)/signal.srate);
% determine intervals to retain
% retain_data_intervals = reshape(find(diff([false sample_mask false])),2,[])';
% retain_data_intervals(:,2) = retain_data_intervals(:,2)-1;
signal.data = signal.data(:,sample_mask);
end

function [mu,sig,alpha,beta] = fit_eeg_distribution(X,min_clean_fraction,max_dropout_fraction,quants,step_sizes,beta)
% Estimate the mean and standard deviation of clean EEG from contaminated data.
% [Mu,Sigma,Alpha,Beta] = fit_eeg_distribution(X,MinCleanFraction,MaxDropoutFraction,FitQuantiles,StepSizes,ShapeRange)
%
% This function estimates the mean and standard deviation of clean EEG from a sample of amplitude
% values (that have preferably been computed over short windows) that may include a large fraction
% of contaminated samples. The clean EEG is assumed to represent a generalized Gaussian component in
% a mixture with near-arbitrary artifact components. By default, at least 25% (MinCleanFraction) of
% the data must be clean EEG, and the rest can be contaminated. No more than 10%
% (MaxDropoutFraction) of the data is allowed to come from contaminations that cause lower-than-EEG
% amplitudes (e.g., sensor unplugged). There are no restrictions on artifacts causing
% larger-than-EEG amplitudes, i.e., virtually anything is handled (with the exception of a very
% unlikely type of distribution that combines with the clean EEG samples into a larger symmetric
% generalized Gaussian peak and thereby "fools" the estimator). The default parameters should be
% fine for a wide range of settings but may be adapted to accomodate special circumstances.
%
% The method works by fitting a truncated generalized Gaussian whose parameters are constrained by
% MinCleanFraction, MaxDropoutFraction, FitQuantiles, and ShapeRange. The alpha and beta parameters
% of the gen. Gaussian are also returned. The fit is performed by a grid search that always finds a
% close-to-optimal solution if the above assumptions are fulfilled.
%
% In:
%   X : vector of amplitude values of EEG, possible containing artifacts
%       (coming from single samples or windowed averages)
%
%   MinCleanFraction : Minimum fraction of values in X that needs to be clean
%                      (default: 0.25)
%
%   MaxDropoutFraction : Maximum fraction of values in X that can be subject to
%                        signal dropouts (e.g., sensor unplugged) (default: 0.1)
%
%   FitQuantiles : Quantile range [lower,upper] of the truncated generalized Gaussian distribution
%                  that shall be fit to the EEG contents (default: [0.022 0.6])
%
%   StepSizes : Step size of the grid search; the first value is the stepping of the lower bound
%               (which essentially steps over any dropout samples), and the second value
%               is the stepping over possible scales (i.e., clean-data quantiles)
%               (default: [0.01 0.01])
%
%   ShapeRange : Range that the clean EEG distribution's shape parameter beta may take (default:
%                1.7:0.15:3.5)
%
% Out:
%   Mu : estimated mean of the clean EEG distribution
%
%   Sigma : estimated standard deviation of the clean EEG distribution
%
%   Alpha : estimated scale parameter of the generalized Gaussian clean EEG distribution (optional)
%
%   Beta : estimated shape parameter of the generalized Gaussian clean EEG distribution (optional)

% assign defaults
if ~exist('min_clean_fraction','var') || isempty(min_clean_fraction)
    min_clean_fraction = 0.25; end
if ~exist('max_dropout_fraction','var') || isempty(max_dropout_fraction)
    max_dropout_fraction = 0.1; end
if ~exist('quants','var') || isempty(quants)
    quants = [0.022 0.6]; end
if ~exist('step_sizes','var') || isempty(step_sizes)
    step_sizes = [0.01 0.01]; end
if ~exist('beta','var') || isempty(beta)
    beta = 1.7:0.15:3.5; end

% sanity checks
if ~isvector(quants) || numel(quants) > 2
    error('Fit quantiles needs to be a 2-element vector (support for matrices deprecated).'); end
if any(quants(:)<0) || any(quants(:)>1)
    error('Unreasonable fit quantiles.'); end
if any(step_sizes<0.0001) || any(step_sizes>0.1)
    error('Unreasonable step sizes.'); end
if any(beta>=7) || any(beta<=1)
    error('Unreasonable shape range.'); end

% sort data so we can access quantiles directly
X = double(sort(X(:)));
n = length(X);

% calc z bounds for the truncated standard generalized Gaussian pdf and pdf rescaler
for b=1:length(beta)    
    zbounds{b} = sign(quants-1/2).*gammaincinv(sign(quants-1/2).*(2*quants-1),1/beta(b)).^(1/beta(b)); %#ok<*AGROW>
    rescale(b) = beta(b)/(2*gamma(1/beta(b)));
end

% determine the quantile-dependent limits for the grid search
lower_min = min(quants);                    % we can generally skip the tail below the lower quantile
max_width = diff(quants);                   % maximum width is the fit interval if all data is clean
min_width = min_clean_fraction*max_width;   % minimum width of the fit interval, as fraction of data

% get matrix of shifted data ranges
X = X(bsxfun(@plus,(1:round(n*max_width))',round(n*(lower_min:step_sizes(1):lower_min+max_dropout_fraction))));
X1 = X(1,:); X = bsxfun(@minus,X,X1);

opt_val = Inf;
% for each interval width...
for m = round(n*(max_width:-step_sizes(2):min_width))
    % scale and bin the data in the intervals
    nbins = round(3*log2(1+m/2));
    H = bsxfun(@times,X(1:m,:),nbins./X(m,:));
    logq = log(histc(H,[0:nbins-1,Inf]) + 0.01);
    
    % for each shape value...
    for b=1:length(beta)
        bounds = zbounds{b};
        % evaluate truncated generalized Gaussian pdf at bin centers
        x = bounds(1)+(0.5:(nbins-0.5))/nbins*diff(bounds);
        p = exp(-abs(x).^beta(b))*rescale(b); p=p'/sum(p);
        
        % calc KL divergences
        kl = sum(bsxfun(@times,p,bsxfun(@minus,log(p),logq(1:end-1,:)))) + log(m);
        
        % update optimal parameters
        [min_val,idx] = min(kl);
        if min_val < opt_val
            opt_val = min_val;
            opt_beta = beta(b);
            opt_bounds = bounds;
            opt_lu = [X1(idx) X1(idx)+X(m,idx)];
        end
    end
end

% recover distribution parameters at optimum
alpha = (opt_lu(2)-opt_lu(1))/diff(opt_bounds);
mu = opt_lu(1)-opt_bounds(1)*alpha;
beta = opt_beta;

% calculate the distribution's standard deviation from alpha and beta
sig = sqrt((alpha^2)*gamma(3/beta)/gamma(1/beta));
end