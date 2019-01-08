function signal = sk_clean_channels(signal,min_corr,ignored_quantile,window_len,max_broken_time,linenoise_aware)
% Remove channels with abnormal data from a continuous data set.
% Signal = clean_channels(Signal,MinCorrelation,IgnoredQuantile,WindowLength,MaxBrokenTime,LineNoiseAware)
%
% This is an automated artifact rejection function which ensures that the data contains no channels
% that record only noise for extended periods of time. If channels with control signals are
% contained in the data these are usually also removed. The criterion is based on correlation: if a
% channel is decorrelated from all others (pairwise correlation < a given threshold), excluding a
% given fraction of most correlated channels -- and if this holds on for a sufficiently long fraction 
% of the data set -- then the channel is removed.
%
% In:
%   Signal          : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%                     with a 0.5Hz - 2.0Hz transition band).
%
%   MinCorrelation  : Minimum correlation between a channel and any other channel (in a short period 
%                     of time) below which the channel is considered abnormal for that time period.
%                     Reasonable range: 0.4 (very lax) to 0.6 (quite aggressive). The default is 0.45. 
%                     
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these to your data.
%   
%   IgnoredQuantile : Fraction of channels that need to have at least the given MinCorrelation value
%                     w.r.t. the channel under consideration. This allows to deal with channels or
%                     small groups of channels that measure the same noise source, e.g. if they are
%                     shorted. If many channels can be disconnected during an experiment and you
%                     have strong noise in the room, you might increase this fraction, but consider
%                     that this a) requires you to decrease the MinCorrelation appropriately and b)
%                     this can make the correlation measure more brittle. Reasonable range: 0.05 (rather
%                     lax) to 0.2 (very tolerant re disconnected/shorted channels).The default is
%                     0.1.
%
%   WindowLength    : Length of the windows (in seconds) for which correlation is computed; ideally
%                     short enough to reasonably capture periods of global artifacts (which are
%                     ignored), but not shorter (for statistical reasons). Default: 2.
% 
%   MaxBrokenTime : Maximum time (either in seconds or as fraction of the recording) during which a 
%                   retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                   (very lax). The default is 0.5.
%
%   LineNoiseAware : Whether the operation should be performed in a line-noise aware manner. If enabled,
%                    the correlation measure will not be affected by the presence or absence of line 
%                    noise (using a temporary notch filter). Default: true.
%
% Out:
%   Signal : data set with bad channels removed
%
% Notes:
%   This function requires the Signal Processing toolbox.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-07-06

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

if ~exist('min_corr','var') || isempty(min_corr) min_corr = 0.45; end
if ~exist('ignored_quantile','var') || isempty(ignored_quantile) ignored_quantile = 0.1; end
if ~exist('window_len','var') || isempty(window_len) window_len = 2; end
if ~exist('max_broken_time','var') || isempty(max_broken_time) max_broken_time = 0.5; end
if ~exist('linenoise_aware','var') || isempty(linenoise_aware) linenoise_aware = true; end


% flag channels
if max_broken_time > 0 && max_broken_time < 1  %#ok<*NODEF>
    max_broken_time = size(signal.data,2)*max_broken_time;
else
    max_broken_time = signal.srate*max_broken_time;
end

signal.data = double(signal.data);
[C,S] = size(signal.data);
window_len = window_len*signal.srate;
wnd = 0:window_len-1;
offsets = 1:window_len:S-window_len;
W = length(offsets);
retained = 1:(C-ceil(C*ignored_quantile));

% optionally ignore both 50 and 60 Hz spectral components...
if linenoise_aware
    Bwnd = design_kaiser(2*45/signal.srate,2*50/signal.srate,60,true);
    if signal.srate <= 130
        B = design_fir(length(Bwnd)-1,[2*[0 45 50 55]/signal.srate 1],[1 1 0 1 1],[],Bwnd);
    else
        B = design_fir(length(Bwnd)-1,[2*[0 45 50 55 60 65]/signal.srate 1],[1 1 0 1 0 1 1],[],Bwnd);
    end
    for c=signal.nbchan:-1:1
        X(:,c) = filtfilt_fast(B,1,signal.data(c,:)'); end
else
    X = signal.data';
end

% for each window, flag channels with too low correlation to any other channel (outside the
% ignored quantile)
flagged = zeros(C,W);
for o=1:W
    sortcc = sort(abs(corrcoef(X(offsets(o)+wnd,:))));
    flagged(:,o) = all(sortcc(retained,:) < min_corr);
end

% mark all channels for removal which have more flagged samples than the maximum number of
% ignored samples
removed_channels = sum(flagged,2)*window_len > max_broken_time;


%  signal.data = signal.data(~removed_channels,:);
 signal.removed_channels = removed_channels;

end


function B = design_fir(N,F,A,nfft,W)
% B = design_fir(N,F,A,nFFT,W)
% Design an FIR filter using the frequency-sampling method.
%
% The frequency response is interpolated cubically between the specified
% frequency points.
%
% In:
%   N : order of the filter
%
%   F : vector of frequencies at which amplitudes shall be defined
%       (starts with 0 and goes up to 1; try to avoid too 
%        sharp transitions)
%
%   A : vector of amplitudes, one value per specified frequency
%
%   nFFT : optionally number of FFT bins to use
%
%   W : optionally the window function to use (default: Hamming)
%
% Out:
%   B : designed filter kernel
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-14

% Copyright (C) Christian Kothe, SCCN, 2013, ckothe@ucsd.edu
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

if nargin < 4 || isempty(nfft)
    nfft = max(512,2^ceil(log(N)/log(2))); end
if nargin < 5
    W = 0.54 - 0.46*cos(2*pi*(0:N)/N); end

% calculate interpolated frequency response
F = interp1(round(F*nfft),A,(0:nfft),'pchip');

% set phase & transform into time domain
F = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B = real(ifft([F conj(F(end-1:-1:2))]));

% apply window to kernel
B = B(1:N+1).*W(:)';
end


function W = design_kaiser(lo,hi,atten,odd)
% Design a Kaiser window for a low-pass FIR filter
%
% In:
%   Lo : normalized lower frequency of transition band
%
%   Hi : normalized upper frequency of transition band
%
%   Attenuation : stop-band attenuation in dB (-20log10(ratio))
%
%   OddLength : whether the length shall be odd
%
% Out:
%   W : Designed window
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-08-17

% Copyright (C) Christian Kothe, SCCN, 2013, ckothe@ucsd.edu
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

% determine beta of the kaiser window
if atten < 21
    beta = 0;
elseif atten <= 50
    beta = 0.5842*(atten-21).^0.4 + 0.07886*(atten-21);
else
    beta = 0.1102*(atten-8.7);
end

% determine the number of points
N = round((atten-7.95)/(2*pi*2.285*(hi-lo)))+1;
if odd && ~mod(N,2)
    N = N+1; end

% design the window
W = window_func('kaiser',N,beta);
end



function w = window_func(name,m,param)
% Design a window for a given window function
%
% In:
%   Name : name of the window, can be any of the following:
%          'bartlett' : Bartlett window
%          'barthann' : Bartlett-Hann window
%          'blackman' : Blackman window
%          'blackmanharris' : Blackman-Harris window
%          'flattop'  : Flat-top window
%          'gauss'    : Gaussian window with parameter alpha (default: 2.5)
%          'hamming'  : Hamming window
%          'hann'     : Hann window
%          'kaiser'   : Kaiser window with parameter beta (default: 0.5)
%          'lanczos'  : Lanczos window
%          'nuttall'  : Blackman-Nuttall window
%          'rect'     : Rectangular window
%          'triang'   : Triangular window
%
%   N : number of points in the window
%
%   Param : window parameter (if any)
%
% Out:
%   W : designed window (column vector)
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2013-08-16

% Copyright (C) Christian Kothe, SCCN, 2011, ckothe@ucsd.edu
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

p = (0:(m-1))/(m-1);
switch name
    case 'bartlett'
        w = 1 - abs(((0:(m-1)) - (m-1)/2)/((m-1)/2));
    case {'barthann','barthannwin'}
        w = 0.62 - 0.48*abs(p-0.5) - 0.38*cos(2*pi*p);
    case 'blackman'
        w = 0.42-0.5*cos(2*pi*p) + 0.08*cos(4*pi*p);
    case 'blackmanharris'        
        w = 0.35875 - 0.48829*cos(2*pi*p) + 0.14128*cos(4*pi*p) - 0.01168*cos(6*pi*p);
    case {'bohman','bohmanwin'}
        w = (1-abs(p*2-1)).*cos(pi*abs(p*2-1)) + (1/pi)*sin(pi*abs(p*2-1));
    case {'flattop','flattopwin'}
        w = 0.2157 - 0.4163*cos(2*pi*p) + 0.2783*cos(4*pi*p) - 0.0837*cos(6*pi*p) + 0.0060*cos(8*pi*p);
    case {'gauss','gausswin'}
        if nargin < 3
            param = 2.5; end        
        w = exp(-0.5*(param*2*(p-0.5)).^2);
    case 'hamming'
        w = 0.54-0.46*cos(2*pi*p);
    case 'hann'
        w = 0.5-0.5*cos(2*pi*p);
    case 'kaiser'
        if nargin < 3
            param = 0.5; end
        w = besseli(0,param*sqrt(1-(2*p-1).^2))/besseli(0,param);
    case 'lanczos'
        w = sin(pi*(2*p-1))./(pi*(2*p-1)); w(isnan(w)) = 1;
    case {'nuttall','nuttallwin'}
        w = 0.3635819 - 0.4891775*cos(2*pi*p) + 0.1365995*cos(4*pi*p) - 0.0106411*cos(6*pi*p);
    case {'rect','rectwin'}
        w = ones(1,m);
    case 'triang'
        w = 1 - abs(((0:(m-1)) - (m-1)/2)/((m+1)/2));
    otherwise
        % fall back to the Signal Processing toolbox for unknown windows
        if nargin < 3
            w = window(name,m);
        else
            w = window(name,m,param);
        end
end

w = w(:);
end


function X = filtfilt_fast(varargin)
% Like filtfilt(), but faster when filter and signal are long (and A=1).
% Y = filtfilt_fast(B,A,X)
%
% Uses FFT convolution (needs fftfilt). The function is faster than filter when approx.
% length(B)>256 and size(X,Dim)>1024, otherwise slower (due size-testing overhead).
%
% Note:
%  Can also be called with four arguments, as Y = filtfilt_fast(N,F,A,X), in which case an Nth order
%  FIR filter is designed that has the desired frequency response A at normalized frequencies F; F
%  must be a vector of numbers increasing from 0 to 1.
%
% See also: 
%   filtfilt, filter
% 
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-07-14

% Copyright (C) Christian Kothe, SCCN, 2010, ckothe@ucsd.edu
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

if nargin == 3
    [B A X] = deal(varargin{:});
elseif nargin == 4
    [N F M X] = deal(varargin{:});
    B = design_fir(N,F,sqrt(M)); A = 1; % note: we use the sqrt() because we run forward and backward
else
    help filtfilt_fast;
    return;
end

if A == 1
    was_single = strcmp(class(X),'single');
    w = length(B); t = size(X,1);    
    % extrapolate
    X = double([bsxfun(@minus,2*X(1,:),X(1+mod(((w+1):-1:2)-1,t),:)); X; bsxfun(@minus,2*X(t,:),X(1+mod(((t-1):-1:(t-w))-1,t),:))]);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % filter, reverse
    X = filter_fast(B,A,X); X = X(length(X):-1:1,:);
    % remove extrapolated pieces
    X([1:w t+w+(1:w)],:) = [];
    if was_single
        X = single(X); end    
else    
    % fall back to filtfilt for the IIR case
    X = filtfilt(B,A,X);
end
end


function [X,Zf] = filter_fast(B,A,X,Zi,dim)
% Like filter(), but faster when both the filter and the signal are long.
% [Y,Zf] = filter_fast(B,A,X,Zi,Dim)
%
% Uses FFT convolution. The function is faster than filter when approx. length(B)>256 and
% size(X,Dim)>1024, otherwise slower (due size-testing overhead).
%
% See also:
%   filter, fftfilt
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-07-09
%
%                           contains fftfilt.m from Octave:
%                           Copyright (C) 1996-1997 John W. Eaton


% Copyright (C) Christian Kothe, SCCN, 2010, ckothe@ucsd.edu
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


if nargin <= 4
    dim = find(size(X)~=1,1); end
if nargin <= 3
    Zi = []; end

lenx = size(X,dim);
lenb = length(B);
if lenx == 0
    % empty X
    Zf = Zi;
elseif lenb < 256 || lenx<1024 || lenx <= lenb || lenx*lenb < 4000000 || ~isequal(A,1)
    % use the regular filter
    if nargout > 1
        [X,Zf] = filter(B,A,X,Zi,dim);
    else
        X = filter(B,A,X,Zi,dim);
    end
else
    was_single = strcmp(class(X),'single');
    % fftfilt can be used
    if isempty(Zi)
        % no initial conditions to take care of
        if nargout < 2
            % and no final ones
            X = unflip(oct_fftfilt(B,flip(double(X),dim)),dim);
        else
            % final conditions needed
            X = flip(X,dim);
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
            X = oct_fftfilt(B,double(X));
            X = unflip(X,dim);
        end
    else
        % initial conditions available
        X = flip(X,dim);
        % get a Zi-informed piece
        tmp = filter(B,1,X(1:length(B),:),Zi,1);
        if nargout > 1
            % also need final conditions
            [dummy,Zf] = filter(B,1,X(end-length(B)+1:end,:),Zi,1); %#ok<ASGLU>
        end
        X = oct_fftfilt(B,double(X));
        % incorporate the piece
        X(1:length(B),:) = tmp;
        X = unflip(X,dim);
    end
    if was_single
        X = single(X); end
end
end

function X = flip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = permute(X,order);
end
end

function X = unflip(X,dim)
if dim ~= 1
    order = 1:ndims(X);
    order = order([dim 1]);
    X = ipermute(X,order);
end
end


function y = oct_fftfilt(b, x, N)
% Copyright (C) 1996, 1997 John W. Eaton
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% -*- texinfo -*-
% @deftypefn {Function File} {} fftfilt (@var{b}, @var{x}, @var{n})
%
% With two arguments, @code{fftfilt} filters @var{x} with the FIR filter
% @var{b} using the FFT.
%
% Given the optional third argument, @var{n}, @code{fftfilt} uses the
% overlap-add method to filter @var{x} with @var{b} using an N-point FFT.
%
% If @var{x} is a matrix, filter each column of the matrix.
% @end deftypefn
%
% Author: Kurt Hornik <Kurt.Hornik@wu-wien.ac.at>
% Created: 3 September 1994
% Adapted-By: jwe

% If N is not specified explicitly, we do not use the overlap-add
% method at all because loops are really slow.  Otherwise, we only
% ensure that the number of points in the FFT is the smallest power
% of two larger than N and length(b).  This could result in length
% one blocks, but if the user knows better ...
transpose = (size(x,1) == 1);

if transpose
    x = x.'; end

[r_x,c_x] = size(x);
[r_b,c_b] = size(b);
if min([r_b, c_b]) ~= 1
    error('octave:fftfilt','fftfilt: b should be a vector'); end

l_b = r_b*c_b;
b = reshape(b,l_b,1);

if nargin == 2
    % Use FFT with the smallest power of 2 which is >= length (x) +
    % length (b) - 1 as number of points ...
    N = 2^(ceil(log(r_x+l_b-1)/log(2)));
    B = fft(b,N);
    y = ifft(fft(x,N).*B(:,ones(1,c_x)));
else
    % Use overlap-add method ...
    if ~isscalar(N)
        error ('octave:fftfilt','fftfilt: N has to be a scalar'); end
    N = 2^(ceil(log(max([N,l_b]))/log(2)));
    L = N - l_b + 1;
    B = fft(b, N);
    B = B(:,ones(c_x,1));
    R = ceil(r_x / L);
    y = zeros(r_x, c_x);
    for r = 1:R
        lo = (r - 1) * L + 1;
        hi = min(r * L, r_x);
        tmp = zeros(N, c_x);
        tmp(1:(hi-lo+1),:) = x(lo:hi,:);
        tmp = ifft(fft(tmp).*B);
        hi = min(lo+N-1, r_x);
        y(lo:hi,:) = y(lo:hi,:) + tmp(1:(hi-lo+1),:);
    end
end

y = y(1:r_x,:);
if transpose
    y = y.'; end

% Final cleanups: if both x and b are real respectively integer, y
% should also be
if isreal(b) && isreal(x)
    y = real(y); end
if ~any(b - round(b))
    idx = ~any(x - round(x));
    y(:,idx) = round(y(:,idx));
end
end