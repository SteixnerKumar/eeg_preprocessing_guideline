%% General Fourier Transform
% 
%   [XFreqRange, YAmplitude] = sk_dofft(Data, Fs, wType)
%
% Output:
%   XFreqRange - Frequency Vector 
%   YAmplitude - Amplitude Vector
%
% Input:
%   Data  - Temporal data vector
%   Fs    - Sampling Frequency
%   wType - Window Type / Algorithm 
%           0 - no Window 
%           2 - Hanning
%           3 - Flattop  
%           4 - Welch algorithm (personal ppreference)
%           5 - Manual algorithm incl. scaling
%           6 - Manual FFT without scaling
%           7 - Inverse FFT agorithm
%           8 - Invserse FFT (unscaled coeff required)
%
% Creator: Saurabh Kumar
%

function [XFreqRange, YAmplitude]= sk_dofft(Data, Fs, wType)
XFreqRange=[];
YAmplitude=[];
Window=[];

if Data==0
    XFreqRange=0;
    YAmplitude=0;
    return;
end;    

Window=WindowType(wType,length(Data));

if ~isempty(Window),Data=Data.*Window'*length(Data)/sum(Window);end

if wType<4,[XFreqRange, YAmplitude]=posFFT(Data,Fs);end

if wType==5,[XFreqRange, YAmplitude]=manFFT(Data,Fs);end

if wType==6,[XFreqRange, YAmplitude]=ManFFT4inv(Data, Fs);end

if wType==7,XFreqRange=invFFT(Data,Fs);end

if wType==8,XFreqRange=invFFTonly(Data);end

if wType==4         % Welch Algorithm
    
    Resolution=floor(length(Data)/Fs);
    if Resolution==0
        L=length(Data); % window length
        nfft=L;
        over=0;
    else
        Resolution=0.1;       
        tmp=Resolution/(Fs/length(Data));   %Ratio Resolution:Sampling Resolution
        if tmp<1,tmp=1;end;            
        L=floor(length(Data)/tmp);      %length of segment
        
        tmp= 2^(nextpow2(L)-2);
        nfft= ceil(L/tmp)*tmp;          % length of fft
        
        over=floor(0.8*L);              % length overlapping
    end;
    
    Window=WindowType(3,L); % flattop
    
    S1=sum(Window);
    S2=sum(Window.^2);
    
    try
        [Pxx,XFreqRange]=pwelch(Data,Window,over,nfft,Fs);
        dF=XFreqRange(2);

        YAmplitude=  (2/sqrt(2)) * sqrt(Pxx/2) * sqrt(nfft*S2./(S1.^2)*dF);  %Convert PSD to g_rms
    catch
        disp('Signal Toolbox Licence Error. No window function is used...');
        [XFreqRange, YAmplitude]=posFFT(Data,Fs);
    end;    
    
end;


Cut=length(Data)/2;
Cut=[];
if ~isempty(Cut)
    XFreqRange=XFreqRange(1:Cut);
    YAmplitude=single(YAmplitude(1:Cut));
end;
%% FFT
function [XFreqRange, YAmplitude]=posFFT(x, Fs)

 %% FFT x-Data, Fs-Sampling Freq
N=2^(nextpow2(length(x)));
%N=length(x);

k=[0:N-1]/N;
XFreqRange=single(k*Fs);
YAmplitude= single(2/sqrt(2)* fft(x,N)/length(x));    %rms



 %% cutoff negative values
cutoff=ceil(N/2);
XFreqRange=XFreqRange(1:cutoff);
YAmplitude=YAmplitude(1:cutoff);
%YAmplitude(1)=0;

%% MANUAL FFT
function [XFreqRange, YAmplitude]=manFFT(x, Fs)

N=length(x);


for k=1:N
    for j=1:N
        z(k,j)=( x(j)*exp(-2*pi*i*(j-1)*(k-1)/N) );
    end;    
    X(k)=sum(z(k,:));
end;

z=[];

k=[0:N-1]/N;
XFreqRange=single(k*Fs);
YAmplitude= single(sqrt(2)*abs(X)/N);    %rms




 %% cutoff negative values
cutoff=ceil(N/2);
XFreqRange=XFreqRange(1:cutoff);
YAmplitude=YAmplitude(1:cutoff);
YAmplitude(1)=0;

%% Windowing
function Window=WindowType(wType,N)

try

    Window=[];    
    switch wType
        case 2
            Window = hann(N).'; %create a hanning window vector
        case 3
            Window=flattopwin(N).';
    end

catch
    Window=[];
end;    

%% inverse FFT
function x=invFFT(x, Fs)

% manual FFT
N=length(x);
for k=1:N
    for j=1:N
        z(k,j)=( x(j)*exp(-2*pi*i*(j-1)*(k-1)/N) );
    end;    
    X(k)=sum(z(k,:))/N;
end;

z=[];

%iFFT
N=length(X);
for j=1:N
    for k=1:N
        z(j,k)=( X(k)*exp(-(-2*pi*i*(j-1)*(k-1)/N)));
    end;    
    x(j)=sum(z(j,:));
end;

z=[];

%% manual FFT for inverse FFT
function [XFreqRange, YAmplitude]=ManFFT4inv(x, Fs)

%% FFT Transform without scaling
N=length(x);
k=[0:N-1]/N;
XFreqRange= (k*Fs);
YAmplitude= fft(x,N);
YAmplitude(1) = 0;

%% inverse FFT
function x=invFFTonly(YAmplitude)

n=length(YAmplitude);
x=real(ifft(YAmplitude(:,1),n));

