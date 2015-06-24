close all
clear all

Fs = 100;
t = 0:1/Fs:10;
y = sin(2*pi*15*t) + sin(2*pi*30*t);
nfft = 512;
Y = fft(y,nfft);
f = Fs*(0:nfft-1)/nfft;
Power = Y.*conj(Y)/nfft;
figure
plot(f,Power)
title('Periodogram')
figure
plot(t,y)
%%
rate = 8;
numBits = 256;
t = 0:1/rate:numBits-1/rate;
h = ones(1,rate);
filterNomalisation = 1/sqrt(rate);
ovTxbit = zeros(rate*numBits,1);
txbit = (rand(numBits,1)<0.5)*2-1;
ovTxbit(1:rate:end) = txbit;
fil = filterNomalisation*conv(h,ovTxbit);
NFFT = 512;
Y = fft(fil,NFFT);
f = rate*(0:NFFT-1)/NFFT;
Power = Y.*conj(Y)/NFFT;

figure
plot(t,ovTxbit)

figure
plot(f,Power)
title('Periodogram')
%%
close all
clear all
% Up Sample
% ---------
Fs = 8;
numBits = 2^16;
ovTxbit = zeros(Fs*numBits,1);
txbit = (rand(numBits,1)<0.5)*2-1;
ovTxbit(1:Fs:end) = txbit;
x= ovTxbit;
M = length(x);
NfftArray = [2^7 2^8 2^9 2^10];

for i = 1:length(NfftArray)
    Nfft = NfftArray(i);
    % Hanning Window
    % --------------
    [ win ] = Hanning( Nfft );
    Ls = length(win);
    
    noverlap = fix(0.5.*Ls);
    Nfft = Ls;
    
    k = (M-noverlap)./(Ls-noverlap);
    k = fix(k);
    
    % Define Overlap
    % --------------
    LminusOverlap = Ls-noverlap;
    xStart = 1:LminusOverlap:k*LminusOverlap;
    xEnd   = xStart+Ls-1;
    
    Sxx = zeros(Nfft,1);                    % Power Spectrum Density accumulator
    
    for ii = 1:k,
        xk = x(xStart(ii):xEnd(ii));
        xw = xk.*win;
        U = win'*win;                       % Compensates for the power of the window.
        Xx = fft(xw,Nfft);
        P = Xx.*conj(Xx)/U;                 % Auto spectrum.
        Sxx = Sxx + P;                      % Sum scaled periodograms
    end
    Sxx = Sxx./k;                           % Average the sum of the periodograms
    
% Generate the one-sided spectrum [Power] if so wanted    
%    if rem(Nfft,2),
%       select = 1:(Nfft+1)/2;  % ODD
%       Sxx_unscaled = Sxx(select,:); % Take only [0,pi] or [0,pi)
%       Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end,:)];  % Only DC is a unique point and doesn't get doubled
%    else
%       select = 1:Nfft/2+1;    % EVEN
%       Sxx_unscaled = Sxx(select,:); % Take only [0,pi] or [0,pi)
%       Sxx = [Sxx_unscaled(1,:); 2*Sxx_unscaled(2:end-1,:); Sxx_unscaled(end,:)]; % Don't double unique Nyquist point
%    end
%    F = F(select);

    Pxx = Sxx./Fs; % Scale by the sampling frequency to obtain the psd
    freq_res = Fs/Nfft;
    F = freq_res*(0:Nfft-1);
    
    figure
    subplot(2,1,1); plot(win)
    title('Hanning Window');
    subplot(2,1,2); plot(F,Pxx)
    title(['PSD Estimate Fs = ', num2str(Fs),'  Overlap percentage = 50%','  NFFT = ',num2str(Nfft)]);
    pause
end
%%
close all
clear all
Fs = 8;
numBits = 2^16;
t  = (0:1/Fs:numBits-1/Fs).'; 
ovTxbit = zeros(Fs*numBits,1);
txbit = (rand(numBits,1)<0.5)*2-1;
ovTxbit(1:Fs:end) = txbit;
x= ovTxbit;
Nx = length(x);
NfftArray = [2^7 2^8 2^9 2^10];

for i = 1:length(NfftArray)
    Nfft = NfftArray(i);
    [ w ] = Hanning( Nfft );
    [Pxx,F] = pwelch(x,w,[],Nfft,Fs,'twosided');
    
    figure
    subplot(2,1,1); plot(w)
    title('Hanning Window');
    subplot(2,1,2); plot(F,Pxx)
    title(['PSD Estimate Fs = ', num2str(Fs),'  Overlap percentage = 50%','  NFFT = ',num2str(Nfft)]);
    pause
end