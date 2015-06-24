close all
clear all
% Up Sample
% ---------
Fs = 8;
No = (1);    % Noise Spectral Density
SingalSize = 2^12;
Noise = (((randn(1,SingalSize))+1i*(randn(1,SingalSize)))/(sqrt(2))); 
Noise=Noise';
x=Noise;
% x = ones(SingalSize,1);
% x =x/length(x);
% varx = var(x);  % unit variance 
% f=1;                            % Input Signal Frequency
% fs=2000;                        % Sampling Frequency
% t=0:f/fs:1-f/fs;              % Time Samples
% x=(sin(2*pi*t)/length(t)).';                  % Generate Sine Wave  

M = length(x);
NfftArray = [2^7 2^8];

for i = 1:length(NfftArray)
    Nfft = NfftArray(i);
    
    % Hanning Window
    % --------------
    [ win ] = Hanning( Nfft );
    Ls = length(win);
    
    noverlap = fix(0.5.*Ls);
    
    k = (M-noverlap)./(Ls-noverlap);
    k = fix(k);
    
    % Define Overlap
    % --------------
    LminusOverlap = Ls-noverlap;
    xStart = 1:LminusOverlap:k*LminusOverlap;
    xEnd   = xStart+Ls-1;
    
    Sxx = zeros(Nfft,1);                    % Power Spectrum Density accumulator
    U = win'*win;                           % Compensates for the power of the window.
    for ii = 1:k,
        xk = x(xStart(ii):xEnd(ii));
        xw = xk.*win;
        Xx = fft(xw,Nfft);
        P = Xx.*conj(Xx)/U;                 % Auto spectrum.
        Sxx = Sxx + P;                      % Sum scaled periodograms
    end
    Sxx = Sxx./k;                           % Average the sum of the periodogram
    Pxx = Sxx./Fs; % Scale by the sampling frequency to obtain the psd

    
    figure
    
    HF = fft(x,length(x)*256);
    F = -Fs/2:Fs/(length(HF)-1):Fs/2;
    subplot(2,1,1); plot(F,10*log10(fftshift(abs(HF))));
    xlabel('Frequency');
    ylabel('Amplitude dB');
     title('FFT of original signal')
    F = -Fs/2:Fs/(length(Pxx)-1):Fs/2;
    subplot(2,1,2); plot(F,10*log10(abs(fftshift(Pxx))))
    xlabel('Frequency');
    ylabel('Amplitude dB');
    title(['PSD Estimate Fs = ', num2str(Fs),'  Overlap percentage = 50%','  NFFT = ',num2str(Nfft)]);
    hold off
    MaxValue = max(10*log10(abs(fftshift(Pxx))))
    pause
end