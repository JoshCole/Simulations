clear all 
close all
Fs= 8;
winRec = ones(1,Fs);
n = length(winRec);
t = (0:n-1)/Fs;
freq_res = Fs/n;
f = freq_res*(0:n-1);

figure
subplot(2,1,1); plot(t,winRec)
xlabel('Frequency');
ylabel('Amplitude');

winRecdB = 10*log10(abs(fftshift(fft(winRec))));
subplot(2,1,2); plot(f,winRecdB)
xlabel('Frequency');
ylabel('Amplitude');

    

[ FilterRRC ] = RootRaiseCosine( 0.22, 12, Fs );
N = length(FilterRRC);
T = (0:N-1)/Fs;
freq_resUp = Fs/N;
F = freq_resUp*(0:N-1);

figure
subplot(2,1,1); plot(T,FilterRRC)
xlabel('Time (Sec)');
ylabel('Amplitude');

FilterRRCdB = 10*log10(abs(fftshift(FilterRRC)));
subplot(2,1,2); plot(F,FilterRRCdB)
xlabel('Frequency');
ylabel('Amplitude');