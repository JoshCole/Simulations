close all
clear all
Fs= 4;

[ FilterRRC ] = RootRaiseCosine( 0.22, 32, Fs );
N = length(FilterRRC);
T = (0:N-1)/Fs;
freq_resUp = Fs/N;
F = freq_resUp*(0:N-1);
 
figure
subplot(2,1,1); plot(T,FilterRRC)
xlabel('Time (Sec)');
ylabel('Amplitude');

FilterRRCdB = 10*log10(abs(fft(FilterRRC)));
subplot(2,1,2); plot(F,FilterRRCdB)
xlabel('Frequency');
ylabel('Amplitude');

figure
h=conv(FilterRRC,FilterRRC);
subplot(2,1,1)
plot(h)

subplot(2,1,2)
plot(h(1:Fs:end))
