clear all 
close all
M = 4;

[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
numBits = k*(2^8);

Txbits = double(rand(1,numBits) > 0.5);
[ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( Txbits, k );
[ TxSymbol ] = SymbolArray(TxSymbolIdx);
numSymbols = length(TxSymbol);
n = numSymbols;
B = 1;
t = 0:n-1;
freq_res = B/n;
f = freq_res*(0:n-1);
Bandwidth = f(end)-f(1);

FsArray = ([2, 4, 8])*B;
for i = 1:length(FsArray)
    Fs = FsArray(i);
    txUpSample = zeros(Fs*numSymbols,1);
    txUpSample(1:Fs:end) = TxSymbol;

    N = length(txUpSample);
    T = (0:N-1)/Fs;
    freq_resUp = Fs/N;
    F = freq_resUp*(0:N-1);
    BandwidthUpSample = F(end);

    figure
    subplot(3,1,1); plot(t,real(TxSymbol))
    axis([0 20 -1 1])
    xlabel('Time (Sec)');
    ylabel('Amplitude');
    title(['Time Domain Transmitted Bits, Bandwidth = ', num2str(Bandwidth),'Hz']);
    
    subplot(3,1,2); plot(T,real(txUpSample))
    axis([0 20 -1 1])
    xlabel('Time (Sec)');
    ylabel('Amplitude');
    title(['Time Domain Over Sampled Fs = ',num2str(Fs), ' Bandwidth = ', num2str(BandwidthUpSample),'Hz']);
    

%     IndexZero = find(txUpSample == 0);
%     IndexNZero = find(txUpSample ~= 0);
%     temp = txUpSample;
%     temp(IndexZero) = [];
%     Ttemp = T;
%     Ttemp(IndexZero) = [];
    subplot(3,1,3); stem(t,real(TxSymbol),'r')
    hold on
    plot(t,real(TxSymbol),'--')
    hold off
    axis([0 20 -1 1]);
    legend('Sampled Signal (Discrete)', 'Generated Signal');
    xlabel('Time (Sec)');
    ylabel('Amplitude');
    title(['Time Domain Over Sampled Fs = ',num2str(Fs), ' Bandwidth = ', num2str(BandwidthUpSample),'Hz']);
    pause

    figure
    subplot(2,1,1); plot(f,10*log10(abs(fft(TxSymbol))))
    % axis([0 20 -1 1])
    xlabel('Frequency (Hz)');
    ylabel('dB');
    title(['Frequency Domain Transmitted Bits, Bandwidth = ', num2str(Bandwidth),'Hz']);

    subplot(2,1,2); plot(F,10*log10(abs(fft(txUpSample))))
    % axis([0 20 -1 1])
    xlabel('Frequency (Hz)');
    ylabel('dB');
    title(['Frequency Domain Over Sampled Fs = ',num2str(Fs), ' Bandwidth = ', num2str(BandwidthUpSample),'Hz']);
    pause
end
close all