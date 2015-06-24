clear all 
close all
AFilter = ...  
[0.0008
    0.0005
   -0.0000
   -0.0006
   -0.0012
   -0.0015
   -0.0011
    0.0000
    0.0017
    0.0033
    0.0039
    0.0029
   -0.0000
   -0.0040
   -0.0076
   -0.0088
   -0.0063
    0.0000
    0.0084
    0.0157
    0.0181
    0.0129
   -0.0000
   -0.0173
   -0.0327
   -0.0387
   -0.0288
    0.0000
    0.0451
    0.0989
    0.1500
    0.1867
    0.2000
    0.1867
    0.1500
    0.0989
    0.0451
    0.0000
   -0.0288
   -0.0387
   -0.0327
   -0.0173
   -0.0000
    0.0129
    0.0181
    0.0157
    0.0084
    0.0000
   -0.0063
   -0.0088
   -0.0076
   -0.0040
   -0.0000
    0.0029
    0.0039
    0.0033
    0.0017
    0.0000
   -0.0011
   -0.0015
   -0.0012
   -0.0006
   -0.0000
    0.0005
    0.0008
];
M = 4;

[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
numBits = k*(2^8);

Txbits = double(rand(1,numBits) > 0.5);
[ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( Txbits, k );
[ TxSymbol ] = SymbolArray(TxSymbolIdx);
numSymbols = length(TxSymbol);
Freq = 1;
B = 0:Freq/(numSymbols-1):Freq;
Bandwidth = B(end);
Fn = B(end);

FsArray = [4, 8];
for i = 1:length(FsArray)
    Fs = FsArray(i);
    txUpSample = zeros(Fs*numSymbols,1);
    txUpSample(1:Fs:end) = TxSymbol;
    figure
    plot(20*log10(abs(fft(txUpSample))))
    pause
    
    n = length(txUpSample);
    ZeroIndex = find(txUpSample==0);
    Bin = Fs/n;
    temp = linspace(-Fs/2,Fs/2,n);
    test =temp(end)-temp(end-1);
    
    fOverSampled = 0:Fs/(n-1):Fs;
    FnOverSampled = fOverSampled(end);
    
    Beta = FnOverSampled/Fn;
    
    f = -Fs/2:Fs/(n-1):Fs/2;
    
    FilterRec = ones(Fs,1)/Fs;
    FilterAZeropad = conv(AFilter,txUpSample);
    FilterRecOutput = conv(FilterRec,txUpSample);
    FilterARect = conv(AFilter,FilterRecOutput);
    newNfft = 1024;
    HF = fft(FilterRec,newNfft);
    F = -Fs/2:Fs/(length(HF)-1):Fs/2;
    
    figure
    subplot(2,1,1);
    plot(F,20*log10(fftshift(abs(HF))));
    
    xlabel('Frequency');
    ylabel('Amplitude');
    title('Frequency Response of Rectangular filter');
    legend(['Beta = ', num2str(Fs)])
    
    Fsn = 8;
    FilterRec2 = ones(Fsn,1)/Fsn;
    HF = fft(FilterRec2,newNfft);
    F = -Fsn/2:Fsn/(length(HF)-1):Fsn/2;
    subplot(2,1,2);
    plot(F,20*log10(fftshift(abs(HF))));
    
    xlabel('Frequency');
    ylabel('Amplitude');
    title('Frequency Response of Rectangular filter');
    legend(['Beta = ', num2str(Fsn)])
    pause
    
    zeroPad = zeros(1024-length(AFilter),1);
    AFilterzeroPad = [AFilter; zeroPad];
    
    HFzeropad = fft(AFilterzeroPad);
    HF = fft(AFilter,newNfft);
    F = -Fs/2:Fs/(length(HF)-1):Fs/2;
    
    figure
%     plot(F,20*log10(fftshift(abs(HF))));
%     hold on
    subplot(2,1,1);  
    
%     t =(0:length(AFilter)-1)/Fs;
    n = length(AFilter);
    t = (-(n-1)/2:(n-1)/2);
    plot(t,AFilter);
    grid on
%     axis([-(n-1)/2 (n-1)/2 -0.05 0.2])
    IndexZ = find(AFilter == 0);
    ZeroCross = AFilter(IndexZ);
    t = t(IndexZ);
    hold on
    plot(t,ZeroCross,'ro');
    xlabel('T');
    ylabel('Amplitude');
    title('Impulse Response of Anologue filter');

    subplot(2,1,2);  
    plot(F,20*log10(fftshift(abs(HFzeropad))));
    
    xlabel('Frequency');
    ylabel('Amplitude');
    title('Frequency Response of Anologue filter');
    pause
    
    tAZeroPad = (0:length(FilterAZeropad)-1)/Fs;
    tARect = (0:length(FilterARect)-1)/Fs;
    
    fAZeroPad = 0:Fs/(length(FilterAZeropad)-1):Fs;
    fARect = 0:Fs/(length(FilterARect)-1):Fs;
    
    fAZeroPadShift = -Fs/2:Fs/(length(FilterAZeropad)-1):Fs/2;
    fARectShift = -Fs/2:Fs/(length(FilterARect)-1):Fs/2;
    
    B = 0:Fs/(length(txUpSample)-1):Fs;
    BandwidthFilter = B(end);
    
    figure 
    [mA, Index] =max(AFilter);
    subplot(2,1,1); 
    
    FilterAZeropad2 = FilterAZeropad;
    FilterAZeropad2(1:Index)=[];
    t = (0:length(FilterAZeropad2)-1)/Fs;
    plot(t,real(FilterAZeropad2))
%     plot(tAZeroPad,real(FilterAZeropad))
    hold on
    Txup =(0:length(txUpSample)-1)/Fs;
    plot(Txup,real(txUpSample),':r');
    axis([0 20 -1 1])
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title(['Analogue Filtered, Over Sampled Fs = ',num2str(Fs)]);

    subplot(2,1,2); 
%     plot(tARect,real(FilterARect))
%     AfilterRect2= FilterARect(Aindex+1:end);
%     t = (0:length(AfilterRect2)-1)/Fs;
    
    
    FilterARect2 = FilterARect;
    FilterARect2(1:Index)=[];
    t = (0:length(FilterARect2)-1)/Fs;
    plot(t,real(FilterARect2))
    hold on
    plot((0:length(FilterRecOutput)-1)/Fs,real(FilterRecOutput),':r');
    
%     samplePoint = AfilterRect2(1:4:numBits);
    t = (0:length(TxSymbol)-1);
    stem(t,real(TxSymbol)/Fs,'g')
    axis([0 20 -1/Fs 1/Fs])
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('Rectangular Filter+Analogue Filtered');
    pause
    FilterAZeropadB = 20*log10(abs(fft(FilterAZeropad)));
    FilterAZeropadBShift = 20*log10(abs(fftshift(fft(FilterAZeropad))));
    
    FilterARectdB = 20*log10(abs(fft(FilterARect)));
    FilterARectdBShift = 20*log10(abs(fftshift(fft(FilterARect))));
   
%     figure
%     subplot(2,1,1); plot(fft(FilterAZeropad)/fft(txUpSample))
%     xlabel('Frequency');
%     ylabel('dB');
%     title(['Analogue Filtered, Over Sampled Fs = ',num2str(Fs)]);
%     subplot(2,1,2); plot(fARect,FilterARectdB)
%     xlabel('Frequency');
%     ylabel('dB');
%     title('Rectangular Filter+Analogue Filtered');
%     pause
%     
    figure
    subplot(2,1,1); plot(fAZeroPad,FilterAZeropadB)
    xlabel('Frequency');
    ylabel('dB');
    title(['Analogue Filtered, Over Sampled Fs = ',num2str(Fs)]);
    subplot(2,1,2); plot(fARect,FilterARectdB)
    xlabel('Frequency');
    ylabel('dB');
    title('Rectangular Filter+Analogue Filtered');
    pause

    figure
    subplot(2,1,1); plot(fAZeroPadShift,FilterAZeropadBShift)
    xlabel('Frequency');
    ylabel('dB');
    title(['Time Domain Transmitted Bits, Bandwidth = ', num2str(BandwidthFilter),'Hz']);
    subplot(2,1,2); plot(fARectShift,FilterARectdBShift)
    xlabel('Frequency');
    ylabel('dB');
    title(['Over Sampled Fs = ',num2str(Fs), ' Bandwidth = ', num2str(BandwidthFilter),'Hz']);
    pause
    
    val = 65;
    [eyeDRect] = eyeDiagram( val , FilterAZeropad );
    [eyeDZeroPad] = eyeDiagram( val , FilterRecOutput );
    
    figure
    plot(eyeDZeroPad,'b')
    title('Eye Diagram Zero Pad Filtered Signal');
    pause
    
    figure
    plot(eyeDRect,'b')
    title('Eye Diagram Rect Filtered Signal');
    pause
    
%     NfftArray = [2^7 2^8 2^9 2^10];
%     
%     for j = 1:2
%         if j ==1
%             x = txUpSample;
%             xTitle = 'Over Sampled Singal   ';
%             tempPlot = txUpSampledB;
%             tempPlot2 = txUpSampledB;
%             tempF = 0:Fs/(n-1):Fs;
%             tempPD = txUpSample;
%         elseif j ==2
%             x = FilterRecOutput;
%             xTitle = 'Rectangular Filtered Singal   ';
%             tempPlot = FilterRecOutputdB;
%             tempPlot2 = FilterRecOutputdBShift;
%             tempF = 0:Fs/(length(FilterRecOutputdB)-1):Fs;
%             tempPD = FilterRecOutput;
%         end
%         for ii = 1:length(NfftArray)
%             
%             Nfft = NfftArray(ii);
%             [ win ] = Hanning( Nfft );
%             [Pxx,F] = pwelch(x,win,[],Nfft,Fs);
%             B =F(end);
%             N = length(Pxx);
%             
%             Fpd = -Fs/2:Fs/(length(tempPD)-1):Fs/2;
%             
%             
%             figure
%             y = fft(tempPD, length(tempPD));
%             y0 = fftshift(y);
%             
%             power = y0.*conj(y0)/length(tempPD);
%             subplot(4,1,1); plot(Fpd,10*log10(abs(power)))
%             title([xTitle,'Periodogram']);
%             
%             subplot(4,1,2); plot(win)
%             axis([0 length(win) 0 1])
%             title([xTitle,'Hanning Window']);
%             
% 
%             subplot(4,1,3); plot(F,10*log10(Pxx),'r')
%             hold on
%             plot(tempF,tempPlot)
%             hold off
%             title([xTitle,'PSD Estimate Fs = ', num2str(Fs),'  Overlap percentage = 50%','  NFFT = ',num2str(Nfft)]);
%             xlabel('Frequency');
%             ylabel('dB');
%             
%             F = -Fs/2:Fs/(N-1):Fs/2;
%             subplot(4,1,4); plot(F,10*log10(fftshift(Pxx)),'r')
%             hold on
%             plot(-Fs/2:Fs/(length(tempPlot)-1):Fs/2,tempPlot2)
%             hold off
%             title([xTitle,'PSD Estimate Fs = ', num2str(Fs),'  Overlap percentage = 50%','  NFFT = ',num2str(Nfft)]);
%             xlabel('Frequency');
%             ylabel('dB');
%             pause
%             
%             figure
%             plot(Fpd,10*log10(abs(power)))
%             hold on
%             plot(F,10*log10(fftshift(Pxx)),'r')
%             hold off
%             pause
%             
%             
%         end
%     end

end