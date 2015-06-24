close all
clear all
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

FsArray = [2];
for i = 1:length(FsArray)
    Fs = FsArray(i);
    txUpSample = zeros(Fs*numSymbols,1);
    txUpSample(1:Fs:end) = TxSymbol;
    
    % RRC Filter
    [ FilterRRC1 ] = (RootRaiseCosine( 0.22, Fs )/Fs)';
    FilterOutputRRC = conv(FilterRRC1,txUpSample);
    N = length(FilterOutputRRC);
    
    % Frequency Response
    
    Nfft = 1024;
    Nh = length(FilterRRC1);
    zeroPad = zeros(Nfft -Nh,1);
     
    h = [FilterRRC1;zeroPad];
    H = fft(h);
    f = -Fs/2:Fs/(Nfft-1):Fs/2;
    HShift = fftshift(abs(H));
    
    figure
    subplot(3,1,1);
    n = length(FilterRRC1);
    t = (-(n-1)/2:(n-1)/2)/Fs; 
    plot(t,FilterRRC1);
    hold on
    grid on
    IndexZ = find(FilterRRC1 == 0);
    ZeroCross = FilterRRC1(IndexZ);
    t = t(IndexZ);
    plot(t,ZeroCross,'ro');
    xlabel('T');
    ylabel('Amplitude');
    title('Impulse Response of RRC Filter');
    
    subplot(3,1,2);
    plot(f,20*log10(HShift));
    xlabel('Frequency');
    ylabel('Amplitude');
    title('Frequency Response of RRC filter');
    legend('Pulse Shaping Filter')
    
    % Filter Peak
   [Peak, index] =max(FilterRRC1);
    index = index-1;
    % Pulse Shaping
    FilterOutputRRC2 = FilterOutputRRC;
    FilterOutputRRC2(1:index)=[]; 
    t = (0:length(FilterOutputRRC2)-1)/Fs;
    
    subplot(3,1,3);plot(t,real(FilterOutputRRC2));
    hold on
    t = (0:length(TxSymbol)-1);
    stem(t,real(TxSymbol)/Fs,'r')
    axis([0 20 -1.25/Fs 1.25/Fs])
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title(['Over Sampled Fs = ',num2str(Fs)]);
    legend('RRC Pulse Shaped Signal','Sampled Signal')
    pause 
        % Applying Analogue Filter
    Filteredsig = conv(AFilter,FilterOutputRRC);
    
    
    figure
    subplot(2,1,1)
    X = fft(FilterOutputRRC);
    f = -Fs/2:Fs/(length(X)-1):Fs/2;
    XShift = fftshift(abs(X));
    plot(f,20*log10(XShift));
    axis([-4 4 -100 50])
    subplot(2,1,2)
    X = fft(Filteredsig);
    f = -Fs/2:Fs/(length(X)-1):Fs/2;
    XShift = fftshift(abs(X));
    plot(f,20*log10(XShift));
    axis([-4 4 -200 50])
    pause
    
    val = 65;
    [eyeD] = eyeDiagram( val , FilterOutputRRC );
    
    figure
    plot(eyeD,'b')
    title('Eye Diagram Rect Filtered Signal');
    pause
    
   [eyeD] = eyeDiagram( val , Filteredsig );
    
    figure
    plot(eyeD,'b')
    title('Eye Diagram Rect Filtered Signal');
    pause
    
   FilterOutputRRC1 = conv(FilterRRC1,FilterOutputRRC);
   FilterOutputRRC4 = conv(FilterRRC1,Filteredsig);
   
   % Filter Peak
   [Peak2, index2] = max(AFilter);
  index2 = index2-1;
  
    figure
    index= 2*index;
    FilterOutputRRC3 = FilterOutputRRC1;
    FilterOutputRRC3(1:index)=[]; 
    t = (0:length(FilterOutputRRC3)-1)/Fs;
    plot(t,real(FilterOutputRRC3));
    hold on
    t = (0:length(TxSymbol)-1);
    stem(t,real(TxSymbol)/Fs,'r')
    
    index2 = index2+index;
    FilterOutputRRC5 = FilterOutputRRC4;
    FilterOutputRRC5(1:index2)=[]; 
    t = (0:length(FilterOutputRRC5)-1)/Fs;
    plot(t,real(FilterOutputRRC5),':g');
    
    axis([0 20 -1.25/Fs 1.25/Fs])
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title(['Over Sampled Fs = ',num2str(Fs)]);
    legend('RRC Pulse Shaped Signal','Sampled Signal', 'Bandlimited Singal')
    pause 
   
   [eyeD] = eyeDiagram( 2*val , FilterOutputRRC1 );
    figure
    plot(eyeD,'b')
    title('Eye Diagram Rect Filtered Signal');
    pause
    
   [eyeD] = eyeDiagram( 2*val , FilterOutputRRC4 );
    
    figure
    plot(eyeD,'b')
    title('Eye Diagram Rect Filtered Signal');
    pause
   
    
end