close all
clear all
% Create Reference Signal            
G = [1 0 0 0 0 1 1];
GF = G(2:end); 
N = length(G)-1; 
M = 2^N-1;
ProcessGainMsequence = 10*log10(M);
Reg = zeros(1,N); 
x = zeros(M,1);
 
% Initalise Reg 
Reg(end)=1;   
 
for n = 1:M
    buffer = mod(sum(GF.*Reg),2);      % Caluculate New Bit
    Reg = [buffer,Reg(1:N-1)];         % Shift Register   
    x(n) = Reg(N);
end
ReferenceSignal = (x.*2-1).';
MeanRef = mean(ReferenceSignal);
VarRef = var(ReferenceSignal);
CyclicPrefix = ReferenceSignal(M-14:M);
Ncp = length(CyclicPrefix);

EsNo_dB(:,1) = linspace(50, 0, 20);      % Noise density dB
% EsNo_dB(:,1) = linspace(0, 10, 20);      % Noise density dB
EsNo_linArray = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_linArray);         % Numb of Iterations
  
M = 4;
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
[ LLRBitArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M, k );
Ntxp = 10^3;

for num = 1:NumEsNo
TxPacket = double(rand(1,Ntxp) > 0.5);
[ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( TxPacket, k );
[ TxSymbol ] = SymbolArray(TxSymbolIdx);    % Gets Symbols to be Transmitted for Symbol Array
Nsymb = length(TxSymbol)/2;
TxSymbolRe = reshape(TxSymbol,2,[]);
Tx = [TxSymbolRe(1,:) CyclicPrefix ReferenceSignal TxSymbolRe(2,:)];
 
% Up Sampled Signal
% -----------------
NyquistRate = 4;
[ TxOVS, NTx] = OverSample( Tx, NyquistRate );
[ TxRefOVS, NTxRef] = OverSample( ReferenceSignal, NyquistRate );      % Varaiance Reduced by a factor of the Nyquist Rate
[ FilterRRC ] = (RootRaiseCosine( 0.22, NyquistRate )); %
TxUnFiltered = Tx;
% 1st RRC Filter
% ----------
TxRRC1 = conv(FilterRRC*NyquistRate,TxOVS);

% Channel Delay
% -------------

Delay = 4*NyquistRate;
SymbolyDelay = Delay/NyquistRate;
h = zeros(1,SymbolyDelay*NyquistRate);
h(randi([1 SymbolyDelay*NyquistRate],1)) = 1; % Delay Unit Power
TxDelay = conv(h,TxRRC1);  
% Convole transmitted signal will channel

    EsNo_lin = EsNo_linArray(num);
    [ Rx, AWGN ,No, NoiseVar, rayFade ] = RayleighFade( TxDelay, Es, EsNo_lin, NyquistRate, 'Rayleigh' );
    
%     f = -NyquistRate/2:NyquistRate/(length(Tx)-1):NyquistRate/2;
%     f2 = -NyquistRate/2:NyquistRate/(length(TxUnFiltered)-1):NyquistRate/2;
%  
%     figure
%     plot(f,20*log10(fftshift(abs(fft(Tx))))); hold on; plot(f,20*log10(fftshift(abs(fft(AWGN)))),'r'); ...
%     plot(f2,20*log10(fftshift(abs(fft(TxUnFiltered)))),'g'); 
%     grid on;
%     legend('Rx Signal', 'AWGN', 'Tx Unfiltered');
    Nrx = length(Rx);
    NoiseVarArray(num) = NoiseVar;

    % Channel and Reciever Noise
    varAWGN = var(AWGN); % Approximately the Nyquist Rate * spectal Noise Density


    % 2nd RRC Filter
    % --------------
    RxRRC2 = conv(FilterRRC,Rx);
    Rxawgn = conv(FilterRRC,AWGN);
    
    varRxawgn = var(Rxawgn); % Approximately Equal to the Spectral Density

%     f = -NyquistRate/2:NyquistRate/(length(RxRRC2)-1):NyquistRate/2;

%     figure
%     subplot(3,1,1)
%     plot(abs(Rx));hold on; plot(abs(Rxawgn),'r')
%     legend('Rx Signal', 'AWGN')
%     subplot(3,1,2)
%     plot(f,20*log10(fftshift(abs(fft(Rx))))); hold on; plot(f,20*log10(fftshift(abs(fft(Rxawgn)))),'r');
%     legend('Rx Signal','AWGN')
    RxAWGN = RxRRC2 + Rxawgn;
%     subplot(3,1,3)
%     plot(f,20*log10(fftshift(abs(fft(Rx)))));hold on; plot(f,20*log10(fftshift(abs(fft(TxOriginal)))),'r');
%     legend('Rx Signal + AWGN', 'Tx')

 
    Nref = length(ReferenceSignal);

    MidPoint = round(Nrx/2); 
    i = MidPoint -round((Nref*NyquistRate)/2);
    j = i+(Nref*NyquistRate)-1;

    % Channel Estimation
    [ TxRefFreq, TxRefTime ] = CyclicCorr( TxRefOVS, TxRefOVS );
    LocalReplica = abs(TxRefTime)/Nref;
    
    TestSig = RxAWGN(i:j);  
    [ ChannelEstimateFreq, ChannelEstimateTime ] = CyclicCorr( TestSig, TxRefOVS );
    ChannelEstimateNorm = abs(ChannelEstimateTime)/Nref;
    [Peak, Index] = max(ChannelEstimateNorm);
 
%     figure  
%     plot(ChannelEstimateNorm); hold on; plot(LocalReplica,'r');

    i = i + Index-1;
    j = i+(Nref*NyquistRate)-1;
    if j<=Nrx
        % Sync Channel Estimation 
        TestSig = RxAWGN(i:j);
        [ SyncFreq, SyncTime ] = CyclicCorr( TestSig, TxRefOVS );
        SyncNorm = SyncTime/Nref;
        SyncNormABS = abs(ChannelEstimateTime)/Nref;
        [Peak, Index] = max(SyncNormABS);

%         figure
%         plot(SyncNormABS); hold on; plot(LocalReplica,'r');
        RefNorm = LocalReplica(1:NyquistRate:end);
        SyncNormDec = SyncNorm(1:NyquistRate:end);

%         figure  
%         plot(abs(SyncNormDec)); hold on; plot(RefNorm,'r');

        % Gain and Phase Alignment
        AlignmentTerm = 1/SyncNormDec(1);

        % Start and End Position 
        Sp = (i-(NyquistRate*Ncp));
        Ep = j;

        % Tx Original and Tx phas and Gain Aligned
        RxOriginal = RxAWGN;
        Rx_Aligned = RxAWGN*AlignmentTerm;  
%         figure
%         plot(abs(Rx_Aligned(i:NyquistRate:Ep)),'--g'); hold on; plot(abs(RxOriginal(i:NyquistRate:Ep)),'r');%plot(abs(AlignmentTerm*TxOriginal(i:NyquistRate:Ep)),':m');
%         legend('Rx Aligned', 'Rx Original')
        meanRx = mean(abs(Rx_Aligned(i:NyquistRate:Ep)));
 

        % Recieved Reference Sig 
        RxRef = Rx_Aligned(Sp:NyquistRate:Ep);
        RxSymbol = [Rx_Aligned((Sp-NyquistRate):-NyquistRate:Sp-NyquistRate*Nsymb) Rx_Aligned(Ep+1:NyquistRate:Ep+NyquistRate*Nsymb)];
        % Recieved Channel Estimated Decimated
        RxChan = ReferenceSignal;
        RxChancmp = Rx_Aligned(i:NyquistRate:Ep);

        subSig = zeros(size(Rx_Aligned));
        subSig(i:NyquistRate:Ep) = RxChan;

        SigDiff = subSig - Rx_Aligned;
        Noise = RxChan-RxChancmp;
%         figure 
%         plot(20*log10(abs(RxChan))); hold on; plot(20*log10(abs(RxChancmp)),'r');plot(20*log10(abs(Noise)),':g');
%         legend('Reference Sig','Rx Aligned', 'Noise Sig')

        VarNoise = var(Noise)/2;
        NoiseVarArrayEstimate(num) = VarNoise;
        
        [ LLRData ] = LLR( LLRBitArray, RxSymbol, k, NoiseVar );
        [RxSoftDecision, RxHardDecision] = DecisionType( LLRData, k, TrCHFrameSize );
        [ RxSymbolIdx, RxTrCHFrameSize ] = SymbolIndex( RxHardDecision, k );
        BitErrorPacket(num) = mean(RxHardDecision(:) ~= TxPacket(:));
        SymbolErrorPacket(num) = mean(RxSymbolIdx(:) ~= TxSymbolIdx(:));
%         figure
%         plot(RxSymbol,'o') 
%         axis([-1 1 -1 1])
%         title(['Estimated Noise Variance: ', num2str(NoiseVar), ' Actual Noise Variance: ' , num2str(VarNoise)])
 
    end 
end 
figure
semilogy(EsNo_dB,SymbolErrorPacket,'b');

grid on

figure
loglog(NoiseVarArray,NoiseVarArrayEstimate)
hold on
loglog(NoiseVarArray,NoiseVarArray,'k')
legend('Noise Var vs Estimate Noise Var','Noise Var')