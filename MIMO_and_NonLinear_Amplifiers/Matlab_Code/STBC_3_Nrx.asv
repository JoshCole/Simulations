close all
clear all
% Number of Transmit Antenna's
NTx = 3;
% Number of Recieve Antenna's
NRx = 1;

% Configure Test
% --------------
EsNo_dB(:,1) = linspace(35, 0, 10); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations

% Constellation
M = 4;
N = 400;
NumberOfPackets = 1000;
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
[ LLRBitArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M, k );
EbNo_lin = EsNo_lin/k;
EbNo_dB = 10*log10(EbNo_lin);       % EbNo Scaling

for i = 1:NumEsNo
    for ii = 1:NumberOfPackets
        % Generate equally likely binary values
        txBits = rand(1,N*k) > 0.5;
        % Symbol idx for symbol array
        [ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( txBits, k );
        % Symbols
        [ sn ] = SymbolArray(TxSymbolIdx)*(1/sqrt(NTx)); 
        s=reshape(sn,4,length(sn)/4);
        
        % Channel 
        % ------- 
        x = randn(NRx, NTx); y = randn(NRx, NTx); % Random Variables X and Y with Gaussian Distribution mean zero variance one
        h =(x + 1i.*y)/sqrt(2);                   % Normalised Ralyiegh random Variable
        x1=[1 0 0 0];
        x2=[0 1 0 0];
        x3=[0 0 1 0];
        x4=[0 0 0 1];
        X1 = [x1; x2; x3;];
        X2 = [-x2; x1; -x4;];
        X3 = [-x3; x4; x1;];
        X4 = [-x4; -x3; x2;];
        H1 = [h*X1; h*X2; h*X3; h*X4];
        H2 = conj([h*X1; h*X2; h*X3;h*X4]);  
        H = [H1;H2];    
        % Effects of Channel
        Hs = H*s;
        hs = reshape(Hs,1,[]);
        
        % AWGN
        No = (Es./EsNo_lin(i));    % Noise Spectral Density
        NoiseVar = No;             % Noise Variance
        n = sqrt(No)*(randn(8,length(Hs))+1i*randn(8,length(Hs)))*(1/sqrt(2));     % Normalised AWGN in respect to Esnorm
        
        % Receiver
        y = Hs+n; 
        
        HMod = pinv(H);
        sHat = HMod*y; 
        sHat = reshape(sHat,1,N);
        HPower = sum(H.'.*conj(H.'),1);
        sHat1 = (H'*H*s+H'*n);
        sHat1 = reshape(sHat1,1,N)./HPower(1);
        
        LLRData  = LLR( LLRBitArray, sHat, k, NoiseVar );
        [RxSoftDecision, RxHardDecision] = DecisionType( LLRData, k, N );
        [ RxSymbolIdx, TrCHFrameSize ] = SymbolIndex( RxHardDecision, k );
        [ rxsn ] = SymbolArray(RxSymbolIdx); 

        BitErrorPacket(ii) = mean(txBits(:) ~= RxHardDecision(:));
        SymbolErrorPacket(ii) = mean(sn(:) ~= rxsn(:));
    end
        BitError(i) = sum(BitErrorPacket)/NumberOfPackets;
        SymbolError(i) = sum(SymbolErrorPacket)/NumberOfPackets;
end
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbNo_lin).^(-0.5)); 

p = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 
pAlamouti = 1/2 - 1/2*(1+2./EbNo_lin).^(-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti)); 

colLit=hsv(3);
figure
semilogy(EsNo_dB,theoryBer_nRx1,'col',colLit(1,:))
hold on
semilogy(EsNo_dB,theoryBerMRC_nRx2,'col',colLit(2,:))
semilogy(EsNo_dB,BitError,'-p','col',colLit(3,:))
grid on
xlabel('EbNo dB')
ylabel('BER')
title(['Spatial Rate = 1/2 Rank = ' num2str(rank(H)), ' Condition = ', num2str(cond(H))])
legend('Therory Ntx = 1, Nrx = 1', 'Therory Ntx = 1, Nrx = 2 MRC', ['Sim Ntx = 4, Nrx = ',num2str(NRx) ,' STBC'])