close all
clear all
% Number of Transmit Antenna's
NTx = 4;
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
NumberOfPackets = 10000;
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
        s=reshape(sn,NTx,length(sn)/NTx);
        
        sTemp = zeros(NTx, 2*N);
                
        c1 = s;                         % [s1;s2;s3;s4]

        c2 = s([2, 1, 4, 3],:);         % [-s2;s1;-s4;s3]
        c2(1,:) = -c2(1,:);             
        c2(3,:) = -c2(3,:);

        c3 = s([3, 4, 1, 2],:);         % [-s3;s4;s1;-s2]
        c3(1,:) = -c3(1,:);
        c3(4,:) = -c3(4,:);

        c4 = s([4, 3, 2, 1],:);         % [-s4;-s3;s2;s1]
        c4(1,:) = -c4(1,:);
        c4(2,:) = -c4(2,:);

        c1conj = conj(c1);              % [s1*;s2*;s3*;s4*]
        c2conj = conj(c2);              % [-s2*;s1*;-s4*;s3*]
        c3conj = conj(c3);              % [-s3*;s4*;s1*;-s2*]
        c4conj = conj(c4);              % [-s4*;-s3*;s2*;s1*]
        
        sTemp(:,1:8:end) = c1;
        sTemp(:,2:8:end) = c2;
        sTemp(:,3:8:end) = c3;
        sTemp(:,4:8:end) = c4;
        sTemp(:,5:8:end) = c1conj;
        sTemp(:,6:8:end) = c2conj;
        sTemp(:,7:8:end) = c3conj;
        sTemp(:,8:8:end) = c4conj;

        
        % Channel 
        % ------- 
        x = randn(NRx, NTx); y = randn(NRx, NTx); % Random Variables X and Y with Gaussian Distribution mean zero variance one
        h =(x + 1i.*y)/sqrt(2);                   % Normalised Ralyiegh random Variable
%         X1 = [x1; x2; x3; x4;];
%         X2 = [x2; -x1; x4; -x3;];
%         X3 = [x3; -x4; -x1; x2;];
%         X4 = [x4; x3; -x2; -x1];


        x1=[1 0 0 0];
        x2=[0 1 0 0];
        x3=[0 0 1 0];
        x4=[0 0 0 1];
        X1 = [x1; x2; x3; x4;];
        X2 = [-x2; x1; -x4; x3;];
        X3 = [-x3; x4; x1; -x2;];
        X4 = [-x4; -x3; x2; x1];
        H1 = [h*X1;h*X2; h*X3;h*X4];
        H2 = conj(H1); 

        H = [H1; H2];    

        % Effects of Channel
        HsTemp = h*sTemp;
        Hs = H*s;
        hs = reshape(Hs,1,[]);
        
        % AWGN
        No = (Es./EsNo_lin(i));    % Noise Spectral Density
        NoiseVar = No;             % Noise Variance
        n = sqrt(No)*(randn(2*NRx*NTx,length(Hs))+1i*randn(2*NRx*NTx,length(Hs)))*(1/sqrt(2));     % Normalised AWGN in respect to Esnorm
        
        % Receiver
        y = Hs+n; 
        
        HMod = pinv(H);
        sHat = HMod*y; 
        sHat = reshape(sHat,1,N)*sqrt(NTx);
        sHat12 = (s+HMod*n);
        sHat12 = reshape(sHat12,1,N)*sqrt(NTx);
        
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
semilogy(EsNo_dB,theoryBerAlamouti_nTx2_nRx1,'col',colLit(3,:))
semilogy(EsNo_dB,BitError,'-p','col',colLit(3,:))
grid on
xlabel('EbNo dB')
ylabel('BER')
title(['Spatial Rate = 1/2 Rank = ' num2str(rank(H)), ' Condition = ', num2str(cond(H))])
legend('Therory Ntx = 1, Nrx = 1', 'Therory Ntx = 1, Nrx = 2 MRC', 'Theory Ntx = 2, Nrx = 1 Alamouti', ['Sim Ntx = 4, Nrx = ',num2str(NRx) ,' Alamouti'])