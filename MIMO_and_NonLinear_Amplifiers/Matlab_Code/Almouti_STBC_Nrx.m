close all
clear all
% Number of Transmit Antenna's
NTx = 2;

NRxArray = [1 2 3 4 5 6];
% Configure Test
% --------------
EsNo_dB(:,1) = linspace(35, 0, 10); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations

% Constellation
M = 4;
N = 100;
NumberOfPackets = 10000;
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
[ LLRBitArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M, k );
EbNo_lin = EsNo_lin/k;
EbNo_dB = 10*log10(EbNo_lin);       % EbNo Scaling
NumRxAntennas = length(NRxArray);
BitError = zeros(NumEsNo,1);
NrxBitError = zeros(NumEsNo,NumRxAntennas);
Leg{1} = ['Therory Ntx = 1, Nrx = 1'];
Leg{2}=['Therory Ntx = 1, Nrx = 2 Dual Diversity'];
Leg{3}= ['Theory Ntx = 2, Nrx = 1 Alamouti'];
count = 4;
for antNum = 1:NumRxAntennas
    % Number of Recieve Antenna's
    NRx = NRxArray(antNum);
    for i = 1:NumEsNo
        for ii = 1:NumberOfPackets
            % Generate equally likely binary values
            txBits = rand(1,N*k) > 0.5;
            % Symbol idx for symbol array
            [ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( txBits, k );
            % Symbols
            [ sn ] = SymbolArray(TxSymbolIdx)*(1/sqrt(NTx)); 
            s=reshape(sn,2,length(sn)/2);
            
            sTemp = zeros(NTx, N);

            c1 = s;                         % [s1;s2]

            c2conj = conj(s([2, 1],:));               % [-s2*;s1*;]
            c2conj(1,:) = -c2conj(1,:);                     

            sTemp(:,1:2:end) = c1;
            sTemp(:,2:2:end) = c2conj;

        

            % Channel 
            % ------- 
            x = randn(NRx, NTx); y = randn(NRx, NTx); % Random Variables X and Y with Gaussian Distribution mean zero variance one
            h =(x + 1i.*y)/sqrt(2);                   % Normalised Ralyiegh random Variable
            x1=[1 0 ];
            x2=[0 1 ];

            X1 = [x1; x2];
            X2 = [-x2; x1];
            H1 = [h*X1];
            H2 = conj([h*X2]);
            H = [H1;H2];    
%             x1=[1 0];
%             x2=[0 1];
%             xR = [x1 ;x2];
%             xConj = [-x2; x1];
%             H1 = h*xR; 
%             H2 = conj(h)*xConj;
%             H3 = [H1;H2];


            % Effects of Channel
            HsTemp = h*sTemp;
            Hs = H*s;
            hs = reshape(Hs,1,[]);


            % AWGN
            No = (Es./EsNo_lin(i));    % Noise Spectral Density
            NoiseVar = No;             % Noise Variance
            n = sqrt(No)*(randn(NRx*NTx,length(Hs))+1i*randn(NRx*NTx,length(Hs)))*(1/sqrt(2));     % Normalised AWGN in respect to Esnorm

            % Receiver
            y = Hs+n; 

            HMod = pinv(H);
            % Singular value decomposition
            [U S V] = svd(H);
            HMod2 = V*inv(S'*S)*S'*U';
            sHat = HMod*y; 
            sHat = reshape(sHat,1,N)*sqrt(NTx);
            HPower = sum(H.'.*conj(H.'),1);
            sHat1 = (H'*H*s+H'*n);
            sHat12 = (s+HMod*n);
            sHat12 = reshape(sHat12,1,N)*sqrt(NTx);
            sHat1 = reshape(sHat1,1,N)./HPower(1);
    %         Htemp = [conj(h(1)) h(2); conj(h(2)) -h(1)];
    %         sHat2 = Htemp*y;
    %         sHat2 = reshape(sHat2,1,N)./HPower(1);

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
    NrxBitError(:,antNum)= BitError;
    Leg{count} = ['Sim Ntx = 2, Nrx = ',num2str(NRx) ,' Alamouti'];
    count = count+1;
end
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbNo_lin).^(-0.5)); 

p = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 
pAlamouti = 1/2 - 1/2*(1+2./EbNo_lin).^(-1/2);
theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti)); 

colLit=hsv(antNum+3);
figure
semilogy(EbNo_dB,theoryBer_nRx1,'col',colLit(1,:))
hold on
semilogy(EbNo_dB,theoryBerMRC_nRx2,'col',colLit(2,:))
semilogy(EbNo_dB,theoryBerAlamouti_nTx2_nRx1,'col',colLit(3,:))
for n = 0:antNum-1
    semilogy(EbNo_dB,NrxBitError(:,n+1),'-p','col',colLit(3+n,:))
end
grid on
xlabel('EbNo dB')
ylabel('BER')
title(['Spatial Rate = 1 Rank = ' num2str(rank(H)), ' Condition = ', num2str(cond(H))])
legend(Leg)