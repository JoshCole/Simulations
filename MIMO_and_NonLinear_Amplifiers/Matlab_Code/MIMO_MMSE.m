close all
clear all

% Number of Transmit Antenna's
NTx = 2;
% Number of Recieve Antenna's
NRx = 2;
NTxArray = [2];
% % Bandwith
% Bandwidth = 5e6;
% ActualBandwidth = 4.5e6;
% % Subcarrier Frequency Spacing
% deltaF = 15e3;
% % Symbol Duration
% T = 1/deltaF;
% % Sample Size
% Nfft = 512;
% Ts = T/Nfft;
% % Number of Subcarriers
% Nsubcarriers = ActualBandwidth/deltaF;
% SubCarrierIndex = [-Nsubcarriers/2:-1 1: Nsubcarriers/2];
% NumberOfSlots = 20;
% NumberOfSymbolsInRB = 7;
% % Number of OFDM Symbols
% NumberOfOFDMSymbols = NumberOfSymbolsInRB*NumberOfSlots;
% 
% % Modulation Scheme
% M = [4];
% NumM = length(M);

% Configure Test
% --------------
EsNo_dB(:,1) = linspace(35, 0, 10); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations

% OFDMsymbol = zeros(1,Nfft);
% Ncp = 36;
% ReferenceSymbol = zeros(1,Nfft);
% Ref = ones(1,Nsubcarriers); 
% ReferenceSymbol(SubCarrierIndex+Nfft/2+1) = Ref;
% 
% powerRef = mean(abs(Ref).^2);
% Lref = ReferenceSymbol;
% ReferenceSymbol = ifft(ReferenceSymbol,Nfft)*sqrt(Nfft);
% 
% RefPrefix = ReferenceSymbol((Nfft-Ncp)+1:Nfft);
% LocalReplica = ReferenceSymbol;
% ReferenceSymbol = [ RefPrefix ReferenceSymbol];
% % ReferenceSymbol = ReferenceSymbol*sqrt(OverSampleRate)*sqrt(Nfft);
% powerReferenceSymbol = mean(abs(ReferenceSymbol).^2);
M = 4;
PacketSize = 256;
NumberOfPackets = 1000;
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
[ LLRBitArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M, k );
EbNo_lin = EsNo_lin/k;
EbNo_dB = 10*log10(EbNo_lin);       % EbNo Scaling
AntennaVariation = length(NTxArray);
figure
colLit=hsv(AntennaVariation);
for Ant = 1:AntennaVariation
    NTx = NTxArray(Ant);
    NRx = NTx;
    BitError = zeros(1,NumEsNo);
    BitErrorPacket = zeros(1,NumberOfPackets);
     for i = 1:NumEsNo
        for ii = 1:NumberOfPackets
            % Generate equally likely binary values
            TxPacket = double(rand(NTx,PacketSize*k) > 0.5);
            % Symbol idx for symbol array
            [ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( TxPacket, k );
            % Symbols
            [ sn ] = SymbolArray(TxSymbolIdx); 
            Nsym = length(sn);
            % Transmit Antennas
    %         sn = reshape(Symbols,NTx,Nsym/NTx);
    %         sn = repmat(Symbols,[],NTx);

            % Channel 
            % ------- 
            x = randn(NRx, NRx); y = randn(NRx, NRx); % Random Variables X and Y with Gaussian Distribution mean zero variance one
            H =(x + 1i.*y)/sqrt(2);       % Normalised Ralyiegh random Variable
            Hs = H*sn;
            
           % AWGN
            N = Nsym;
            No = (Es./EsNo_lin(i));    % Noise Spectral Density
            NoiseVar = No;             % Noise Variance
            n = sqrt(No)*(randn(NRx,N)+1i*randn(NRx,N))*(1/sqrt(2));     % Normalised AWGN in respect to Esnorm

            % MMSE
            W_H =  inv(H'*H + No*eye(NTx))*H';

            y = Hs+n;
            sHat = W_H*y;
    %         if ii ==1
    %             [U,S,V] = svd(H)
    %             conditionNumber = cond(H)
    %             colList = hsv(NRx);
    %             figure
    %             for cl = 1:NRx
    %             plot(y(cl,:),'p', 'col',colList(cl,:))
    %             hold on
    %             end
    %             
    %             figure
    %             for cl = 1:NRx
    %             plot(yHat(cl,:),'p', 'col',colList(cl,:))
    %             hold on
    %             end
    %             axis([-1 1 -1 1]);
    %             pause
    %         end

            Rxbits = zeros(size(TxPacket)); 
            for ntx = 1:NTx
                 LLRData  = LLR( LLRBitArray, sHat(ntx,:), k, NoiseVar );
                % Soft or Hard Decision
                % ---------------------
                [RxSoftDecision, RxHardDecision] = DecisionType( LLRData, k, PacketSize );
                Rxbits(ntx,:) = RxHardDecision;
            end
            [ RxSymbolIdx, TrCHFrameSize ] = SymbolIndex( Rxbits, k );
            [ rxsn ] = SymbolArray(RxSymbolIdx); 
            BitErrorPacket(ii) = mean(TxPacket(:) ~= Rxbits(:));
            SymbolErrorPacket(ii) = mean(sn(:) ~= rxsn(:));
        end
        BitError(i) = sum(BitErrorPacket)/NumberOfPackets;
        SymbolError(i) = sum(SymbolErrorPacket)/NumberOfPackets;
     end

    
    %  semilogy(EbNo_dB,bit_err_theo)
    %  hold on

     theoryBer = 0.5.*(1-1*(1+1./EbNo_lin).^(-0.5)); 

%      semilogy(EbNo_dB,theoryBer,'r')
%      hold on
     semilogy(EbNo_dB,BitError,'col',colLit(Ant,:))
     hold on
     leg{Ant} =['MMSE Simulation Nrx = Ntx = ',num2str(NTx)];

%       figure
%     %  semilogy(EbNo_dB,bit_err_theo)
%     %  hold on
%      semilogy(EsNo_dB,SymbolError)
%      legend('MMSE Simulation')
%      grid on
end
     legend(leg)
     grid on
     xlabel('EbNo dB')
     ylabel('BER')
