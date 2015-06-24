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
Nref = length(ReferenceSignal);
Ncp = length(CyclicPrefix);
ReferenceSignalPosition = 'midamble';
[ ReferenceSignalFreq, ReferenceSignalTime ] = CyclicCorr( ReferenceSignal,ReferenceSignal );

NyquistRate = 2;

[ FilterRRC ] = RootRaiseCosine( 0.22, 24, NyquistRate );
                    
M = [4];

NumM = length(M);         % Number of M-arys to be tested

% Configure Test
% --------------
EsNo_dB(:,1) = linspace(35, 0, 20); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations
EsNo_dBArray = zeros(NumEsNo,NumM);
EbNo_dBArray = zeros(NumEsNo,NumM); 

% Simulation configuration
% ------------------------
Parity = 'gCRC24A';
Decision = 'SoftDecision';

TTI = 4;

PacketSize = 2^8;
NumberOfPackets = 100;2^16;

% Coding 
% ------
Coding = 'Coded';


% Get Polynomial Generator Data, LLR Symbols, Packet Size
% -------------------------------------------------------
[ genPoly, NumberOfGenPoly, CodingRate, Memory, numStates, currentStates, statesBinary, numBranches, BranchLogic ] = GetGenPolyData;
TrCHSize = NumberOfGenPoly*(256);     

% Nodes Setup
% ------------
[ inputValue, numInput, PreviousState, Value, LogicValue, Decode ] = NodeSetup( Decision, numStates );
[ BranchNum, NextState ] = BranchLogicValues( numStates, inputValue, numInput, statesBinary, genPoly, NumberOfGenPoly, Memory, zeros(numStates,NumberOfGenPoly,numInput ), zeros(numStates,Memory,numInput ) );

 

%--------------------------------------------------------------------------
[ NewPacketSize, TrCHSize,  ParitySize ] = GetPacketSize( Coding, PacketSize, TrCHSize, Parity, Memory );

for m = 1:NumM
    % Get Symbol Array Data
    % ---------------------

    [ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M(m) );
    [ LLRBitArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M(m), k );
    EbNo_lin = EsNo_lin/k;
    EbNo_dB = 10*log10(EbNo_lin);       % EbNo Scaling
    EsNo_dBArray(:,m) = EsNo_dB;
    EbNo_dBArray(:,m) = EbNo_dB;
    % Theoretical Data
    % ---------------- 
    switch M(m)
        case 1
            bit_err_theo(:,1) = (1/k)*erfc(sqrt((k/Esnorm)*EbNo_lin));
            simb_error_theo(:,1) = erfc(sqrt((k/Esnorm)*EbNo_lin)); 
        case 2
%             bit_err_theo(:,1) =
%             (1/(2*k))*erfc(sqrt((k/Esnorm)*EbNo_lin));
            bit_err_theo = 0.5.*(1-sqrt((EbNo_lin)./(EbNo_lin+1)));
            simb_error_theo(:,1) = (1/2)*erfc(sqrt((k/Esnorm)*EbNo_lin));
        case 4
            bit_err_theo = 0.5.*(1-1*(1+1./EbNo_lin).^(-0.5)); 
%               bit_err_theo = 0.5.*(1-sqrt((EbNo_lin)./(EbNo_lin+1)));
%             bit_err_theo(:,1) =
%             ((1/k)*erfc(sqrt((k/Esnorm)*EbNo_lin))-(1/(4*k))*erfc(sqrt((k/Es)*EbNo_lin)).^2);
            simb_error_theo(:,1) = erfc(sqrt((k/Esnorm)*EbNo_lin))-(1/4)*erfc(sqrt((k/Es)*EbNo_lin)).^2;
        case 16
%             bit_err_theo(:,1) = (3/8)*erfc(sqrt((k/Esnorm)*EbNo_lin));%((3/2)*(3/16)*erfc(sqrt((k/Es)*EbNo_lin))-(9/16)*(21/64)*erfc(sqrt((k/Es)*EbNo_lin)).^2);
            simb_error_theo(:,1) = ((3/2)*erfc(sqrt((k/Esnorm)*EbNo_lin))-(9/16)*erfc(sqrt((k/Es)*EbNo_lin)).^2);
                                    bit_err_theo = 3/8 * ( 1 - sqrt(2/5*EbNo_lin./(1+2/5*EbNo_lin)) ) ...
                            + 1/4 * ( 1 - sqrt(18/5*EbNo_lin./(1+18/5*EbNo_lin)) ) ...
                            - 1/8 * ( 1 - sqrt(10*EbNo_lin./(1+10*EbNo_lin)) );
    end
    MBitErrorTheo(:,m) = bit_err_theo;
    MSymbolErrorTheo(:,m) = simb_error_theo;
    BitError = zeros(1,NumEsNo);
    SymbolError = zeros(1,NumEsNo);
    CRCER = zeros(1,NumEsNo);
     for i = 1:NumEsNo
        BitErrorPacket = zeros(1,NumberOfPackets);
        SymbolErrorPacket = zeros(1,NumberOfPackets);
        totalCRCErrors = zeros(1,NumberOfPackets);
        RayVariables = zeros(NumberOfPackets,1);
        NoiseVariance = zeros(NumberOfPackets,1);
        EstimatedNV = zeros(NumberOfPackets,1);
        PacketErrorStop = 0;
        PacketCount = 0;
        for ii = 1:NumberOfPackets
            close all
            TxPacket = double(rand(1,NewPacketSize) > 0.5);
            meanTx = mean(TxPacket); 
            varTx = var(TxPacket);
              
            [TxPacket, TxPacketLength] = GetCRC('Transmitter', TxPacket, Parity);
            xi = EncodePacket( TxPacket, genPoly, NumberOfGenPoly, Memory, Coding );
            [ yi, Xi, Yi ] = Interleaver1st( xi, TTI );
            [ yi, DeltaNij , Nij , Fi, TrCHArray ] = GetTrCHParameters( yi, TrCHSize );
            [ Eini , Eplus , Eminus, DeltaNi ] = GetRateMatchParametersNonTurbo ( DeltaNij , Nij , Fi );
            [ TxBits ] = RateMatching( yi, Yi, DeltaNij,  Eini , Eplus , Eminus, TrCHArray, Fi );
            [ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( TxBits, k );

            % Modulation, Transmit over Medium, Add Burst Error, Compute LLR
            % --------------------------------------------------------------
            [ TxSymbol ] = SymbolArray(TxSymbolIdx);    % Gets Symbols to be Transmitted for Symbol Array
            Nsymb = length(TxSymbol);
 
            switch ReferenceSignalPosition
                 case 'midamble'
                     TxSymbolRe = reshape(TxSymbol,[],2);
                     Nre = length(TxSymbolRe(:,1));
                     TxSymbol_Plus_Ref = [TxSymbolRe(:,1).',CyclicPrefix, ReferenceSignal,TxSymbolRe(:,2).'];
            end 
            numSymbols = length(TxSymbol_Plus_Ref);

            % Up Sampled Signal 
            % -----------------
            [ TxUpSample, NTx] = OverSample( TxSymbol_Plus_Ref, NyquistRate );
            [ TxRefOVS, NTxRef] = OverSample( ReferenceSignal, NyquistRate );

            % RRC Filter
            % ----------
            TxRRC = conv(FilterRRC*NyquistRate,TxUpSample);
            halfNsym = Nsymb/2;
            
            [rrcMax,rrcIndex] = max(FilterRRC);
            iref = (Nre*NyquistRate)+rrcIndex*2+1;
            iref = iref+ (Ncp*NyquistRate)-1;
            % Channel Delay
            % -------------
            SymbolyDelay = 15;
            Delay = 1;%randi(SymbolyDelay,1); 
            h = zeros(1,SymbolyDelay*NyquistRate);
            h(Delay*NyquistRate) = 1; % Delay Unit Power
            SymbolDelay = Delay-1; 
 
           % Convole transmitted signal will channel
            TxDelay = conv(h,TxRRC);
  

            % Reciever Noise, Rayleigh Fade
            % -------------------
            [ RxSymbol, AWGN ,No, NoiseVar, rayFade ] = RayleighFade( TxDelay, Es, EsNo_lin(i), NyquistRate, 'Rayleigh1' );
            varAWGN = var(AWGN)/2; % Approximately the Nyquist Rate * spectal Noise Density
            RayVariables(ii) = rayFade;
        
            % 2nd RRC Filter
            % --------------
            RxRRC2 = conv(FilterRRC,RxSymbol);
            Rxawgn = conv(FilterRRC,AWGN);
            %RxAWGN = TxUpSample + Rxawgn(1:length(TxUpSample));
            RxAWGN = RxRRC2+Rxawgn;
            varRxawgn = var(Rxawgn)/2; 
            NvPacket(ii) = var(rayFade+Rxawgn)/2;
            f =-NyquistRate/2:NyquistRate/(length(RxRRC2)-1):NyquistRate/2;
            f2 =-NyquistRate/2:NyquistRate/(length(TxDelay)-1):NyquistRate/2;
            f3 =-NyquistRate/2:NyquistRate/(length(AWGN)-1):NyquistRate/2;
            
%             figure
%             subplot(2,1,1)
%              plot(f2,20*log10(fftshift(abs(fft(TxDelay))))); hold on; plot(f3,20*log10(fftshift(abs(fft(AWGN)))),'r');
% %             plot(f,abs(RxRRC2));hold on; plot(f,abs(Rxawgn),'r')
%             title(['Noise Spectral Density: ', num2str(No)],'FontSize',14);
%             legend('Rx Signal', ['AWGN, Variance: ',num2str(varAWGN)],'FontSize',14) 
%             xlabel('Frequency', 'FontSize',14);
%             ylabel('Power dB', 'FontSize',14);
%             set(gca, 'FontSize',14); 
%             subplot(2,1,2)
%             plot(f,20*log10(fftshift(abs(fft(RxRRC2))))); hold on; plot(f,20*log10(fftshift(abs(fft(Rxawgn)))),'r');plot(f,20*log10(fftshift(abs(fft(RxAWGN)))),'k');
%             legend('Rx Signal',['AWGN, Variance: ',num2str(varRxawgn)],'Recieved Signal', 'FontSize',14)
%             title(['Noise Spectral Density: ', num2str(No)], 'FontSize',14);
%             xlabel('Frequency', 'FontSize',14);
%             ylabel('Power dB', 'FontSize',14);
%             set(gca, 'FontSize',14); 
             
            Nrx = length(RxAWGN);
            MidPoint = round(Nrx/2); 
%             iref = MidPoint -round((Nref*NyquistRate)/2);
            jref = iref+(Nref*NyquistRate)-1;
            
            [ TxRefFreq, TxRefTime ] = CyclicCorr( TxRefOVS, TxRefOVS );
            LocalReplica = abs(TxRefTime)/Nref;
            [ ChannelEstimateFreq, ChannelEstimateTime ] = CyclicCorr( RxAWGN(iref:jref),TxRefOVS );
            ChannelEstimateFreqdB = 20.*log10(abs(ChannelEstimateFreq));
            ChannelEstimateNorm = abs(ChannelEstimateTime)/Nref;
            [Peak, Index] = max(ChannelEstimateNorm);
            if Index~= 1
                iref = iref + Index-1;
            end
            jref = iref+(Nref*NyquistRate)-1;
            
%             figure  
%             subplot(2,1,1)
%             plot(abs(TxRefTime),'r'); hold on; plot(abs(ChannelEstimateTime));
%             title(['Symbol Delay: ', num2str(SymbolDelay)], 'FontSize',14)
%             xlabel('Index', 'FontSize',14);
%             ylabel('Magnitude', 'FontSize',14);
%             legend('Local Replica','Channel Estimate', 'FontSize',14)
%             set(gca, 'FontSize',14); 
  
        % Cyclic Correlation Frequency Domain
        if jref<=length(RxAWGN)

            [ SyncFreq, SyncTime ] = CyclicCorr( RxAWGN(iref:jref),TxRefOVS );
            SyncNorm = SyncTime/Nref;
            SyncNormABS = abs(SyncTime)/Nref;

            RefNorm = LocalReplica(1:NyquistRate:end);
            SyncNormDec = SyncNorm(1:NyquistRate:end);
%             subplot(2,1,2)
%             plot(LocalReplica,'r');  hold on; plot(SyncNormABS);
%             title(['Symbol Delay: ', num2str(SymbolDelay)], 'FontSize',14)
%             xlabel('Index', 'FontSize',14);
%             ylabel('Magnitude', 'FontSize',14);
%             legend('Norm Local Replica','Norm Channel Estimate', 'FontSize',14)
%             set(gca, 'FontSize',14); 
%                 figure  
%                 plot(abs(SyncNormDec)); hold on; plot(RefNorm,'r');
            AlignmentTerm = 1/SyncNormDec(1);

            % Tx Original and Tx phas and Gain Aligned
            % Start and End Position 
            Pstart = (iref-(NyquistRate*Ncp));
            Pend = jref;

            Refstart = iref;
            Refend = jref;
            RxSymStart = Pstart-NyquistRate;
            RxSymEnd = Pstart - NyquistRate*Nsymb/2;

%                 RxSymbol_Unaligned = [RxAWGN(RxSymEnd:NyquistRate:RxSymStart) RxAWGN(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];
%                 figure
%                 subplot(2,1,1) 
%                 plot(RxSymbol_Unaligned,'o') 
%                 legend('Rx Signal', 'FontSize',14)
%                 xlabel('I', 'FontSize',14);
%                 ylabel('Q', 'FontSize',14);
%                 title('Unaligned Constellation Points', 'FontSize',14)
%                 set(gca, 'FontSize',14);
            Rx_Aligned = RxAWGN/SyncNormDec(1);   

              % Recieved Reference Sig 
              RxRef = Rx_Aligned(Refstart:NyquistRate:Refend);
            if Pend+NyquistRate*Nsymb/2<=Nrx

                  RxSymbol = [Rx_Aligned(RxSymEnd:NyquistRate:RxSymStart) Rx_Aligned(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];

                  Noise = ReferenceSignal - RxRef;
                  VarNoise = var(Noise)/2;
                  NvEstimatePacket(ii) = VarNoise;
                  
%                   subplot(2,1,2)
%                   plot(RxSymbol,'o') 
%                   legend('Rx Signal Aligned', 'FontSize',14) 
%                   xlabel('I', 'FontSize',14);
%                   ylabel('Q', 'FontSize',14);
%                   title('Aligned Constellation Points', 'FontSize',14)
%                   set(gca, 'FontSize',14);
%                   axis([-1 1 -1 1])
            else
                RxSymbol = RxAWGN(1:Nsymb);
 
            end
        else
            RxSymbol = RxAWGN(1:Nsymb);
        end 
            
%             figure
%             plot(RxSymbol,'o');
             
            [ LLRData ] = LLR( LLRBitArray, RxSymbol, k, No );

            % Soft or Hard Decision
            % ---------------------
            [RxSoftDecision, RxHardDecision] = DecisionType( LLRData, k, TrCHFrameSize );
            Ni = Nij;  % Original Number of bits Per Frame
             switch Decision
                 case 'HardDecision'
                     RxMessage = RxHardDecision;
                 case 'SoftDecision'
                     [ RxMessage ] = DeRateMatching( RxSoftDecision, Yi, DeltaNij,  Eini , Eplus , Eminus, zeros(Ni,Fi), Fi );
             end

            % DeInterLeaver
            % -------------
            [ RxMessage ] = DeInterleaver1st( RxMessage, Fi, Ni , TTI );

            if strcmp(Coding,'Coded')
                [ BlockCode ] = BlockCodeWord( RxMessage, NumberOfGenPoly,  TxPacketLength );

                [ BMetric, BlockCode ] = BranchMetric( Decision, BlockCode, BranchLogic, numBranches, NumberOfGenPoly,TxPacketLength );

                % Initialise Node
                % --------------
                Stages  = TxPacketLength+1;
                TrellisValue = repmat(Value,1,Stages);
                TrellisPreviousState = repmat(PreviousState,1,Stages);
                TrellisDecode = repmat(Decode,1,Stages);
                TrellisValue(1,1) = 0;
 
                % Taverse Trellis Map
                % -------------------
                [ TrellisValue, TrellisPreviousState ,TrellisDecode ] = TrellisMap( Decision, BMetric , BranchNum, BlockCode, NextState, TrellisValue, TrellisPreviousState, TrellisDecode, currentStates, Stages, numStates ,numInput );
                [ RxMessage] = DecodePacket( Decision, TrellisValue, TrellisPreviousState ,TrellisDecode, Stages,BranchNum , TxPacket );
            end
            % Check For Errors
            % ----------------
            [ RxSymbolIdx, RxTrCHFrameSize ] = SymbolIndex( RxHardDecision, k );
            BitErrorPacket(ii) = mean(RxHardDecision(:) ~= TxBits(:));
            NewBitErrorPacket(ii) = mean(RxMessage(:) ~= TxPacket(:));
            SymbolErrorPacket(ii) = mean(RxSymbolIdx(:) ~= TxSymbolIdx(:));
            if ~isnan(Parity)
                [ DecodedMessage, crcError ] = GetCRC('Reciever', RxMessage, Parity);
                totalCRCErrors(ii) = crcError;
            elseif BitErrorPacket(ii)>0
                totalCRCErrors(ii) = 1;
            end
            % Stop current test if enough packet errors have occured
            % ------------------------------------------------------
            PacketErrorStop = PacketErrorStop + totalCRCErrors(ii);
            PacketCount = PacketCount+1;
%             if PacketErrorStop >50
%                 break
%             end
        end
        CRCERTotal(1,i)= sum(totalCRCErrors);
        CRCER(1,i) = sum(totalCRCErrors)/PacketCount;
        BitError(1,i) = sum(BitErrorPacket)/PacketCount;
        SymbolError(1,i) = sum(SymbolErrorPacket)/PacketCount;
        NoiseVarArray(i) = sum(NvPacket)/PacketCount;
        NewBitError(1,i) = sum(NewBitErrorPacket)/PacketCount;
%         NoiseVarArrayEstimate(i)=sum(NvEstimatePacket)/PacketCount;
        rayTemp(i)  = mean(abs(RayVariables).^2);
     end
        % Determine Theoretical Block Error Rates
        % ---------------------------------------
        N = TxPacketLength;
        BLER_theo(:,1)  = (1-(1-bit_err_theo).^N);
        MBLERTheo(:,m)= BLER_theo(:,1);
        MBitError(:,m) = BitError(1,:);
        MBLER(:,m)= CRCER(1,:);
end
        
figure
semilogy(EbNo_dBArray,MBitErrorTheo);

hold on
semilogy(EbNo_dBArray,MBitError,':o');
bit_err_theo(:,1) =((1/k)*erfc(sqrt((k/Esnorm)*EbNo_lin))-(1/(4*k))*erfc(sqrt((k/Es)*EbNo_lin)).^2);
semilogy(EbNo_dB,bit_err_theo,'r');
semilogy(EbNo_dBArray,MBLER,':*'); 

semilogy(EbNo_dB,NewBitError,'k-o');

grid on
xlabel('Eb/No, dB', 'FontSize',14);
ylabel('Bit Error Rate', 'FontSize',14);
legend('Thereory Bit Error BPSK','Thereory Bit Error QPSK','Thereory Bit Error 16-QAM', 'Bit Error BPSK','Bit Error QPSK','Bit Error 16-QAM', 'BLER BPSK','BLER QPSK','BLER 16-QAM', 'FontSize',14) 
% legend('Thereory Bit Error QPSK','Bit Error QPSK','BLER QPSK','FontSize',14) 
% legend('Thereory Bit Error BPSK','Thereory Bit Error QPSK','Thereory Bit
% Error 16-QAM', 'Bit Error BPSK','Bit Error QPSK','Bit Error 16-QAM','Thereory BLER BPSK','Thereory BLER QPSK','Thereory BLER 16-QAM', 'BLER BPSK','BLER QPSK','BLER 16-QAM', 'FontSize',14) 
title('Bit Error and BLER Using Channel Estimate', 'FontSize',14)
set(gca, 'FontSize',14); 
% % axis([-6.020599913279623 35 10^-4 10^-0 ])

% figure
% subplot(2,1,1)
% semilogy(EbNo_dBArray,NoiseVarArrayEstimate)
% hold on
% semilogy(EbNo_dBArray,NoiseVarArray,'k')
% legend('Noise Var Estimate','Noise Var')
% title('Noise Variance Estimate for LLR', 'FontSize',14)
% axis([min(EbNo_dBArray) max(EbNo_dB) min(NoiseVarArrayEstimate) max(NoiseVarArrayEstimate) ])
% % axis([min(NoiseVarArray) max(NoiseVarArray) min(NoiseVarArrayEstimate) max(NoiseVarArrayEstimate) ])
% xlabel('EbNo dB', 'FontSize',14);
% ylabel('Noise Variance', 'FontSize',14);
% set(gca, 'FontSize',14); 
% grid on
% % flipdim(h,1));
% 
% Error = abs(NoiseVarArray -NoiseVarArrayEstimate)*100;
% subplot(2,1,2)
% plot(EbNo_dBArray,Error)
% title('Noise Variance Error', 'FontSize',14)
% xlabel('Error Percentage', 'FontSize',14);
% ylabel('EbNo dB', 'FontSize',14);
% grid on
% 
% figure
% plot(0:length(RayVariables)-1,20*log10(abs(RayVariables)),'b');
% title('Rayleigh Fade')
% grid on

