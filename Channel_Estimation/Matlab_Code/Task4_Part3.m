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

NyquistRate = 4;

[ FilterRRC ] = RootRaiseCosine( 0.22, 12, NyquistRate );
                    
M = [4];

NumM = length(M);         % Number of M-arys to be tested

% Configure Test
% --------------
EsNo_dB(:,1) = 10;%linspace(0, 35, 20); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations
EsNo_dBArray = zeros(NumEsNo,NumM);
EbNo_dBArray = zeros(NumEsNo,NumM);

% Simulation configuration
% ------------------------
Parity = 'gCRC24A';
Decision = 'HardDecision';

TTI = 4;

PacketSize = 2^8;
NumberOfPackets= 100;%2^20./PacketSize;

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
            bit_err_theo(:,1) = (1/k)*erfc(sqrt((k/Es)*EbNo_lin));
            simb_error_theo(:,1) = erfc(sqrt((k/Es)*EbNo_lin)); 
        case 2
            bit_err_theo(:,1) = (1/(2*k))*erfc(sqrt((k/Es)*EbNo_lin));
            simb_error_theo(:,1) = (1/2)*erfc(sqrt((k/Es)*EbNo_lin));
        case 4
            p  = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
            bit_err_theo = p.^2.*(1+2*(1-p)); 
%             L =2;
%             for k = 1:L
%                 a = (0.5*(1-sqrt(EbNo_lin./(1+EbNo_lin)))).^L;
%                 b = (0.5*(1-sqrt(EbNo_lin./(1+EbNo_lin)))).^k;
%             end
%             p  = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
%             bit_err_theo = p.^2.*(1+2*(1-p)); 
            bit_err_theo(:,1) = ((1/k)*erfc(sqrt((k/Es)*EbNo_lin))-(1/(4*k))*erfc(sqrt((k/Es)*EbNo_lin)).^2);
            simb_error_theo(:,1) = erfc(sqrt((k/Es)*EbNo_lin))-(1/4)*erfc(sqrt((k/Es)*EbNo_lin)).^2;
        case 16
            bit_err_theo(:,1) = (3/8)*erfc(sqrt((k/Es)*EbNo_lin));%((3/2)*(3/16)*erfc(sqrt((k/Es)*EbNo_lin))-(9/16)*(21/64)*erfc(sqrt((k/Es)*EbNo_lin)).^2);
            simb_error_theo(:,1) = ((3/2)*erfc(sqrt((k/Es)*EbNo_lin))-(9/16)*erfc(sqrt((k/Es)*EbNo_lin)).^2);
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
            powerRRC = mean(abs(TxRRC).^2);
            
            halfNsym = Nsymb/2;
            
            [rrcMax,rrcIndex] = max(FilterRRC);
%             iref = (Nref*NyquistRate)+rrcIndex*2+1;
%             iref = iref+ (Ncp*NyquistRate)-1;
%             iref = iref+6;
            iref = (Nref*NyquistRate*3)+rrcIndex*2+1;
            iref = iref+ (Ncp*NyquistRate)-1;
            iref = iref+14;
            % Channel Delay
            % -------------
            SymbolyDelay = 15;
            Delay = randi(SymbolyDelay,1); 
            h = zeros(1,SymbolyDelay*NyquistRate);
            h(Delay*NyquistRate) = 1; % Delay Unit Power
            SymbolDelay = Delay-1; 
            
            % Convole transmitted signal will channel
            TxDelay = conv(h,TxRRC);
            powerDelay = mean(abs(TxDelay).^2);
            
            % Reciever Noise, Rayleigh Fade
            % -------------------
            [ RxSymbol, AWGN ,No, NoiseVar, rayFade ] = RayleighFade( TxDelay, Es, EsNo_lin(i), NyquistRate, 'Rayleigh' );
            % 2nd RRC Filter
            % --------------
            RRC2 = conv(FilterRRC,RxSymbol);
            powerRRC2 = mean(abs(RRC2).^2);
            Rxawgn = conv(FilterRRC,AWGN);
            RxAWGN = RRC2+Rxawgn;
            powerRxAWGN = mean(abs(RxAWGN).^2);
            varRxawgn = var(Rxawgn); 
            
%             MidPoint = round(Nrx/2); 
            iref1 = iref;
%             iref = MidPoint -round((Nref*NyquistRate)/2);
            Nrx = length(RxAWGN);
            jref1 = iref1+(Nref*NyquistRate)-1;
            % Linear Correlation Channel Estimate
%             y = conv(RxAWGN, conj(flipdim(TxRefOVS,2)));
% 
%             yNorm = abs(y); 
%             t= linspace(1,824,length(y));
%             [Peak, Index] = max(yNorm);
%             figure
%             plot(t,abs(yNorm))  
             
            [ TxRefFreq, TxRefTime ] = CyclicCorr( TxRefOVS, TxRefOVS );
            LocalReplica = abs(TxRefTime)/Nref;
            [ ChannelEstimateFreq, ChannelEstimateTime ] = CyclicCorr( RxAWGN(iref1:jref1),TxRefOVS );
            ChannelEstimateFreqdB = 20.*log10(abs(ChannelEstimateFreq));
            ChannelEstimateNorm = abs(ChannelEstimateTime)/Nref;
            [Peak, Index] = max(ChannelEstimateNorm);
            if Index~= 1
                iref1 = iref1 + Index-1;
            end
            jref1 = iref1+(Nref*NyquistRate)-1;
%             figure  
%             subplot(2,1,1)
%             plot(abs(TxRefTime),'r'); hold on; plot(abs(ChannelEstimateTime));
%             title(['Symbol Delay: ', num2str(SymbolDelay)])
%             xlabel('Index');
%             ylabel('Magnitude');
%             legend('Local Replica','Channel Estimate')

        % Cyclic Correlation Frequency Domain
        if jref1<=length(RxAWGN)
            [ SyncFreq, SyncTime ] = CyclicCorr( RxAWGN(iref1:jref1),TxRefOVS );
            SyncNorm = SyncTime/Nref;
            SyncNormABS = abs(SyncTime)/Nref;
             [Peak2, Index2] = max(SyncNormABS);
            RefNorm = LocalReplica(1:NyquistRate:end);
            SyncNormDec = SyncNorm(1:NyquistRate:end);
%             subplot(2,1,2)
%             plot(LocalReplica,'r');  hold on; plot(SyncNormABS);
%             title(['Symbol Delay: ', num2str(SymbolDelay)])
%             xlabel('Index');
%             ylabel('Magnitude');
%             legend('Norm Local Replica','Norm Channel Estimate')

        
            % Tx Original and Tx phas and Gain Aligned
            % Start and End Position 
            Pstart = (iref1-(NyquistRate*Ncp));
            Pend = jref1;

            Refstart = iref1;
            Refend = jref1;
            RxSymStart = Pstart-NyquistRate;
            RxSymEnd = Pstart - NyquistRate*Nsymb/2;

%                 RxSymbol_Unaligned = [RxAWGN(RxSymEnd:NyquistRate:RxSymStart) RxAWGN(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];
%                 figure
%                 subplot(2,1,1) 
%                 plot(RxSymbol_Unaligned,'o') 
%                 legend('Rx Signal')
%             Rx_Aligned = RxAWGN*AlignmentTerm;   
            Rx_Aligned = RxAWGN;%(conj(SyncNormDec(1))*RxAWGN)/(SyncNormDec(1)*conj(SyncNormDec(1)));
            h1 =SyncNormDec(1);
              % Recieved Reference Sig 
              RxRef = Rx_Aligned(Refstart:NyquistRate:Refend);
            if Pend+NyquistRate*Nsymb/2<=Nrx

                  RxSymbol1 = [Rx_Aligned(RxSymEnd:NyquistRate:RxSymStart) Rx_Aligned(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];

                  Noise = ReferenceSignal - RxRef;
                  VarNoise = var(Noise);
                   NoiseArray(1,:) = RxRef;
                  NvPacketEstimate(ii) = VarNoise;
                   
%                   subplot(2,1,2)
%                   plot(RxSymbol,'o') 
%                   legend('Rx Signal Aligned') 
%                   axis([-1 1 -1 1])
            else
                RxSymbol1 = RxAWGN(1:Nsymb);
 
            end
        else
            RxSymbol1 = RxAWGN(1:Nsymb);
        end 
            Antenna1 = RxSymbol1;
            % Reciever Noise, Rayleigh Fade
            % -------------------
            [ RxSymbol, AWGN ,No, NoiseVar, rayFade ] = RayleighFade( TxDelay, Es, EsNo_lin(i), NyquistRate, 'Rayleigh1' );
             % 2nd RRC Filter
            % --------------
            RRC2 = conv(FilterRRC,RxSymbol);
            powerRRC2 = mean(abs(RRC2).^2);
            Rxawgn = conv(FilterRRC,AWGN);
            RxAWGN = RRC2+Rxawgn;
            powerRxAWGN = mean(abs(RxAWGN).^2);
            varRxawgn = var(Rxawgn); 
            
%             MidPoint = round(Nrx/2); 
            iref2 = iref;
%             iref = MidPoint -round((Nref*NyquistRate)/2);
            Nrx = length(RxAWGN);
            jref2 = iref2+(Nref*NyquistRate)-1;
            % Linear Correlation Channel Estimate
%             y = conv(RxAWGN, conj(flipdim(TxRefOVS,2)));
% 
%             yNorm = abs(y); 
%             t= linspace(1,824,length(y));
%             [Peak, Index] = max(yNorm);
%             figure
%             plot(t,abs(yNorm)) 
             
            [ TxRefFreq, TxRefTime ] = CyclicCorr( TxRefOVS, TxRefOVS );
            LocalReplica = abs(TxRefTime)/Nref;
            [ ChannelEstimateFreq, ChannelEstimateTime ] = CyclicCorr( RxAWGN(iref2:jref2),TxRefOVS );
            ChannelEstimateFreqdB = 20.*log10(abs(ChannelEstimateFreq));
            ChannelEstimateNorm = abs(ChannelEstimateTime)/Nref;
            [Peak, Index] = max(ChannelEstimateNorm);
            if Index~= 1
                iref2 = iref2 + Index-1;
            end
            jref2 = iref2+(Nref*NyquistRate)-1;
%             figure  
%             subplot(2,1,1)
%             plot(abs(TxRefTime),'r'); hold on; plot(abs(ChannelEstimateTime));
%             title(['Symbol Delay: ', num2str(SymbolDelay)])
%             xlabel('Index');
%             ylabel('Magnitude');
%             legend('Local Replica','Channel Estimate')

        % Cyclic Correlation Frequency Domain

        if jref2<=length(RxAWGN)
            [ SyncFreq, SyncTime ] = CyclicCorr( RxAWGN(iref2:jref2),TxRefOVS );
            SyncNorm = SyncTime/Nref;
            SyncNormABS = abs(SyncTime)/Nref;
             [Peak2, Index2] = max(SyncNormABS);
            RefNorm = LocalReplica(1:NyquistRate:end);
            SyncNormDec = SyncNorm(1:NyquistRate:end);
%             subplot(2,1,2)
%             plot(LocalReplica,'r');  hold on; plot(SyncNormABS);
%             title(['Symbol Delay: ', num2str(SymbolDelay)])
%             xlabel('Index');
%             ylabel('Magnitude');
%             legend('Norm Local Replica','Norm Channel Estimate')
%                 figure  
%                 plot(abs(SyncNormDec)); hold on; plot(RefNorm,'r');
%             AlignmentTerm = 1/SyncNormDec(1);

            % Tx Original and Tx phas and Gain Aligned
            % Start and End Position 
            Pstart = (iref2-(NyquistRate*Ncp));
            Pend = jref2;

            Refstart = iref2;
            Refend = jref2;
            RxSymStart = Pstart-NyquistRate;
            RxSymEnd = Pstart - NyquistRate*Nsymb/2;

%                 RxSymbol_Unaligned = [RxAWGN(RxSymEnd:NyquistRate:RxSymStart) RxAWGN(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];
%                 figure
%                 subplot(2,1,1) 
%                 plot(RxSymbol_Unaligned,'o') 
%                 legend('Rx Signal')
%             Rx_Aligned = RxAWGN*AlignmentTerm;   
               Rx_Aligned = RxAWGN;%(conj(SyncNormDec(1))*RxAWGN)/(SyncNormDec(1)*conj(SyncNormDec(1)));
               h2 =SyncNormDec(1);

              % Recieved Reference Sig 
              RxRef = Rx_Aligned(Refstart:NyquistRate:Refend);
            if Pend+NyquistRate*Nsymb/2<=Nrx

                  RxSymbol2 = [Rx_Aligned(RxSymEnd:NyquistRate:RxSymStart) Rx_Aligned(Pend+1:NyquistRate:Pend+NyquistRate*Nsymb/2)];

                  Noise = ReferenceSignal - RxRef;
                  VarNoise = var(Noise);
                  NoiseArray(2,:) = RxRef;
                  NvPacketEstimate(ii) = VarNoise;
                   
%                   subplot(2,1,2)
%                   plot(RxSymbol,'o') 
%                   legend('Rx Signal Aligned') 
%                   axis([-1 1 -1 1])
            else
                RxSymbol2 = RxAWGN(1:Nsymb);
 
            end
        else
            RxSymbol2 = RxAWGN(1:Nsymb);
        end 
            Antenna2 = RxSymbol2; 
           
            y(1,:) = Antenna1;
            y(2,:) = Antenna2;
            h = [];
            h(1,:) = h1;
            h(2,:) = h2;
            h =repmat(h,1,length(y));
            temp(1,:) =conj(h1)*y(1,:);
            temp(2,:) =conj(h2)*y(2,:);
            sumtemp= sum(temp,1);
            temp2(1,:) = h1*conj(h1);
            temp2(2,:) = h2*conj(h2);
            sumtemp2= sum(temp2,1);
            yHat =  sumtemp/sumtemp2; 
            
            TempN(1,:) =conj(h1)*NoiseArray(1,:);
            TempN(2,:) =conj(h2)*NoiseArray(2,:);
            sumtempN= sum(TempN,1);
            nHat =  sumtempN/sumtemp2; 
            NoTest = var(ReferenceSignal-nHat);
           NvPacketEstimate(ii) = NoTest;
            [ LLRData ] = LLR( LLRBitArray, yHat, k, No );

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
            if PacketErrorStop >50
                break
            end
        end
        NoiseVarArray(i) = No;
        CRCERTotal(1,i)= sum(totalCRCErrors);
        CRCER(1,i) = sum(totalCRCErrors)/PacketCount;
        BitError(1,i) = sum(BitErrorPacket)/PacketCount;
        NewBitError(1,i) = sum(NewBitErrorPacket)/PacketCount;
        SymbolError(1,i) = sum(SymbolErrorPacket)/PacketCount;
        NoiseVarArrayEstimate(i) = sum(NvPacketEstimate)/PacketCount;
%         NoiseVarArray(i) = sum(NvPacket)/PacketCount;
     end
        % Determine Theoretical Block Error Rates
        % ---------------------------------------
        N = TxPacketLength;
        BLER_theo(:,1)  = (1-(1-bit_err_theo).^N);
        MBLERTheo(:,m)= BLER_theo(:,1);
end
bit_err_theo = 0.5.*(1-sqrt((EbNo_lin)./(EbNo_lin+1)));

theoryBer = 0.5.*(1-sqrt(EbNo_lin./(EbNo_lin+1)));
N = TxPacketLength;
BLER_theo(:,1)  = (1-(1-theoryBer).^N);
theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbNo_lin).^(-0.5)); 
p = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
theoryBer_nRx2 = p.^2.*(1+2*(1-p));
% bit_err_theo = 3/8 * ( 1 - sqrt(2/5*EbNo_lin./(1+2/5*EbNo_lin)) ) ...
% + 1/4 * ( 1 - sqrt(18/5*EbNo_lin./(1+18/5*EbNo_lin)) ) ...
% - 1/8 * ( 1 - sqrt(10*EbNo_lin./(1+10*EbNo_lin)) );
% theoryBer = bit_err_theo;
% 
% p  = 1/2 - bit_err_theo;
% p  = 1/2 - 1/2*(1+1./EbNo_lin).^(-1/2);
% bit_err_theo = p.^2.*(1+2*(1-p)); 
figure
semilogy(EbNo_dBArray,bit_err_theo,'m');
hold on
semilogy(EbNo_dBArray,theoryBer_nRx2,'b');
semilogy(EbNo_dBArray,NewBitError,'r:o');
semilogy(EbNo_dBArray,BitError,'k');
% semilogy(EbNo_dB,BLER_theo,':b');
semilogy(EbNo_dB,CRCER,'rp-');
axis([min(EbNo_dB) max(EbNo_dB) 10^-6 10^-0 ])
grid on
xlabel('Eb/No, dB', 'FontSize',14);
ylabel('Bit Error Rate', 'FontSize',14);
set(gca, 'FontSize',14)
title('The effects of Diversity in a Rayleigh Fading Channel', 'FontSize',14);
legend('Thereory BIT Error Rayleigh Fade', 'Theoretical Diveristy Bit Error R = 2', 'Simulated Diversity Bit Error R = 2','Simulated Diversity BLER R = 2', 'FontSize',14);
% legend('Thereory BIT Error AWGN', 'Bit Error', 'Theory Bit Error Rayleigh
% Fade', 'BLER Theory') 

% figure
% loglog(NoiseVarArray,NoiseVarArrayEstimate)
% hold on
% loglog(NoiseVarArray,NoiseVarArray,'k')
% legend('Noise Var vs Estimate Noise Var','Noise Var')
% 
% grid on
figure
subplot(2,1,1)
semilogy(EbNo_dBArray,NoiseVarArrayEstimate)
hold on
semilogy(EbNo_dBArray,NoiseVarArray,'k')
legend('Noise Var Estimate','Noise Var')
title('Noise Variance Estimate for LLR', 'FontSize',14)
% axis([min(EbNo_dBArray) max(EbNo_dB) min(NoiseVarArrayEstimate)
% max(NoiseVarArrayEstimate) ])
% axis([min(NoiseVarArray) max(NoiseVarArray) min(NoiseVarArrayEstimate)
% max(NoiseVarArrayEstimate) ])
xlabel('EbNo dB', 'FontSize',14);
ylabel('Noise Variance', 'FontSize',14);
set(gca, 'FontSize',14); 
grid on
% flipdim(h,1));

Error = abs(NoiseVarArray -NoiseVarArrayEstimate);
subplot(2,1,2)
semilogy(EbNo_dBArray,Error)
title('Noise Variance Error', 'FontSize',14)
xlabel('EbNo dB', 'FontSize',14);
ylabel('Error Percentage', 'FontSize',14);
grid on