clear all
close all
parity = 'gCRC8';
groupBits = 4;
interleaver = 0;
blSize = 4; % Block noise size

genPoly1 = [1 0 1 1 0 1 1 1 1];
genPoly2 = [1 1 0 1 1 0 0 1 1];
genPoly3 = [1 1 1 1 0 1 0 1 1];

genPoly = [genPoly1;genPoly2;genPoly3];

mem = length(genPoly)-1;
ttiIndex = 4;

EbNo_dB = linspace(0, 10, 11);      % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin = 10.^(EbNo_dB / 10);      % Eb/No values in linear scale
% Map binary Bits
% ---------------
N = 8e2;                            % Number of symbols in simulation
packetSize = 8;
modulationScheme = '16-QAM';
switch modulationScheme
    case 'BPSK'
        L = 1;                              % Distance between bits
        M = 2^L;
        Es = 1;                             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-1,1];
    case 'QPSK'
        L = 2;                              % Distance between bits
        M = 2^L;
        Es = 2;                             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-1-1i,1-1i,-1+1i,1+1i];
    case '16-QAM'
        A= 1;
        L = 4;                              % Distance between bits
        M = 2^L;
        Es = ((2*(M-1))/3)*A^2;             % Energy per symbol
        % Symbol definition
        % ------------
        symbArray =[-3+3i,-3+1i,-3-3i,-3-1i,-1+3i,-1+1i,-1-3i,-1-1i,3+3i,3+1i,3-3i,3-1i,1+3i,1+1i,1-3i,1-1i];
end
k = log2(M);                        % bits per symbol       
Eb = Es / k;                        % Energy per bit
states = 2^mem;                     % Number of states in trellis, dependent on encoding rate i.e 1/3
statesDec = 0:states-1;
[ statesBinary ] = viterbiDec2Bi(statesDec, mem);
numOfframes = N/packetSize;
BER = zeros(1,length(EbNo_lin));
BLER = zeros(1,length(EbNo_lin));
CRCER = zeros(1,length(EbNo_lin));
% Generate XOR Gates
% ------------------
[ xorIndex, numGenPoly, genPolyElements ] = xorRegister( genPoly );
for numEbNo = 1:length(EbNo_lin)
    frameErrors = zeros(1,numOfframes);
    totalCRCErrors = zeros(1,numOfframes);
    bitErrors = zeros(1,numOfframes);
    numFrames = 1;
    numErrors = 0;
    for numFrames = 1:numOfframes
        
        bitStreamSize = k*packetSize;
        txPacket = (rand(1,bitStreamSize) > 0.5);  
        
        % Encode data stream
        % ------------------
        [ dataStreamcrc ] = crc_TS36212('Transmitter', txPacket, parity);
        nBits = length(dataStreamcrc);
        [ codeWord ] = encodeData( dataStreamcrc, genPoly, numGenPoly );
%         register = zeros(1,genPolyElements);
%         [ codeWord, codeWordBlock, nBits, State ] = encodeDataStream( register, dataStreamcrc, xorIndex, numGenPoly, mem, 'Encode Data' );
        if interleaver==1
            % Interleaver
            % -----------
            [ codeWord ] = interLeaver1st( codeWord,ttiIndex );
        end
        
        [frameRow,frameCol]=size(codeWord);
        txFrame = reshape(codeWord,1,k,[]);
        [m, n, newFrameSize] = size(txFrame);
        txIdx = Bi2Dec(txFrame,k)+1;

        % Modulation
        % ------------
        % From look up table determine symbol to be transmitted
        txSymb = symbArray(txIdx);   

        % Define AWGN normalised
        % ----------------------
        AWGN = ((randn(1,newFrameSize))+1i*(randn(1,newFrameSize)))/sqrt(2);
        noiseVar = Eb./(EbNo_lin(numEbNo));
        Noise = sqrt(noiseVar).*AWGN;

        % Transmitt over Medium
        % ---------------------
        rxSymb = txSymb + Noise;
        
%         % Add Block Noise
%         % ---------------
%         midPoint = (newFrameSize/2);
%         rxSymb(1,midPoint:midPoint+blSize-1)=(randn(1,blSize)+1i*randn(1,blSize))/sqrt(2);
        
        % Demodulation
        % ------------
        % Euclidean distance between recieved signal and symbol array
        % Minimum Euclidean distance and associated index
        
%         lenSymbArray = length(symbArray);
%         euclideanDistanceAll = NaN(1, lenSymbArray,newFrameSize);
%         for j = 1:lenSymbArray
%             euclideanDistanceAll(:,j,:) = abs(symbArray(j)-rxSymb).^2;
%         end
%         
%         [minED,rxIdx] = min(euclideanDistanceAll);

        r = rxSymb;     % Recieved Vector
        prob_log = zeros(length(r),k,2);

        yb0(1,:) = [-3+3i,-3+1i,-3-3i,-3-1i,-1+3i,-1+1i,-1-3i,-1-1i];
        yb0(2,:) = [3+3i,3+1i,3-3i,3-1i,1+3i,1+1i,1-3i,1-1i];
        yb1(1,:) = [-3+3i,-3+1i,-3-3i,-3-1i,3+3i,3+1i,3-3i,3-1i];
        yb1(2,:) = [-1+3i,-1+1i,-1-3i,-1-1i,1+3i,1+1i,1-3i,1-1i];
        yb2(1,:) = [-3+3i,-3+1i,-1+3i,-1+1i,3+3i,3+1i,1+3i,1+1i];
        yb2(2,:) = [-3-3i,-3-1i,-1-3i,-1-1i,3-3i,3-1i,1-3i,1-1i];
        yb3(1,:) = [-3+3i,-3-3i,-1+3i,-1-3i,3+3i,3-3i,1+3i,1-3i];
        yb3(2,:) = [-3+1i,-3-1i,-1+1i,-1-1i,3+1i,3-1i,1+1i,1-1i];

        for i = 1:length(r)
            for ii = 1:k/2 % for number of real or imaginary bits

                yprob_b0(i,:,ii) = ((1/(sqrt(2*pi*noiseVar)))*(exp(-(1/2)*((abs(r(i)-yb0(ii,:)).^2)/(noiseVar^2)))));
                yprob_b1(i,:,ii) = ((1/(sqrt(2*pi*noiseVar)))*(exp(-(1/2)*((abs(r(i)-yb1(ii,:)).^2)/(noiseVar^2)))));
                yprob_b2(i,:,ii) = ((1/(sqrt(2*pi*noiseVar)))*(exp(-(1/2)*((abs(r(i)-yb2(ii,:)).^2)/(noiseVar^2)))));
                yprob_b3(i,:,ii) = ((1/(sqrt(2*pi*noiseVar)))*(exp(-(1/2)*((abs(r(i)-yb3(ii,:)).^2)/(noiseVar^2)))));

            end
            % p(r|y=1)/p(r|y=0)
            % ----------------
            yb0_sum(i,:) = sum(yprob_b0(i,:,2))/sum(yprob_b0(i,:,1));
            yb1_sum(i,:) = sum(yprob_b1(i,:,2))/sum(yprob_b1(i,:,1));
            yb2_sum(i,:) = sum(yprob_b2(i,:,2))/sum(yprob_b2(i,:,1));
            yb3_sum(i,:) = sum(yprob_b3(i,:,2))/sum(yprob_b3(i,:,1));

        end

        % [log(abs(b0_log)), r.']
        llr1 = log((yb0_sum));
        llr2 = log((yb1_sum));
        llr3 = log((yb2_sum));
        llr4 = log((yb3_sum));
        % Hard Decision Check
        % -------------------
        outArray = zeros(k,length(rxSymb));
        for i = 1:length(rxSymb)
            if llr1(i)>=0
                outArray(1,i)=1;
            end
            if llr2(i)>=0
                outArray(2,i)=1;
            end
            if llr3(i)>=0
                outArray(3,i)=1;
            end
            if llr4(i)>=0
                outArray(4,i)=1;
            end
        end
        % Define LLR Array
        % ----------------
        cmptc = squeeze(txFrame);
        LLR = zeros(k,length(rxSymb));
        LLR(1,:) = llr1;
        LLR(2,:) = llr2;
        LLR(3,:) = llr3;
        LLR(4,:) = llr4;
        
        rxFrame = reshape(squeeze(outArray),frameRow,frameCol);
        rows = frameCol/numGenPoly;
        for i = 0:rows-1
            cdIndex = i*numGenPoly;
            for ii = 1:numGenPoly
                codeBlock(i+1,ii) = codeWord(cdIndex+ii);
            end
        end
        % Deinterleaving
        % --------------
        if interleaver==1
            newLLR = reshape(squeeze(LLR),frameRow,1,frameCol);
        end
        newLLR = reshape(squeeze(LLR),[],1,numGenPoly);
        
%         rxIdx = squeeze(rxIdx -1);
%         rxFrame = Dec2Bi(rxIdx,k);
%         rxFrame = reshape(squeeze(rxFrame),frameRow,frameCol);
        
        % Define Trellis map for Decoder
        % ------------------------------
        stages  = nBits+1;          % Number of stages in trellis

        % Compute Branch Logic
        % --------------------
        inputValue = [0 1];
        numInput = length(inputValue);
        branchValues = zeros(states,numGenPoly,numInput );
        nextState = zeros(states,mem,numInput );
        for i = 1:states
            for ii = 1:numInput
                register = [inputValue(ii), statesBinary(i,:)];
                [ branchLogic, branchLogicBlock, nBits, State ] = encodeDataStream( register, rxFrame, xorIndex, numGenPoly, mem, 'Trellis Map' );
                branchValues(i,:,ii) = branchLogic;
                nextState(i,:,ii) = State;
            end
        end

        [ nextState ] = viterbiBi2Dec(nextState);
        [ BranchNum ] = viterbiBi2Dec(branchValues);
        currentStates = statesDec + 1;

        % Define Nodes
        % ------------
        PreviousState = NaN(states,1);
        Value = NaN(states,1);
        Logic = NaN(states,1);
        Surviver = NaN;
        Decode = NaN;
        Visit = zeros(states,1);
        hammDistance = zeros(states,numInput,stages-1 );
        strucStages = struct('States',{currentStates},'BranchLogic',{branchValues},'NextState',{nextState},'PreviousState',{PreviousState},'Visit',{Visit},'Value',{Value},'CodeWord',{[NaN NaN]},'SurviverPath',{Surviver},'Decode',{Logic},'SurviverDecode',{Decode},'HammDistance',{hammDistance});
        strucStages = repmat(strucStages,1,stages);

        % LLR Branch Metric
        [m,n]=size(newLLR);
        llrTest= NaN(n,m);
        llrStates=NaN(states,m);
        for i = 1:states
            for ii = 1:numGenPoly
                if statesBinary(i,ii)==0
                   temp = -newLLR(:,:,ii);
                else
                    temp = +newLLR(:,:,ii);
                end
                llrTest(ii,:)=temp;
            end
            llrStates(i,:) = sum(llrTest);
        end
        for i = 1:stages-1
            currentCodeWord = repmat(codeBlock(i,:),states,1);
            for ii = 1:numInput
                hammDistance(:,ii,i) = sum(branchValues(:,:,ii)~=currentCodeWord,numInput);
                lrrM(:,ii,i)=llrStates(BranchNum(:,ii),i);
            end
            strucStages(i).HammDistance = hammDistance(:,:,i);
%             strucStages(i).CodeWord = codeWordBlock(i,:);
            strucStages(i).llrMetric = lrrM(:,:,i);

        end

        % Initialise Node
        % --------------
        strucStages(1).Value(1) = 0;
        strucStages(1).Visit(1) = 1;
        softHard ='Soft';

        for i = 1:stages-1
            for ii = 1:states
                if isnan(strucStages(i).Value(ii)) == 0
                    for j = 1:numInput
                        newState = strucStages(i).NextState(ii,j);
                        strucStages(i+1).Visit(newState) = strucStages(i+1).Visit(newState) + 1;
                        if strucStages(i+1).Visit(newState)~= 2
                            switch softHard
                                case 'Hard'
                                    strucStages(i+1).Value(newState) = strucStages(i).Value(ii) + strucStages(i).HammDistance(ii,j);
                                case 'Soft'
                                    strucStages(i+1).Value(newState) = strucStages(i).Value(ii) + strucStages(i).llrMetric(ii,j);
                            end
                            strucStages(i+1).PreviousState(newState) = strucStages(i).States(ii);
                            strucStages(i+1).Decode(newState) = j-1;
                        else
                            switch softHard
                                case 'Hard'
                                    accMetric = strucStages(i).Value(ii) + strucStages(i).HammDistance(ii,j);
                                    if accMetric <= strucStages(i+1).Value(newState)
                                        strucStages(i+1).Value(newState) = accMetric;
                                        strucStages(i+1).PreviousState(newState) = strucStages(i).States(ii);
                                        strucStages(i+1).Decode(newState) = j-1;
                                    end
                            case 'Soft'
                                accMetric = strucStages(i).Value(ii) + strucStages(i).llrMetric(ii,j);
                                if accMetric >= strucStages(i+1).Value(newState)
                                    strucStages(i+1).Value(newState) = accMetric;
                                    strucStages(i+1).PreviousState(newState) = strucStages(i).States(ii);
                                    strucStages(i+1).Decode(newState) = j-1;
                                end
                            end

                        end
                    end
                end
            end
        end
        decodedMsg = NaN(1, stages-1);
        for i = stages:-1:2
            switch softHard
                case 'Hard'
                    [minValue, stateIndex] = min(strucStages(i).Value);
                case 'Soft'
                    [maxValue, stateIndex] = max(strucStages(i).Value);
            end
            strucStages(i-1).SurviverPath = strucStages(i).PreviousState(stateIndex);
            strucStages(i-1).SurviverDecode = strucStages(i).Decode(stateIndex);
            decodedMsg(i-1) = strucStages(i-1).SurviverDecode;
        end
        [ crcError ] = crc_TS36212('Reciever', decodedMsg, parity);

        totalCRCErrors(numFrames) = crcError;
%         if sum(sum(txFrame ~= rxFrame))>= 1;
%             bitErrors(numFrames) = sum(sum(txFrame ~= rxFrame));
%             frameErrors(numFrames) = 1;
%         end
%         numFrames =numFrames+1;
    end
%     BER(1,numEbNo) = sum(bitErrors)/(numOfframes*newFrameSize);
%     BLER(1,numEbNo) = sum(frameErrors)/numOfframes;
    CRCER(1,numEbNo) = sum(totalCRCErrors)/numOfframes;
end
switch modulationScheme
    case 'BPSK'
        bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));
        Title = title('Block error probability curve for BPSK modulation');
    case 'QPSK'
        bit_err_theo = (1/k)*erfc(sqrt(k*EbNo_lin)/sqrt(2));
        Title = title('Block error probability curve for QPSK Gray Code modulation');
    case '16-QAM'
        bit_err_theo = (1/k)*(3/2)*erfc(sqrt((k/10)*EbNo_lin));
        Title = title('Block error probability curve for 16-QAM Gray Code modulation');
end

% Plot Results
% ---------------
figure(1)
semilogy(EbNo_dB,bit_err_theo,'b.-');
hold on
semilogy(EbNo_dB,CRCER,'ko-');
grid on

legend('Theory Bit Error', 'Packet size 10');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');

% Plot Results
% ---------------
figure(2)
gaindB = bit_err_theo./CRCER;
plot(EbNo_dB,gaindB,'b.-');
grid on

legend('Coding Gain dB');
xlabel('Eb/No, dB');
ylabel('Gain');