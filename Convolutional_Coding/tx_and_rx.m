function [ BER, BLER, CRCER ] = tx_and_rx( N, k, packetSize, Eb, EbNo_lin, symbArray, parity, genPoly, memory )
numOfframes = N/packetSize;
BER = zeros(1,length(EbNo_lin));
BLER = zeros(1,length(EbNo_lin));
CRCER = zeros(1,length(EbNo_lin));
% Generate XOR Gates
% ------------------
[ xorIndex, numGenPoly, genPolyElements ] = xorRegister( genPoly );

% Define Trellis map for Decoder
% ------------------------------

states = 2^memory;              % Number of states in trellis, dependent on encoding rate i.e 1/3
statesDec = 0:states-1;
[ statesBinary ] = viterbiDec2Bi(statesDec, memory);
zeroPadding = zeros(1, genPolyElements-1-numGenPoly);
currentStates = statesDec + 1;

for i = 1:length(EbNo_lin)
    frameErrors = zeros(1,numOfframes);
    totalCRCErrors = zeros(1,numOfframes);
    bitErrors = zeros(1,numOfframes);
    for ii = 1:numOfframes
        % Random Data (binary to decimal)
        % -------------------------------
        bitStreamSize = k*packetSize;
        txPacket = (rand(1,bitStreamSize) > 0.5);  
        [ dataStreamcrc, newFrameSize ] = crc_TS36212('Transmitter', txPacket, parity);
        
        % Encode data stream
        % ------------------
        register = zeros(1,genPolyElements);
        [ codeWord, codeWordBlock, nBits, State ] = encodeDataStream( register, dataStreamcrc, xorIndex, numGenPoly, memory, 'Encode Data' );

        % Map binary Bits
        % ---------------
        txFrame = codeWordBlock;
        txFrame2 = reshape(codeWord,k,[]);
        [m, newFrameSize] = size(txFrame2);
        txIdx(:,1) = Bi2Dec(txFrame2,k);
        
        % Modulation
        % ------------
        % From look up table determine symbol to be transmitted
        txSymb = symbArray(txIdx+1);     
        
        % Define AWGN normalised
        % ----------------------
        AWGN = ((randn(1,(newFrameSize)))+1i*(randn(1,(newFrameSize))))/sqrt(2);
        Noise = sqrt(Eb./(EbNo_lin(i))).*AWGN;
        
        % Transmitt over Medium
        % ---------------------
        rxSymb = txSymb;% + Noise;
        
    % Demodulation
    % ------------
    % Euclidean distance between recieved signal and symbol array
    % Minimum Euclidean distance and associated index
    rxFrame = false( k, newFrameSize );
    for j = 1:(newFrameSize)
        euclideanDistanceAll = abs(symbArray-rxSymb(j)).^2;
        [minED,edIdx] = min(euclideanDistanceAll);
        rxFrame(:,j) = Dec2Bi(edIdx-1, k);
    end
    rxFrame2 = reshape(rxFrame,1,[]);
        % Compute Branch Logic
        % --------------------
        stages  = nBits+1;          % Number of stages in trellis
        inputValue = [0 1];
        numInput = length(inputValue);
        branchValues = zeros(states,numGenPoly,numInput );
        nextState = zeros(states,memory,numInput );
        for k = 1:states
            for kk = 1:numInput
                register = [inputValue(kk), statesBinary(k,:), zeroPadding];
                [ branchLogic, branchLogicBlock, nBits, State ] = encodeDataStream( register, rxFrame2, xorIndex, numGenPoly, memory, 'Trellis Map' );
                branchValues(k,:,kk) = branchLogic;
                nextState(k,:,kk) = State;
            end
        end
        [ nextState ] = viterbiBi2Dec(nextState);
        
        % Define Nodes
        % ------------
        PreviousState = NaN(states,1);
        Value = inf(states,1);
        Logic = NaN(states,1);
        Surviver = NaN;
        Decode = NaN;
        Visit = zeros(states,1);
        hammDistance = zeros(states,numInput,stages-1 );
        strucStages = struct('States',{currentStates},'BranchLogic',{branchValues},'NextState',{nextState},'PreviousState',{PreviousState},'Visit',{Visit},'Value',{Value},'CodeWord',{[NaN NaN]},'SurviverPath',{Surviver},'Decode',{Logic},'SurviverDecode',{Decode},'HammDistance',{hammDistance});
        strucStages = repmat(strucStages,1,stages);
        index = 1;
        for k = 1:stages-1
            currentCodeWord = repmat(rxFrame2(index:index+numGenPoly-1),states,1);
            for kk = 1:numInput
                hammDistance(:,kk,k) = sum(branchValues(:,:,kk)~=currentCodeWord,numInput);
            end
            index = index+numGenPoly;
            strucStages(k).HammDistance = hammDistance(:,:,k);
            strucStages(k).CodeWord = codeWordBlock(k,:);
        end
        % Initialise Node
        % --------------
        strucStages(1).Value(1) = 0;
        strucStages(1).Visit(1) = 1;
        for k = 1:stages-1
            for kk = 1:states
                if strucStages(k).Value(kk) ~= inf
                    for j = 1:numInput
                        newState = strucStages(i).NextState(kk,j);
                        strucStages(k+1).Visit(newState) = strucStages(i+1).Visit(newState) + 1;
                        if strucStages(k+1).Visit(newState)~= 2
                            strucStages(k+1).Value(newState) = strucStages(k).Value(kk) + strucStages(k).HammDistance(kk,j);
                            strucStages(k+1).PreviousState(newState) = strucStages(k).States(kk);
                            strucStages(k+1).Decode(newState) = j-1;
                        else
                            accMetric = strucStages(k).Value(kk) + strucStages(k).HammDistance(kk,j);
                            if accMetric <= strucStages(k+1).Value(newState)
                                strucStages(k+1).Value(newState) = accMetric;
                                strucStages(k+1).PreviousState(newState) = strucStages(k).States(kk);
                                strucStages(k+1).Decode(newState) = j-1;
                            end
                            
                        end
                    end
                end
            end
        end
        decodedMsg = NaN(1, stages-1);
        for k = stages:-1:2
            [minValue, stateIndex] = min(strucStages(k).Value);
            strucStages(k-1).SurviverPath = strucStages(k).PreviousState(stateIndex);
            strucStages(k-1).SurviverDecode = strucStages(k).Decode(stateIndex);
            decodedMsg(k-1) = strucStages(k-1).SurviverDecode;
        end
%         convDecode(rxFrame, length(genPoly));
        [ crcError, newFrameSize ] = crc_TS36212('Reciever', decodedMsg, parity);
        totalCRCErrors(ii) = crcError;
        if sum(sum(txFrame ~= rxFrame))>= 1;
            bitErrors(ii) = sum(sum(txFrame ~= rxFrame));
            frameErrors(ii) = 1;
        end
    end
    BER(1,i) = sum(bitErrors)/(numOfframes*newFrameSize);
    BLER(1,i) = sum(frameErrors)/numOfframes;
    CRCER(1,i) = sum(totalCRCErrors)/numOfframes;
end
end

