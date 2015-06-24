clear all
dataStream = [1 0 1 1 1 0 0 0];
parity = 'gCRC8';

% genPoly1 = [1 1 1];
% genPoly2 = [1 0 1];
% 
% genPoly = [genPoly1;genPoly2];

genPoly1 = [1 0 1 1 0 1 1];
genPoly2 = [1 1 1 1 0 0 1];
genPoly3 = [1 1 1 0 1 0 1];

genPoly = [genPoly1;genPoly2;genPoly3];

% genPoly1 = [1 1 1];
% genPoly2 = [1 0 1];
% genPoly3 = [1 1 1];
% 
% genPoly = [genPoly1;genPoly2;genPoly3];
mem = length(genPoly)-1;
% Generate XOR Gates
% ------------------
[ xorIndex, numGenPoly, genPolyElements ] = xorRegister( genPoly );

% Encode data stream
% ------------------
 [ dataStreamcrc, newFrameSize ] = crc_TS36212('Transmitter', dataStream, parity);
register = zeros(1,genPolyElements);
[ codeWord, codeWordBlock, nBits, State ] = encodeDataStream( register, dataStreamcrc, xorIndex, numGenPoly, mem, 'Encode Data' );

EbNo_dB = 8;      % Energy per bit to noise power spectral density ratio (Eb/No)
EbNo_lin = 10.^(EbNo_dB / 10);      % Eb/No values in linear scale
% Map binary Bits
% ---------------
L = 2;                              % Distance between bits
M = 2^L;
Es = 2;                             % Energy per symbol
k = log2(M);                        % bits per symbol
Eb = Es / k;                        % Energy per bit

txFrame = codeWordBlock;
txFrame2 = reshape(codeWord,k,[]);
[m, newFrameSize] = size(txFrame2);
txIdx(:,1) = Bi2Dec(txFrame2,k);

% Modulation
% ------------
% From look up table determine symbol to be transmitted
symbArray =[-1-1i,1-1i,-1+1i,1+1i];
txSymb = symbArray(txIdx+1);   

% Define AWGN normalised
% ----------------------
AWGN = ((randn(1,(newFrameSize)))+1i*(randn(1,(newFrameSize))))/sqrt(2);
Noise = sqrt(Eb./(EbNo_lin)).*AWGN;

% Transmitt over Medium
% ---------------------
rxSymb = txSymb + Noise;
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
rxFrame3 = reshape(rxFrame2,[],numGenPoly);
% Define Trellis map for Decoder
% ------------------------------
stages  = nBits+1;          % Number of stages in trellis
states = 2^mem;              % Number of states in trellis, dependent on encoding rate i.e 1/3
statesDec = 0:states-1;
[ statesBinary ] = viterbiDec2Bi(statesDec, mem);

% Compute Branch Logic
% --------------------
inputValue = [0 1];
numInput = length(inputValue);
branchValues = zeros(states,numGenPoly,numInput );
nextState = zeros(states,mem,numInput );
for i = 1:states
    for ii = 1:numInput
        register = [inputValue(ii), statesBinary(i,:)];
        [ branchLogic, branchLogicBlock, nBits, State ] = encodeDataStream( register, rxFrame2, xorIndex, numGenPoly, mem, 'Trellis Map' );
        branchValues(i,:,ii) = branchLogic;
        nextState(i,:,ii) = State;
    end
end

[ nextState ] = viterbiBi2Dec(nextState);
currentStates = statesDec + 1;

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
for i = 1:stages-1
    currentCodeWord = repmat(rxFrame2(index:index+numGenPoly-1),states,1);
    for ii = 1:numInput
        hammDistance(:,ii,i) = sum(branchValues(:,:,ii)~=currentCodeWord,numInput);
    end
    index = index+numGenPoly;
    strucStages(i).HammDistance = hammDistance(:,:,i);
    strucStages(i).CodeWord = codeWordBlock(i,:);
end

% Initialise Node
% --------------
strucStages(1).Value(1) = 0;
strucStages(1).Visit(1) = 1;

for i = 1:stages-1
    for ii = 1:states
        if strucStages(i).Value(ii) ~= inf
            for j = 1:numInput
                newState = strucStages(i).NextState(ii,j);
                strucStages(i+1).Visit(newState) = strucStages(i+1).Visit(newState) + 1;
                if strucStages(i+1).Visit(newState)~= 2
                    strucStages(i+1).Value(newState) = strucStages(i).Value(ii) + strucStages(i).HammDistance(ii,j);
                    strucStages(i+1).PreviousState(newState) = strucStages(i).States(ii);
                    strucStages(i+1).Decode(newState) = j-1;
                else
                    accMetric = strucStages(i).Value(ii) + strucStages(i).HammDistance(ii,j);
                    if accMetric <= strucStages(i+1).Value(newState)
                        strucStages(i+1).Value(newState) = accMetric;
                        strucStages(i+1).PreviousState(newState) = strucStages(i).States(ii);
                        strucStages(i+1).Decode(newState) = j-1;
                    end
                    
                end
            end
        end
    end
end
decodedMsg = NaN(1, stages-1);
for i = stages:-1:2
    [minValue, stateIndex] = min(strucStages(i).Value);
    strucStages(i-1).SurviverPath = strucStages(i).PreviousState(stateIndex);
    strucStages(i-1).SurviverDecode = strucStages(i).Decode(stateIndex);
    decodedMsg(i-1) = strucStages(i-1).SurviverDecode;
end
[ crcError, newFrameSize ] = crc_TS36212('Reciever', decodedMsg, parity);