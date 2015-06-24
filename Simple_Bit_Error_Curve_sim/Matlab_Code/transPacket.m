function [ BER, BLER, CRCER ] = transPacket( N, k, packetSize, Eb, EbNo_lin, symbArray, parity )
numOfframes = N/packetSize;
BER = zeros(1,length(EbNo_lin));
BLER = zeros(1,length(EbNo_lin));
CRCER = zeros(1,length(EbNo_lin));

for i = 1:length(EbNo_lin)
    frameErrors = zeros(1,numOfframes);
    totalCRCErrors = zeros(1,numOfframes);
    bitErrors = zeros(1,numOfframes);
    for ii = 1:numOfframes
        % Random Data (binary to decimal)
        % -------------------------------
        bitStreamSize = k*packetSize;
        txPacket = (rand(1,bitStreamSize) > 0.5);  
        [ txFrame, newFrameSize ] = crc_TS36212('Transmitter', txPacket, parity);
        
        % Map binary Bits
        % ---------------
        txFrame = reshape(txFrame,k,[]);
        txIdx = Bi2Dec(txFrame,k);
        
        % Modulation
        % ------------
        % From look up table determine symbol to be transmitted
        txSymb = symbArray(txIdx+1);     
        
        % Define AWGN normalised
        % ----------------------
        AWGN = ((randn(1,(newFrameSize/k)))+1i*(randn(1,(newFrameSize/k))))/sqrt(2);
        Noise = sqrt(Eb./(EbNo_lin(i))).*AWGN;
        
        % Transmitt over Medium
        % ---------------------
        rxSymb = txSymb + Noise;
        
        % Demodulation
        % ------------
        % Euclidean distance between recieved signal and symbol array
        % Minimum Euclidean distance and associated index
        rxFrame = false(k, (newFrameSize/k));
        for j = 1:(newFrameSize/k)
            euclideanDistanceAll = abs(symbArray-rxSymb(j)).^2;
            [minED,edIdx] = min(euclideanDistanceAll);
            rxFrame(:,j) = Dec2Bi(edIdx-1, k);
        end
        [ crcError, newFrameSize ] = crc_TS36212('Reciever', rxFrame, parity);
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

