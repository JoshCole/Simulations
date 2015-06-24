function [ BlockCode ] = BlockCodeWord( Data, NumberOfGenPoly, txPacketLength )
% BlockCodeWord
%
% Reshapes matrix to form code words for each stage of the viterbi
% decodeder
%
% Usage :
%               [ BlockCode ] = BlockCodeWord( Data, NumberOfGenPoly,
%               txCRCPacketLength )
%
% Where         Data                = Recieved Packet
%               NumberOfGenPoly     = Number of polynomials used in Coding
%               txPacketLength      = Transmitted Packet Length

rows = txPacketLength;
rxBlock = zeros(rows,NumberOfGenPoly);      % Empty array for Block Code
for i = 0:rows-1
    Index = i*NumberOfGenPoly;
    for ii = 1:NumberOfGenPoly              % For each row add 1 to NumberOfGenPoly Bits
        rxBlock(i+1,ii) = Data(Index+ii);
    end
end
BlockCode = reshape(rxBlock,rows,1,NumberOfGenPoly);

end

