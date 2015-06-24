function [ NewPacketSize, TrCHSize,  ParitySize ] = GetPacketSize( Coding, PacketSize,TrCHSize,  Parity, Memory )
% GetPacketSize
%
% Calculates Size of packet required to obtain PacketSize after the CRC has
% been added
%
% Usage :
%               [ NewPacketSize, ParitySize ] = GetPacketSize( PacketSize,
%               Parity )
%
% Where         PacketSize = is the required packet size
%               Parity = is the type of CRC Polynomial used

% Parity Size
% -----------
switch Parity
    case 'gCRC24A'
        ParitySize = 24;
    case 'gCRC24B'
        ParitySize = 24;
    case 'gCRC16'
        ParitySize = 16;
    case 'gCRC8'
        ParitySize = 8;
    otherwise
        ParitySize = 0;
end

if PacketSize<=ParitySize+Memory
    NewPacketSize = PacketSize+Memory;
    if strcmp(Coding,'Uncoded')                 % Check Test
        TrCHSize = PacketSize+Memory+ParitySize;
    end
else
    NewPacketSize = PacketSize - (ParitySize);  % Calculates the required packet size
    if strcmp(Coding,'Uncoded')                 % Check Test
        TrCHSize = PacketSize;
    end
end

end

