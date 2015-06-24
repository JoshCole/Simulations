clear all
close all
Packet = [1 0 1 1 0 0 1];
[TxPacket, TxPacketLength] = GetCRC('Transmitter', Packet, 'gCRC8');
[ DecodedMessage, crcError ] = GetCRC('Reciever', TxPacket, 'gCRC8');