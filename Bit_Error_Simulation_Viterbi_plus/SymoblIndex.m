function [ TxSymbolIdx, newFrameSize ] = SymoblIndex( Packet, k )
% SymoblIndex
%
% This function calculates the symbols to be transmited i.e ((0:k-1)+1)
% where K is the number of bit in the M-ary.  This gives the Indexing
% number for the symbol array not the deciamal equivent, that is calculated
% in function Bi2Dec
%
% Usage :
%
% [ TxSymbolIdx, newFrameSize ] = SymoblIndex( Packet, k )
%
% Where         Packet      = Is the Packet to be transmitted
%
%				k           = Number of Bits for symbol definition

TxBinarySymbols = reshape(Packet,1,k,[]);       % Configure packet to represent each symbol i.e for k = 2 00; 01; 10 etc
[m, n, newFrameSize] = size(TxBinarySymbols); 
TxSymbolIdx = Bi2Dec(TxBinarySymbols,k)+1;      % Calculate Index for symbol array 
end

