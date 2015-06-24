function [ TxPacket ] = AddBitErrorSizem( TxPacket, m )
% AddBitErrorSizem
%
% This function Adds Transmitted symbol errors
%
% Usage :
%
% [ RxSymbol ] = AddBurstError( RxSymbol, Esnorm, No, ErrorBurstSymbolSize,
% NumberOfBurstError, TrCHFrameSize )
%
% Where         TxPacket       = Transmitted Packet
%
%				m              = Number of errors to be introduced


if m == 0
    return
else
    LengthTxPac = length(TxPacket);

    lenUpdate = m;
    if m>1
        unqArray = randi([1 LengthTxPac],[1 lenUpdate]);
        BitToChange = unique(unqArray);
        while lenUpdate
            BitToChange = unique([BitToChange,unique(unqArray)]);
            lenUpdate = m -length(BitToChange);
            unqArray = randi([1 LengthTxPac],[1 lenUpdate]);
            
        end
    else
            BitToChange = randi([1 LengthTxPac],[1 m]);
    end
    
    for i = 1:m
        if TxPacket(BitToChange(i))==0
            TxPacket(BitToChange(i)) = 1;
        else
            TxPacket(BitToChange(i))= 0;
        end
    end
end

end

