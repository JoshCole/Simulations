close all
clear all
M = 4;
PacketSize = 10^6;
% Configure Test
% --------------
EsNo_dB(:,1) = linspace(50, 0, 10); % Noise density dB
EsNo_lin = 10.^(EsNo_dB / 10);      % Noise density linear   
NumEsNo = length(EsNo_lin);         % Numb of Iterations
[ k, Es, Esnorm, Eb, Ebnorm, SymbolArray ] = GetSymbolArrayData( M );
EsNo_lin(1) = inf;
for i = 1:NumEsNo
    [ TxPacket ] = double(rand(1,PacketSize) > 0.5);
    
    [ TxSymbolIdx, TrCHFrameSize ] = SymbolIndex( TxPacket, k );

    [ TxSymbol ] = SymbolArray(TxSymbolIdx);    % Gets Symbols to be Transmitted for Symbol Array

    [ RxSymbol, Noise , No, NoiseVar ] = TransmitSymbol( TxSymbol, Es, Esnorm, EsNo_lin(i) );
    figure
    plot(RxSymbol,'o')
    axis([-1 1 -1 1])
end