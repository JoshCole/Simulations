function [ k, Es, Esnorm, Eb, Ebnorm, SymbArray ] = GetSymbolArrayData( M )
% GetSymbolArrayData
%
% Calculates k, Es, Esnorm, Eb, SymbArray
%
% Usage :
%               [ k, Es, Eb, SymbArray ] = GetSymbolArrayData( M )
%
% Where         M = Number of Symbols in M-ary system             

% Symbol Definition
% -----------------
k = log2(M);                              % Number of bits per symbol
switch M
    case 1
        SymbArray =[1,1i,-1,-1i];
        k = 2;
    case 2
        SymbArray =[-1,1];
    case 4
        SymbArray =[-1-1i,1-1i,-1+1i,1+1i];
    case 16
        SymbArray =[-3+3i,-3+1i,-3-3i,-3-1i,-1+3i,-1+1i,-1-3i,-1-1i,3+3i,3+1i,3-3i,3-1i,1+3i,1+1i,1-3i,1-1i];
end
Esnorm = mean(abs(SymbArray).^2);         % Energy per symbol, the sum(Sr.^2+Si.^2)/(Number of symbols)
Ebnorm = Esnorm/k;
SymbArray = SymbArray./sqrt(Esnorm);      % Normalised symbol array
% Es = 1;                                   % Unit Energy per symbol after
% normalisation
Es = Esnorm;
Eb = Esnorm/k;                                % Energy per bit
end

