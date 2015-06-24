function [ BitProbArray, Logic, BitArraySize ] = LLRSymbolArray( SymbolArray, M, k )
% LLRSymbolArray
%
% Creates bit probablity array, for each bit b0,b1,b2.......bn, determines
% the symbols related for a logic 0 or 1 to occur, hence BitProbArray =
% zeros(k, BitArraySize, Logic)
%
% Usage :
%               [ BitProbArray, Logic, BitArraySize ] = LLRSymbolArray(
%               SymbolArray, M, k )
%
% Where         SymbolArray = Normalised Symbols
%               M = Number of Symbols in M-ary system    
%               k = Number of Bits in Symbol

if M ==1
    M=4;
end
Logic = 2;                      % Number of Logic States [0,1]
BitArraySize = M/2;             % Number of symbols for each logic state for each bit
BitProbArray = zeros(k, BitArraySize, Logic);
Index = fliplr(2.^( 0:k));      % Index for splitting Symbol array, Calculates bn.....b2,b1,b0
BitIndex1 = k;
for i = 1:k                     % Loops for Number of bits
    bitArrayTest = reshape(SymbolArray,[],Index(i));        % Breaks Symbol array into relevent sgments i.e for b0 the array is split into two
    ValueIndex =  1;
    BitLogicZero= cell(1,Index(i)/2);
    BitLogicOne= cell(1,Index(i)/2);
    for ii = 1:2:Index(i)  
        BitLogicZero{ValueIndex} = bitArrayTest(:,ii);      % bn.....b2,b1,b0 for Logic 0
        BitLogicOne{ValueIndex} = bitArrayTest(:,ii+1);     % bn.....b2,b1,b0 for Logic 1
        ValueIndex = ValueIndex+1;
    end
    BitIndex2 = 1;
    for j = 1:Index(i)/2                                    % Combines information order (Every other) for BitProbArray
        for jj = 1:length(BitLogicZero{j})
            BitProbArray(BitIndex1,BitIndex2,1) = BitLogicZero{j}(jj);
            BitProbArray(BitIndex1,BitIndex2,2) = BitLogicOne{j}(jj);
            BitIndex2 = BitIndex2+1;
        end
    end
    BitIndex1 = BitIndex1 - 1;

end


end

